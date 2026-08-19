
#include <UmfpackLinSOE02.h>
#include <UmfpackSolver02.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Graph.h>
#include <ID.h>
#include <Matrix.h>
#include <OPS_Globals.h>
#include <Vertex.h>
#include <VertexIter.h>

#include <algorithm>
#include <utility>

UmfpackLinSOE02::UmfpackLinSOE02(UmfpackSolver02 &theSolver,
                                 std::size_t maxScatterBytes)
  : LinearSOE(theSolver, LinSOE_TAGS_UmfpackLinSOE02),
    X(), B(), Ap(), Ai(), Ax(), scatterCache(), scatterCacheBytes(0),
    maxScatterCacheBytes(maxScatterBytes)
{
    scatterCache.max_load_factor(0.5F);
    theSolver.setLinearSOE(*this);
}

UmfpackLinSOE02::UmfpackLinSOE02()
  : LinearSOE(LinSOE_TAGS_UmfpackLinSOE02),
    X(), B(), Ap(), Ai(), Ax(), scatterCache(), scatterCacheBytes(0),
    maxScatterCacheBytes(DefaultScatterCacheSize)
{
    scatterCache.max_load_factor(0.5F);
}

UmfpackLinSOE02::~UmfpackLinSOE02() = default;

int
UmfpackLinSOE02::getNumEqn() const
{
    return X.Size();
}

int
UmfpackLinSOE02::setSize(Graph &theGraph)
{
    const int size = theGraph.getNumVertex();
    if (size < 0) {
        opserr << "UmfpackLinSOE02::setSize() - negative system size\n";
        return -1;
    }

    std::size_t estimatedNnz = 0;
    VertexIter &vertices = theGraph.getVertices();
    Vertex *vertex = nullptr;
    while ((vertex = vertices()) != nullptr) {
        estimatedNnz += static_cast<std::size_t>(vertex->getAdjacency().Size()) + 1;
    }

    Ap.clear();
    Ai.clear();
    Ax.clear();
    scatterCache.clear();
    scatterCache.rehash(0);
    scatterCacheBytes = 0;

    Ap.reserve(static_cast<std::size_t>(size) + 1);
    Ai.reserve(estimatedNnz);
    Ap.push_back(0);

    std::vector<int> rows;
    for (int col = 0; col < size; ++col) {
        vertex = theGraph.getVertexPtr(col);
        if (vertex == nullptr) {
            opserr << "WARNING UmfpackLinSOE02::setSize() - vertex " << col
                   << " is not in graph\n";
            return -1;
        }

        const ID &adjacency = vertex->getAdjacency();
        rows.clear();
        rows.reserve(static_cast<std::size_t>(adjacency.Size()) + 1);
        rows.push_back(vertex->getTag());
        for (int i = 0; i < adjacency.Size(); ++i) {
            rows.push_back(adjacency(i));
        }

        std::sort(rows.begin(), rows.end());
        const std::vector<int>::iterator last = std::unique(rows.begin(), rows.end());
        Ai.insert(Ai.end(), rows.begin(), last);
        Ap.push_back(static_cast<int>(Ai.size()));
    }

    Ax.assign(Ai.size(), 0.0);
    X.resize(size);
    X.Zero();
    B.resize(size);
    B.Zero();

    LinearSOESolver *solver = this->getSolver();
    if (solver == nullptr) {
        return -1;
    }

    return solver->setSize();
}

int
UmfpackLinSOE02::addA(const Matrix &m, const ID &id, double fact)
{
    if (fact == 0.0) {
        return 0;
    }

    const int idSize = id.Size();
    if (idSize != m.noRows() || idSize != m.noCols()) {
        opserr << "UmfpackLinSOE02::addA() - Matrix and ID have incompatible sizes\n";
        return -1;
    }

    if (idSize == 0) {
        return 0;
    }

    if (idSize < MinimumScatterSize) {
        this->assembleUncached(m, id, fact);
        return 0;
    }

    const std::uint64_t idHash = this->hashID(id);
    const auto cachedRange = scatterCache.equal_range(idHash);
    for (auto cached = cachedRange.first; cached != cachedRange.second; ++cached) {
        if (cached->second.matches(id)) {
            this->assembleCached(m, cached->second, fact);
            return 0;
        }
    }

    const std::size_t entries = static_cast<std::size_t>(idSize) * idSize;
    const std::size_t scatterBytes = sizeof(Scatter)
                                   + (entries + static_cast<std::size_t>(idSize)) * sizeof(int)
                                   + ScatterCacheEntryOverhead;
    if (this->canCache(scatterBytes)) {
        Scatter scatter;
        this->buildScatter(id, scatter);
        const auto storedScatter = scatterCache.emplace(idHash, std::move(scatter));
        scatterCacheBytes += storedScatter->second.storageSize() + ScatterCacheEntryOverhead;
        this->assembleCached(m, storedScatter->second, fact);
    } else {
        this->assembleUncached(m, id, fact);
    }

    return 0;
}

int
UmfpackLinSOE02::addB(const Vector &v, const ID &id, double fact)
{
    if (fact == 0.0) {
        return 0;
    }

    const int idSize = id.Size();
    if (idSize != v.Size()) {
        opserr << "UmfpackLinSOE02::addB() - Vector and ID have incompatible sizes\n";
        return -1;
    }

    const int size = B.Size();
    if (fact == 1.0) {
        for (int i = 0; i < idSize; ++i) {
            const int row = id(i);
            if (row >= 0 && row < size) {
                B(row) += v(i);
            }
        }
    } else if (fact == -1.0) {
        for (int i = 0; i < idSize; ++i) {
            const int row = id(i);
            if (row >= 0 && row < size) {
                B(row) -= v(i);
            }
        }
    } else {
        for (int i = 0; i < idSize; ++i) {
            const int row = id(i);
            if (row >= 0 && row < size) {
                B(row) += fact * v(i);
            }
        }
    }

    return 0;
}

int
UmfpackLinSOE02::setB(const Vector &v, double fact)
{
    if (fact == 0.0) {
        B.Zero();
        return 0;
    }

    if (v.Size() != B.Size()) {
        opserr << "UmfpackLinSOE02::setB() - incompatible vector sizes\n";
        return -1;
    }

    const int size = B.Size();
    if (fact == 1.0) {
        for (int i = 0; i < size; ++i) {
            B(i) = v(i);
        }
    } else if (fact == -1.0) {
        for (int i = 0; i < size; ++i) {
            B(i) = -v(i);
        }
    } else {
        for (int i = 0; i < size; ++i) {
            B(i) = fact * v(i);
        }
    }

    return 0;
}

void
UmfpackLinSOE02::zeroA()
{
    std::fill(Ax.begin(), Ax.end(), 0.0);
}

void
UmfpackLinSOE02::zeroB()
{
    B.Zero();
}

const Vector &
UmfpackLinSOE02::getX()
{
    return X;
}

const Vector &
UmfpackLinSOE02::getB()
{
    return B;
}

double
UmfpackLinSOE02::normRHS()
{
    return B.Norm();
}

void
UmfpackLinSOE02::setX(int loc, double value)
{
    if (loc >= 0 && loc < X.Size()) {
        X(loc) = value;
    }
}

void
UmfpackLinSOE02::setX(const Vector &x)
{
    if (x.Size() == X.Size()) {
        X = x;
    }
}

int
UmfpackLinSOE02::sendSelf(int, Channel &)
{
    return 0;
}

int
UmfpackLinSOE02::recvSelf(int, Channel &, FEM_ObjectBroker &)
{
    return 0;
}

bool
UmfpackLinSOE02::Scatter::matches(const ID &id) const
{
    if (size != id.Size()) {
        return false;
    }

    for (int i = 0; i < id.Size(); ++i) {
        if (data[static_cast<std::size_t>(i)] != id(i)) {
            return false;
        }
    }

    return true;
}

std::size_t
UmfpackLinSOE02::Scatter::storageSize() const
{
    return sizeof(Scatter)
         + data.capacity() * sizeof(int);
}

bool
UmfpackLinSOE02::canCache(std::size_t bytes) const
{
    return bytes <= maxScatterCacheBytes
        && scatterCacheBytes <= maxScatterCacheBytes - bytes;
}

std::uint64_t
UmfpackLinSOE02::hashID(const ID &id) const
{
    std::uint64_t hash = 14695981039346656037ULL;
    for (int i = 0; i < id.Size(); ++i) {
        hash ^= static_cast<std::uint32_t>(id(i));
        hash *= 1099511628211ULL;
    }
    hash ^= static_cast<std::uint32_t>(id.Size());
    return hash * 1099511628211ULL;
}

void
UmfpackLinSOE02::buildScatter(const ID &id, Scatter &scatter)
{
    const int idSize = id.Size();
    const std::size_t entries = static_cast<std::size_t>(idSize) * idSize;
    scatter.size = idSize;
    scatter.data.resize(static_cast<std::size_t>(idSize) + entries);
    int *equation = scatter.data.data();
    int *offset = equation + idSize;
    std::fill(offset, scatter.data.data() + scatter.data.size(), -1);

    for (int i = 0; i < idSize; ++i) {
        equation[i] = id(i);
    }

    const int size = X.Size();
    const int *rows = Ai.data();
    for (int j = 0; j < idSize; ++j) {
        const int col = equation[j];
        if (col < 0 || col >= size) {
            continue;
        }

        const int begin = Ap[static_cast<std::size_t>(col)];
        const int end = Ap[static_cast<std::size_t>(col) + 1];
        const int *first = rows + begin;
        const int *last = rows + end;
        int *column = offset + static_cast<std::size_t>(j) * idSize;
        for (int i = 0; i < idSize; ++i) {
            const int row = equation[i];
            if (row < 0 || row >= size) {
                continue;
            }

            const int *entry = std::lower_bound(first, last, row);
            if (entry != last && *entry == row) {
                column[i] = static_cast<int>(entry - rows);
            }
        }
    }
}

void
UmfpackLinSOE02::assembleCached(const Matrix &m, const Scatter &scatter, double fact)
{
    const int idSize = scatter.size;
    const int *offset = scatter.data.data() + idSize;
    double *values = Ax.data();

    if (fact == 1.0) {
        for (int j = 0; j < idSize; ++j) {
            const double *matrixColumn = &m(0, j);
            for (int i = 0; i < idSize; ++i) {
                const int slot = *offset++;
                if (slot >= 0) {
                    values[slot] += matrixColumn[i];
                }
            }
        }
    } else if (fact == -1.0) {
        for (int j = 0; j < idSize; ++j) {
            const double *matrixColumn = &m(0, j);
            for (int i = 0; i < idSize; ++i) {
                const int slot = *offset++;
                if (slot >= 0) {
                    values[slot] -= matrixColumn[i];
                }
            }
        }
    } else {
        for (int j = 0; j < idSize; ++j) {
            const double *matrixColumn = &m(0, j);
            for (int i = 0; i < idSize; ++i) {
                const int slot = *offset++;
                if (slot >= 0) {
                    values[slot] += fact * matrixColumn[i];
                }
            }
        }
    }
}

void
UmfpackLinSOE02::assembleUncached(const Matrix &m, const ID &id, double fact)
{
    const int idSize = id.Size();
    const int size = X.Size();
    const int *rows = Ai.data();
    double *values = Ax.data();

    if (fact == 1.0) {
        for (int j = 0; j < idSize; ++j) {
            const int col = id(j);
            if (col < 0 || col >= size) {
                continue;
            }

            const int *first = rows + Ap[static_cast<std::size_t>(col)];
            const int *last = rows + Ap[static_cast<std::size_t>(col) + 1];
            const double *matrixColumn = &m(0, j);
            for (int i = 0; i < idSize; ++i) {
                const int row = id(i);
                if (row < 0 || row >= size) {
                    continue;
                }

                const int *entry = std::lower_bound(first, last, row);
                if (entry != last && *entry == row) {
                    values[entry - rows] += matrixColumn[i];
                }
            }
        }
    } else if (fact == -1.0) {
        for (int j = 0; j < idSize; ++j) {
            const int col = id(j);
            if (col < 0 || col >= size) {
                continue;
            }

            const int *first = rows + Ap[static_cast<std::size_t>(col)];
            const int *last = rows + Ap[static_cast<std::size_t>(col) + 1];
            const double *matrixColumn = &m(0, j);
            for (int i = 0; i < idSize; ++i) {
                const int row = id(i);
                if (row < 0 || row >= size) {
                    continue;
                }

                const int *entry = std::lower_bound(first, last, row);
                if (entry != last && *entry == row) {
                    values[entry - rows] -= matrixColumn[i];
                }
            }
        }
    } else {
        for (int j = 0; j < idSize; ++j) {
            const int col = id(j);
            if (col < 0 || col >= size) {
                continue;
            }

            const int *first = rows + Ap[static_cast<std::size_t>(col)];
            const int *last = rows + Ap[static_cast<std::size_t>(col) + 1];
            const double *matrixColumn = &m(0, j);
            for (int i = 0; i < idSize; ++i) {
                const int row = id(i);
                if (row < 0 || row >= size) {
                    continue;
                }

                const int *entry = std::lower_bound(first, last, row);
                if (entry != last && *entry == row) {
                    values[entry - rows] += fact * matrixColumn[i];
                }
            }
        }
    }
}
