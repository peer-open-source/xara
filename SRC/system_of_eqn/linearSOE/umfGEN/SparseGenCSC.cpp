
#include <SparseGenCSC.h>
#include <UmfpackSolver02.h>
//#include <UmfpackGenLinSolver.h>
//#define UmfpackSolver02 UmfpackGenLinSolver

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

SparseGenCSC::SparseGenCSC(UmfpackSolver02 &theSolver,
                                 std::size_t maxScatter)
  : LinearSOE(theSolver, -1),
    X(), B(), Ap(), Ai(), Ax(),
    cache(maxScatter)
{
  theSolver.setLinearSOE(*this);
}

SparseGenCSC::SparseGenCSC()
  : LinearSOE(-1),
    X(), B(), Ap(), Ai(), Ax(),
#if 1
    cache(DefaultScatterCacheSize)
#else
    scatterCache(), scatterCacheBytes(0),
    maxScatterCacheBytes(DefaultScatterCacheSize)
#endif
{

}

SparseGenCSC::~SparseGenCSC() = default;

int
SparseGenCSC::getNumEqn() const
{
  return X.Size();
}

int
SparseGenCSC::setSize(Graph &theGraph)
{
  const int size = theGraph.getNumVertex();
  if (size < 0) {
    opserr << "SparseGenCSC::setSize() - negative system size\n";
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
  cache.reset();

  Ap.reserve(static_cast<std::size_t>(size) + 1);
  Ai.reserve(estimatedNnz);
  Ap.push_back(0);

  std::vector<int> rows;
  for (int col = 0; col < size; ++col) {
    vertex = theGraph.getVertexPtr(col);
    if (vertex == nullptr) {
      opserr << "WARNING SparseGenCSC::setSize() - vertex " << col
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
SparseGenCSC::addA(const Matrix &m, const ID &id, double fact)
{
  if (fact == 0.0) {
    return 0;
  }

  const int idSize = id.Size();
  if (idSize != m.noRows() || idSize != m.noCols()) {
    opserr << "SparseGenCSC::addA() - Matrix and ID have incompatible sizes\n";
    return -1;
  }

  if (idSize == 0) {
    return 0;
  }

  if (idSize < CellCache::MinimumScatterSize) {
    this->assembleUncached(m, id, fact);
    return 0;
  }

  const std::uint64_t idHash = cache.hash(id);
  const auto cachedRange = cache.scatterCache.equal_range(idHash);
  for (auto cached = cachedRange.first; cached != cachedRange.second; ++cached) {
    if (cached->second.matches(id)) {
      this->assembleCached(m, cached->second, fact);
      return 0;
    }
  }

  const std::size_t entries = static_cast<std::size_t>(idSize) * idSize;
  const std::size_t scatterBytes = sizeof(CellCache::Scatter)
                                  + (entries + static_cast<std::size_t>(idSize)) * sizeof(int)
                                  + cache.ScatterCacheEntryOverhead;
  if (cache.canCache(scatterBytes)) {
    CellCache::Scatter scatter;
    cache.buildScatter(id, scatter, X.Size(), Ap, Ai);
    const auto storedScatter = cache.scatterCache.emplace(idHash, std::move(scatter));
    cache.scatterCacheBytes += storedScatter->second.storageSize() + cache.ScatterCacheEntryOverhead;
    this->assembleCached(m, storedScatter->second, fact);
  } else {
    this->assembleUncached(m, id, fact);
  }

  return 0;
}

int
SparseGenCSC::addB(const Vector &v, const ID &id, double fact)
{
    if (fact == 0.0) {
        return 0;
    }

    const int idSize = id.Size();
    if (idSize != v.Size()) {
        opserr << "SparseGenCSC::addB() - Vector and ID have incompatible sizes\n";
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
SparseGenCSC::setB(const Vector &v, double fact)
{
    if (fact == 0.0) {
        B.Zero();
        return 0;
    }

    if (v.Size() != B.Size()) {
        opserr << "SparseGenCSC::setB() - incompatible vector sizes\n";
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
SparseGenCSC::zeroA()
{
    std::fill(Ax.begin(), Ax.end(), 0.0);
}

void
SparseGenCSC::zeroB()
{
    B.Zero();
}

const Vector &
SparseGenCSC::getX()
{
    return X;
}

const Vector &
SparseGenCSC::getB()
{
    return B;
}

double
SparseGenCSC::normRHS()
{
    return B.Norm();
}

void
SparseGenCSC::setX(int loc, double value)
{
    if (loc >= 0 && loc < X.Size()) {
        X(loc) = value;
    }
}

void
SparseGenCSC::setX(const Vector &x)
{
    if (x.Size() == X.Size()) {
        X = x;
    }
}

int
SparseGenCSC::sendSelf(int, Channel &)
{
    return 0;
}

int
SparseGenCSC::recvSelf(int, Channel &, FEM_ObjectBroker &)
{
    return 0;
}

bool
SparseGenCSC::CellCache::Scatter::matches(const ID &id) const
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
SparseGenCSC::CellCache::Scatter::storageSize() const
{
    return sizeof(Scatter)
         + data.capacity() * sizeof(int);
}

bool
SparseGenCSC::CellCache::canCache(std::size_t bytes) const
{
    return bytes <= maxScatterCacheBytes
        && scatterCacheBytes <= maxScatterCacheBytes - bytes;
}

std::uint64_t
SparseGenCSC::CellCache::hash(const ID &id) const
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
SparseGenCSC::CellCache::buildScatter(const ID &id, Scatter &scatter, int size, const std::vector<int>& Ap, const std::vector<int>& Ai)
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

//const int size = X.Size();
  const int *rows = Ai.data();
  for (int j = 0; j < idSize; ++j) {
    const int col = equation[j];
    if (col < 0 || col >= size) {
      continue;
    }

    const int begin  = Ap[static_cast<std::size_t>(col)];
    const int end    = Ap[static_cast<std::size_t>(col) + 1];
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
SparseGenCSC::assembleCached(const Matrix &m, const CellCache::Scatter &scatter, double fact)
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
SparseGenCSC::assembleUncached(const Matrix &m, const ID &id, double fact)
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
