
// Description: This file contains the class definition for
// SparseGenCSC. It stores a sparse matrix in the column-compressed
// format required by UMFPACK.
//
#pragma once

#include <LinearSOE.h>
#include <Vector.h>

#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

class UmfpackSolver02;
class Channel;
class FEM_ObjectBroker;
class Graph;
class ID;
class Matrix;

class SparseGenCSC : public LinearSOE
{
public:
    static constexpr std::size_t DefaultScatterCacheSize = 256u * 1024u * 1024u;

    // LinearSOE takes ownership of the supplied solver.
    explicit SparseGenCSC(UmfpackSolver02 &,
                             std::size_t maxScatterBytes = DefaultScatterCacheSize);
    SparseGenCSC();
    ~SparseGenCSC() override;

    SparseGenCSC(const SparseGenCSC &) = delete;
    SparseGenCSC &operator=(const SparseGenCSC &) = delete;
    SparseGenCSC(SparseGenCSC &&) = delete;
    SparseGenCSC &operator=(SparseGenCSC &&) = delete;

    int getNumEqn() const override;
    int setSize(Graph &) override;
    int addA(const Matrix &, const ID &, double fact = 1.0) override;
    int addB(const Vector &, const ID &, double fact = 1.0) override;
    int setB(const Vector &, double fact = 1.0) override;

    void zeroA() override;
    void zeroB() override;

    const Vector &getX() override;
    const Vector &getB() override;
    double normRHS() override;

    void setX(int loc, double value) override;
    void setX(const Vector &x) override;

    //
    //

    const std::vector<int>& getPointers()   const { return Ap; }
    const std::vector<int>& getRowIndices() const { return Ai; }
    const std::vector<double>& getValues()  const { return Ax; }

    // Interface required for FORTRAN solvers like MUMPS
    int*    rawPointers()   { return Ap.data(); }
    int*    rawRowIndices() { return Ai.data(); }
    double* rawValues()  { return Ax.data(); }

    int sendSelf(int commitTag, Channel &) override;
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &) override;

    friend class UmfpackSolver02;

private:
    struct CellCache {
    public:
      CellCache(std::size_t maxScatterBytes)
      :scatterCache(), scatterCacheBytes(0),
       maxScatterCacheBytes(maxScatterBytes)
      {
        scatterCache.max_load_factor(0.5F);
      }

      struct Scatter {
        int size = 0;
        std::vector<int> data;

        bool matches(const ID &) const;
        std::size_t storageSize() const;
      };

      static constexpr int MinimumScatterSize = 3;
      static constexpr std::size_t ScatterCacheEntryOverhead = 8u * sizeof(void *);

      void reset() {
        scatterCache.clear();
        scatterCache.rehash(0);
        scatterCacheBytes = 0;
      }
      bool canCache(std::size_t) const;
      std::uint64_t hash(const ID &) const;
      void buildScatter(const ID &, Scatter &, int size, const std::vector<int>& Ap, const std::vector<int>& Ai);

      std::unordered_multimap<std::uint64_t, Scatter> scatterCache;
      std::size_t scatterCacheBytes;
      std::size_t maxScatterCacheBytes;
    } cache;

    void assembleCached(const Matrix &, const CellCache::Scatter &, double);
    void assembleUncached(const Matrix &, const ID &, double);

    Vector X;
    Vector B;
    std::vector<double> Ax; // Ax[j] is the value of the j-th non-zero.
    // CSC format
    std::vector<int> Ap; // Ap[i] is the index of the first non-zero in column i
    std::vector<int> Ai; // Ai[j] is the row index of the j-th non-zero
    std::vector<int> colA;

};
