//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
// 
// Desctiption: MatrixND is a fixed-size matrix class that is suitable for
// stack-allocation.
//
//
//
//===----------------------------------------------------------------------===//
//
// Claudio Perez
//
#pragma once
#include <math.h>
#include <assert.h>
#include <array>

#include "VectorND.h"
#include "Matrix.h"
#include "Vector.h"
#include <utility/Unroll.h>

// #define RESTRICT

#if defined(__GNUC__) || defined(__clang__)
  #define RESTRICT __restrict__
#elif defined(_MSC_VER)
  #define RESTRICT __restrict
#else
  #define RESTRICT
#endif


namespace OpenSees {

template <index_t NR, index_t NC, typename T=double>
struct alignas(64) MatrixND {
  static constexpr int BlasSize = 13*13;
  std::array<T, NR*NC> values;

  //
  // Indexing
  //
  [[gnu::always_inline]] inline constexpr       T& operator()(int i, int j)       noexcept { return values[j*NR + i]; }
  [[gnu::always_inline]] inline constexpr const T& operator()(int i, int j) const noexcept { return values[j*NR + i]; }

  inline constexpr T* data() noexcept { return values.data(); } // &(*this)(0,0); }
  inline constexpr const T* data() const noexcept { return values.data(); }

  // Convert to regular Matrix class
  #ifndef ALLOW_IMPLICIT_MATRIX
  explicit operator Matrix() { return Matrix(&(*this)(0,0), NR, NC);}
  explicit operator const Matrix() const { return Matrix(&(*this)(0,0), NR, NC);}
  #else
  operator Matrix() { return Matrix(&(*this)(0,0), NR, NC);}
  operator const Matrix() const { return Matrix(&(*this)(0,0), NR, NC);}
  #endif
  inline constexpr void zero() noexcept;

  constexpr T
  trace() const noexcept
  {
    static_assert(NR == NC);
    T sum = 0.0;
    for (index_t i = 0; i < NR; ++i) {
      sum += (*this)(i,i);
    }
    return sum;
  }

  const double norm() const noexcept;

  constexpr double determinant() const ;

  inline constexpr MatrixND<NC, NR> transpose() const noexcept;

  inline constexpr MatrixND<NR,NC,T>& addDiagonal(const double vol) noexcept;

  template <class MatT>
    void addMatrix(const MatT& A, const double scale);

  void addTranspose(const MatrixND<NC,NR>& A, const double scale);

  template <int nk> constexpr
    void setMatrixProduct(const MatrixND<NR,nk,T>& A, const MatrixND<nk,NC,T>& B, double scale) noexcept;

  template <int nk> constexpr
    void addMatrixProduct(const MatrixND<NR,nk,T>& A, const MatrixND<nk,NC>& B, double scale) noexcept;

  template <int nk> constexpr
    void addMatrixTransposeProduct(double scale, const MatrixND<nk,NR,T>& B,
                                                 const MatrixND<nk,NC,T>& C) noexcept;

  template <int nk> inline constexpr
    void setMatrixTransposeProduct(const MatrixND<nk, NR, T>& B, const MatrixND<nk, NC, T>& C) noexcept;

  template <int nk> constexpr
    void addMatrixTransposeProduct(const MatrixND<nk, NR, T>& B, const MatrixND<nk, NC, T>& C) noexcept;

  template <class VecA, class VecB>  constexpr MatrixND<NR,NC,T>& 
    addTensorProduct(const VecA& V, const VecB& W, const double scale) noexcept;

  template <class MatT, int nk> void 
    addMatrixProduct(const MatrixND<NR, nk, T> &, const MatT&, double scale);

  template <class MatT, int nk> void 
    addMatrixProduct(double, const MatrixND<NR, nk, T> &, const MatT&, double scale);

  // += A'B
  template <class MatT, int nk> void
    addMatrixTransposeProduct(double thisFact, const MatrixND<nk, NR, T> &, const MatT&, double scale);

  // += A'BA
  template <int nk> int 
    addMatrixTripleProduct(double thisFact, 
                           const MatrixND<nk, NR, T> &, 
                           const MatrixND<nk, nk, T>&, 
                           double scale);

  // += A'BC
  template <int nk, int nl> int 
    addMatrixTripleProduct(double thisFact, 
                           const MatrixND<nk,NR> &A, 
                           const MatrixND<nk,nl> &B, 
                           const MatrixND<nl,NC> &C, double otherFact);

  template <typename F> void map(F func) const;
  template <typename F> void map(F func, MatrixND<NR,NC,T>& destination);

  template<class VecT> inline constexpr MatrixND<NR,NC,T>& addSpin(const VecT& V) noexcept;
  template<class VecT> inline constexpr MatrixND<NR,NC,T>& addSpin(const VecT& V, double scale) noexcept;
  template<class VecT> inline constexpr MatrixND<NR,NC,T>& addSpinSquare(const VecT& V, double scale) noexcept;
  template<class VecT> inline constexpr void addSpinProduct(const VecT& a, const VectorND<NR,T>& b, double scale) noexcept;

  template<class VecT> inline constexpr void 
    addMatrixSpinProduct(const MatrixND<NR,NC,T>& A, const VecT& b, double scale) noexcept;

  template<class MatT> inline constexpr void 
    addSpinMatrixProduct(const VectorND<NR,T>& a, const MatT& B, double scale) noexcept;

  int invert(MatrixND<NR, NC, T> &) const;

  int solve(const VectorND<NR> &V, VectorND<NR> &res) const noexcept;
  template<index_t n>
  int solve(const MatrixND<n, n>& M, MatrixND<n, n>& X) const noexcept;


  template <int nb>
  inline constexpr MatrixND<NR, nb, T>
  bun(const MatrixND<nb,1,T>& b, double scale) const noexcept
  {
    static_assert(NC == 1, "MatrixND::bun: second argument must be a column vector.");

    MatrixND<NR, nb, T> prod;

    for (index_t j = 0; j < nb; ++j)
      for (index_t i = 0; i < NR; ++i)
        prod(i,j) = (*this)(i,0) * b(j,0) * scale;

    return prod;
  }
 

  template <int row0, int row1, int col0, int col1>
  inline MatrixND<row1-row0,col1-col0>
  extract() const noexcept
  {
    MatrixND<row1-row0,col1-col0> m;
    for (int i=0; i<row1-row0; i++)
      for (int j=0; j<col1-col0; j++)
        m(i,j) = (*this)(row0+i, col0+j);
    return m;
  }

  template<int er, int ec>
  inline constexpr MatrixND<er,ec>
  extract(int row0, int col0) const noexcept
  {
    MatrixND<er,ec> m;
    for (int i=0; i<er; i++)
      for (int j=0; j<ec; j++)
        m(i,j) = (*this)(row0+i, col0+j);
    return m;
  }

  template <int init_row, int init_col, int nr, int nc> 
  inline constexpr void
  insert(const MatrixND<nr, nc, double> &M) noexcept;

  template <int init_row, int init_col, int nr, int nc> 
  inline constexpr void
  insert(const MatrixND<nr, nc, double> &M, double fact) noexcept;

  template <int nr, int nc> 
  inline constexpr void
  insert(const MatrixND<nr, nc, double> &M, int init_row, int init_col, double fact) noexcept;

  template <int nr, int nc> 
  inline constexpr void
  assemble(const MatrixND<nr, nc, double> &M, int init_row, int init_col, double fact) noexcept;

  template <int nr> inline void
  assemble(const VectorND<nr> &v, int init_row, int init_col, double fact) noexcept;
  
  template<int nr, int nc> inline void
  assembleTranspose(const MatrixND<nr, nc, double> &M, int init_row, int init_col, double fact) noexcept;

  template<int nr> void
  assembleTranspose(const VectorND<nr> &v, int init_row, int init_col, double scale) noexcept;


  //
  // Operators
  //

  inline constexpr MatrixND &
  operator=(const Matrix &other) noexcept
  {
    for (index_t j = 0; j < NC; ++j) {
      for (index_t i = 0; i < NR; ++i) {
        (*this)(i,j) = other(i,j);
      }
    }
    return *this;
  }

  constexpr MatrixND &
  operator+=(const double value) noexcept {
    for (index_t j = 0; j < NC; ++j) {
      for (index_t i = 0; i < NR; ++i) {
        (*this)(i,j) += value;
      }
    }
    return *this;
  }

  constexpr MatrixND &
  operator+=(const MatrixND &other) noexcept {
    for (index_t j = 0; j < NC; ++j) {
      for (index_t i = 0; i < NR; ++i) {
        (*this)(i,j) += other(i,j);
      }
    }
    return *this;
  }
  
  constexpr MatrixND &
  operator-=(const MatrixND &other) noexcept
  {
    for (index_t j = 0; j < NC; ++j)
      for (index_t i = 0; i < NR; ++i)
        (*this)(i,j) -= other(i,j);

    return *this;
  }

  inline constexpr MatrixND &
  operator*=(T const scalar) noexcept
  {
    for (index_t j = 0; j < NC; ++j)
      for (index_t i = 0; i < NR; ++i)
        (*this)(i,j) *= scalar;

    return *this;
  }

  inline constexpr MatrixND &
  operator/=(T const scalar) noexcept
  {
    for (index_t j = 0; j < NC; ++j)
      for (index_t i = 0; i < NR; ++i)
        (*this)(i,j) /= scalar;

    return *this;
  }

  inline constexpr VectorND<NC>
  operator^(const VectorND<NR> &V) const noexcept
  {
    VectorND<NC> result;

    const double *dataPtr = &(*this)(0,0);
    Repeat<NC>([&](auto i_) {
      constexpr int i = i_.value;
      result[i] = 0.0;
      for (int j=0; j<NR; j++)
        result[i] += *dataPtr++ * V[j];
    });

    return result;
  }

  //
  // Friends 
  //
  friend constexpr MatrixND
  operator+(MatrixND left, const MatrixND &right) noexcept {
    left += right; 
    return left;
  }

  friend constexpr MatrixND
  operator-(MatrixND left, const MatrixND &right) noexcept {
    left -= right; 
    return left;
  }
  
  friend constexpr MatrixND // scalar * Matrix
  operator*(T scalar, MatrixND mat) noexcept {
    mat *= scalar;
    return mat;
  }

  friend constexpr MatrixND // Matrix * scalar
  operator*(MatrixND mat, T scalar) noexcept {
    mat *= scalar;
    return mat;
  }

  template <index_t J>
  inline constexpr friend MatrixND<NR, J>
  operator*(const MatrixND<NR, NC> &left, const MatrixND<NC, J> &right) noexcept {
    MatrixND<NR, J> prod;
    // if constexpr (NR*NC > 13*13) {
    //   prod.setMatrixProduct(left, right, 1);
    // }
    // else
      for (index_t i = 0; i < NR; ++i) {
        for (index_t j = 0; j < J; ++j) {
          prod(i, j) = 0.0;
          for (index_t k = 0; k < NC; ++k) {
            prod(i, j) += left(i,k) * right(k,j);
          }
        }
      }
    return prod;
  }

  template <index_t J>
  friend  MatrixND<NR, J>
  operator*(const MatrixND<NR, NC> &left, const Matrix &right) noexcept {
    MatrixND<NR, J> prod;
    for (index_t i = 0; i < NR; ++i) {
      for (index_t j = 0; j < J; ++j) {
        prod(i, j) = 0.0;
        for (index_t k = 0; k < NC; ++k) {
          prod(i, j) += left(i,k) * right(k,j);
        }
      }
    }
    return prod;
  }

  // Matrix*Vector
  inline constexpr friend  VectorND<NR>
  operator*(const MatrixND<NR, NC> &left, const VectorND<NC> &right) noexcept {
    VectorND<NR> prod;
    for (index_t i = 0; i < NR; ++i) {
      prod[i] = 0.0;
      for (index_t k = 0; k < NC; ++k) {
        prod[i] += left(i,k) * right[k];
      }
    }
    return prod;
  }

  friend  VectorND<NR>
  operator*(const MatrixND<NR, NC> &left, const Vector &right) noexcept {
    VectorND<NR> prod;
    for (index_t i = 0; i < NR; ++i) {
      prod[i] = 0.0;
      for (index_t k = 0; k < NC; ++k) {
        prod[i] += left(i,k) * right(k);
      }
    }
    return prod;
  }

  template <index_t K>
  inline constexpr friend MatrixND<NC,K>
  operator^(const MatrixND<NR, NC> &left, const MatrixND<NR, K> &right) {
    MatrixND<NC, K> prod;
    if constexpr (NR*NC > 16) {
      prod.setMatrixTransposeProduct(left, right);
    }
    else 
    {
      for (index_t i = 0; i < NC; ++i) {
        for (index_t j = 0; j < K; ++j) {
          prod(i, j) = 0.0;
          for (index_t k = 0; k < NR; ++k) {
            prod(i, j) += left(k,i) * right(k,j);
          }
        }
      }
    }
    return prod;
  }

  friend constexpr MatrixND
  operator/(MatrixND mat, T scalar) {
    mat /= scalar;
    return mat;
  }
};

} // namespace OpenSees
#include "MatrixND.tpp"
