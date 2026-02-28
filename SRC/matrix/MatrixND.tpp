//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
//
// Claudio Perez
//
#pragma once
#include "MatrixND.h"

#include "routines/xblas.h"
#include "routines/cmx.h"


//  NOTE: Currently MATRIX_BRANCHING is NECCESSARY to avoid undefined behavior
//  when using functions like addMatrixTripleProduct function on uninitialized
//  matrices.
#define MATRIX_BRANCHING

namespace OpenSees {

template <int nr, int nc, typename T=double>
MatrixND(const T (&)[nc][nr])->MatrixND<nr, nc, T>;


template <index_t nr, index_t nc, typename T>
inline constexpr void
MatrixND<nr, nc, T>::zero() noexcept
{
  values.fill(0.0);
}

template <index_t nr, index_t nc, typename T>
const double
MatrixND<nr, nc, T>::norm() const noexcept
{
  double sum = 0.0;
  for (index_t j = 0; j < nc; ++j) {
    for (index_t i = 0; i < nr; ++i) {
      const double val = static_cast<double>((*this)(i,j));
      sum += val * val;
    }
  }
  return std::sqrt(sum);
}

#if 0
template <index_t nr, index_t nc, typename T>
constexpr double
MatrixND<nr, nc, T>::determinant() const
{
  static_assert(nr == nc, "Matrix must be square");
  static_assert(nr > 1 && nr < 4, "Matrix must be between 2x2 and 3x3");
  if constexpr (nr == 2) {
    return values[0][0] * values[1][1] - values[0][1] * values[1][0];
  }
  if constexpr (nr == 3) {
    return values[0][0] * (values[1][1] * values[2][2] - values[1][2] * values[2][1]) -
           values[0][1] * (values[1][0] * values[2][2] - values[1][2] * values[2][0]) +
           values[0][2] * (values[1][0] * values[2][1] - values[1][1] * values[2][0]);
  }
}
#endif

template <index_t nr, index_t nc, typename T>
constexpr MatrixND<nc, nr>
MatrixND<nr, nc, T>::transpose() const noexcept
{
  MatrixND<nc, nr> result = {};
  for (index_t j = 0; j < nc; ++j) {
    for (index_t i = 0; i < nr; ++i) {
      result(j,i) = (*this)(i,j);
    }
  }
  return result;
}


template <index_t nr, index_t nc, typename T>
template<typename F> inline void
MatrixND<nr, nc, T>::map(F func) const
{
  for (int i=0; i<nr; i++)
    for (int j = 0; j<nc; j++)
      func((*this)(i,j));
}


template <index_t nr, index_t nc, typename T>
template<typename F> inline void
MatrixND<nr, nc, T>::map(F func, MatrixND<nr,nc,T>& destination)
{
  for (int i=0; i<nr; i++)
    for (int j = 0; j<nc; j++)
      destination(i,j) = func((*this)(i,j));
}

//
//
//
template <index_t NR, index_t NC, typename T>
template <int init_row, int init_col, int nr, int nc> 
inline constexpr void
MatrixND<NR,NC,T>::insert(const MatrixND<nr, nc, double> &M) noexcept
{
  constexpr int final_row = init_row + nr - 1;
  constexpr int final_col = init_col + nc - 1;
  static_assert((init_row >= 0) && (final_row < NR) && (init_col >= 0) && (final_col < NC), 
                "MatrixND::insert: init_row, init_col, nr, nc out of bounds");

  Repeat<nc>([&](auto i_) {
      constexpr int i = i_.value;
      constexpr int pos_Cols = init_col + i;
    Repeat<nr>([&](auto j_) {
      constexpr int j = j_.value;
      constexpr int pos_Rows = init_row + j; 
      (*this)(pos_Rows,pos_Cols) = M(j,i);
      });
  });
}


template <index_t NR, index_t NC, typename T>
template <int init_row, int init_col, int nr, int nc> 
inline constexpr void
MatrixND<NR,NC,T>::insert(const MatrixND<nr, nc, double> &M, double fact) noexcept
{

  constexpr int final_row = init_row + nr - 1;
  constexpr int final_col = init_col + nc - 1;
  static_assert((init_row >= 0) && (final_row < NR) && (init_col >= 0) && (final_col < NC), 
                "MatrixND::insert: init_row, init_col, nr, nc out of bounds");

  for (int i=0; i<nc; i++) {
      int pos_Cols = init_col + i;
      for (int j=0; j<nr; j++) {
        int pos_Rows = init_row + j; 
        (*this)(pos_Rows,pos_Cols) = M(j,i)*fact;
      }
  }
}

template <index_t NR, index_t NC, typename T>
template <int nr, int nc> 
inline constexpr void
MatrixND<NR,NC,T>::insert(const MatrixND<nr, nc, double> &M, int init_row, int init_col, double fact) noexcept
{
  [[maybe_unused]] int final_row = init_row + nr - 1;
  [[maybe_unused]] int final_col = init_col + nc - 1; 
  assert((init_row >= 0) && (final_row < NR) && (init_col >= 0) && (final_col < NC));

  for (int i=0; i<nc; i++) {
    int pos_Cols = init_col + i;
    for (int j=0; j<nr; j++) {
      int pos_Rows = init_row + j; 
      (*this)(pos_Rows,pos_Cols) = M(j,i)*fact;
    }
  }
}


template <index_t NR, index_t NC, typename T>
template <int nr, int nc> 
inline constexpr void
MatrixND<NR,NC,T>::assemble(const MatrixND<nr, nc, double> &M, int init_row, int init_col, double fact) noexcept
{
  [[maybe_unused]] int final_row = init_row + nr - 1;
  [[maybe_unused]] int final_col = init_col + nc - 1; 
  assert((init_row >= 0) && (final_row < NR) && (init_col >= 0) && (final_col < NC));

  for (int i=0; i<nc; i++) {
    int pos_Cols = init_col + i;
    for (int j=0; j<nr; j++) {
      (*this)(init_row + j, pos_Cols) += M(j,i)*fact;
    }
  }
}


template <index_t NR, index_t NC, typename T>
template <int nr> inline void
MatrixND<NR,NC,T>::assemble(const VectorND<nr> &v, int init_row, int init_col, double fact) noexcept
{

  [[maybe_unused]] int final_row = init_row + nr - 1;
  assert((init_row >= 0) && (final_row < NR));

  for (int j=0; j<nr; j++)
    (*this)(init_row+j, init_col) += v[j]*fact;
}


template <index_t NR, index_t NC, typename T>
template<int nr, int nc> void
MatrixND<NR,NC,T>::assembleTranspose(const MatrixND<nr, nc, double> &M, int init_row, int init_col, double fact) noexcept
{ 
  {
    [[maybe_unused]] int final_row = init_row + nc - 1;
    [[maybe_unused]] int final_col = init_col + nr - 1; 
    assert((init_row >= 0) && (final_row < NR) && (init_col >= 0) && (final_col < NC));
  }

  for (int i=0; i<nr; i++) {
    int pos_Cols = init_col + i;
    for (int j=0; j<nc; j++) {
      int pos_Rows = init_row + j; 
      (*this)(pos_Rows,pos_Cols) += M(i,j)*fact;
    }
  }
}

template <index_t NR, index_t NC, typename T>
template<int nr> void
MatrixND<NR,NC,T>::assembleTranspose(const VectorND<nr> &v, int init_row, int init_col, double scale) noexcept
{ 
  {
    [[maybe_unused]] int final_col = init_col + nr - 1; 
    assert((init_row >= 0) && (init_col >= 0) && (final_col < NC));
  }

  for (int i=0; i<nr; i++)
    (*this)(init_row, init_col+i) += v[i]*scale;
}

//
// Solve
//
template <index_t nr, index_t nc, typename T> inline int
MatrixND<nr,nc,T>::invert(MatrixND<nr,nc,T> &M) const
{
  static_assert(nr == nc, "Matrix must be square");
  static_assert(std::is_same_v<T,double>, "Only double storage is supported");

  int status = -1;
  if constexpr (nr == 2) {
    cmx_inv2(&(*this)(0,0), &M(0,0), &status);
    return status;
  }
  if constexpr (nr == 3) {
    cmx_inv3(&(*this)(0,0), &M(0,0), &status);
    return status;
  }
  if constexpr (nr == 4) {
    cmx_inv4(&(*this)(0,0), &M(0,0), &status);
    return status;
  }
  if constexpr (nr == 5) {
    cmx_inv5(&(*this)(0,0), &M(0,0), &status);
    return status;
  }
  if constexpr (nr == 6) {
    cmx_inv6(&(*this)(0,0), &M(0,0), &status);
    return status;
  }

  // Use Lapack
  M.zero();
  M.addDiagonal(1.0); // Identity matrix 

  MatrixND<nr, nc, T> work = *this;

  int pivot_ind[nr];
  int nrhs = nr;
  int info = 0;
  int m = nr;
  int n = nc;
  DGESV(&m, 
        &nrhs,
        &work(0,0), &m,
        pivot_ind, 
        &M(0,0), 
        &n, 
        &info);

  if (info != 0)
    status = -std::abs(info);
  return status;
}


template <index_t NR, index_t NC, typename T>
int
MatrixND<NR,NC,T>::solve(const VectorND<NR> &V, VectorND<NR> &res) const noexcept
{
  static_assert(std::is_same_v<T,double>, "Only double storage is supported");
  static_assert(NR == NC);

  if constexpr (NR < 6) {
    MatrixND<NR,NC,T> Ainv;
    int status = this->invert(Ainv);
    if (status != 0)
      return status;

    res = Ainv * V;
    return 0;
  }

  MatrixND<NR,NC> work = *this;
  int pivot_ind[NR];
  int nrhs = 1;
  int nr = NR;
  int nc = NC;
  int info = 0;
  res = V; // X will be overwritten with the solution
  DGESV(&nr, &nrhs, &work(0,0), &nr, pivot_ind, res.values, &nc, &info);
  return -abs(info);
}


template <index_t NR, index_t NC, typename T>
template<index_t n>
int 
MatrixND<NR,NC,T>::solve(const MatrixND<n, n>& M, MatrixND<n, n>& X) const noexcept
{
  static_assert(std::is_same_v<T,double>, "Only double storage is supported");
  static_assert(NR == NC, "Matrix must be square.");
  static_assert(n == NR, "RHS row-count must match A.");

  if constexpr (NR < 6) {
    MatrixND<NR,NC,T> Ainv;
    int status = this->invert(Ainv);
    if (status != 0)
      return status;

    X = Ainv * M;
    return 0;
  }

  MatrixND<NR,NC,T> work = *this;               // copy of A to be factorised
  int ipiv[NR]{};

  int n_eq  = NR;               // order of the system
  int nrhs  = n;                // number of RHS columns
  int lda   = NR;               // leading dim of A
  int ldb   = NR;               // leading dim of X
  int info  = 0;

  X = M;                               // copy RHS, DGESV overwrites
  DGESV(&n_eq, &nrhs,
        &work(0,0), &lda,
        &ipiv[0],
        &X(0,0), &ldb,
        &info);

  return -std::abs(info);
}

template <index_t nr, index_t nc, typename T>
constexpr MatrixND<nr, nc, T>& 
MatrixND<nr, nc, T>::addDiagonal(const double diag) noexcept
{
  for (int i=0; i<nr; i++)
    (*this)(i,i) += diag;

  return *this;
}


template <index_t nr, index_t nc, typename T>
template <class MatT> inline
void MatrixND<nr, nc, T>::addMatrix(const MatT& A, const double scale)
{
  for (int i=0; i<nr; i++)
    for (int j=0; j<nc; j++)
      (*this)(i,j) += A(i,j)*scale;
}

template <index_t nr, index_t nc, typename T>
inline void
MatrixND<nr, nc, T>::addTranspose(const MatrixND<nc,nr>& A, const double scale)
{
  for (int i=0; i<nr; i++)
    for (int j=0; j<nc; j++)
      (*this)(i,j) += A(j,i)*scale;
}


template <index_t nr, index_t nc, typename T>
template <class VecA, class VecB>
constexpr inline
MatrixND<nr, nc, T>&
MatrixND<nr, nc, T>::addTensorProduct(const VecA& a, const VecB& b, const double scale) noexcept
{
  // Chrystal's bun order
  for (int j=0; j<nc; j++)
    for (int i=0; i<nr; i++)
      (*this)(i,j) += a[i]*b[j]*scale;
  return *this;
}


template <index_t nr, index_t nc, typename T> 
template <index_t nk> 
constexpr void
MatrixND<nr, nc, T>::addMatrixProduct(const MatrixND<nr, nk, T>& A, 
                                      const MatrixND<nk, nc>& B,
                                      double scale) noexcept
{
  static_assert(std::is_same_v<T,double>, "Only double storage is supported");

  if constexpr (nr*nk < BlasSize) {
    // want: this += A * B * scale
    const double *ckjPtr  = &B(0,0);
    for (int j=0; j<nc; j++) {
      double * RESTRICT aijPtrA = &values[j*nr];
      for (int k=0; k<nk; k++) {
        double tmp = *ckjPtr++ * scale;
        double * aijPtr = aijPtrA;
        const double * RESTRICT bikPtr = &A.values[k*nr];
        for (int i=0; i<nr; i++)
          *aijPtr++ += *bikPtr++ * tmp;
      }
    }
  }
  else
  {
    int m = nr,
        n = nc,
        k = nk;
    double one = 1.0;
    DGEMM("N", "N", &m, &n, &k, &scale, 
                                const_cast<double*>(&A(0,0)), &m,
                                const_cast<double*>(&B(0,0)), &k,
                                &one,   &(*this)(0,0),        &m);
  }
}

template <index_t nr, index_t nc, typename T> 
template <index_t nk> constexpr
void
MatrixND<nr, nc, T>::setMatrixProduct(const MatrixND<nr, nk, T>& B, 
                                      const MatrixND<nk, nc, T>& C,
                                      double scale) noexcept
{
  this->zero();
  if constexpr (nr*nk < BlasSize) {
    // want: this = B * C * scale
    const double * RESTRICT ckjPtr  = C.data();
    for (int j=0; j<nc; j++) {
      double * aijPtrA = data() + j*nr;
      for (int k=0; k<nk; k++) {
        double tmp = (*ckjPtr++) * scale;
        double * aijPtr = aijPtrA;
        const double * RESTRICT bikPtr = B.data() + k*nr;
        for (int i=0; i<nr; i++)
          *aijPtr++ += *bikPtr++ * tmp;
      }
    }
  }
  else
  {
    static_assert(std::is_same_v<T,double>, "Only double storage is supported");
  
    int m = nr,
        n = nc,
        k = nk;
    double zero = 0.0;
    DGEMM("N", "N", &m, &n, &k, &scale, 
                                const_cast<double*>(&B(0,0)), &m,
                                const_cast<double*>(&C(0,0)), &k,
                                &zero,   &(*this)(0,0),       &m);
  }
}

template <index_t nr, index_t nc, typename T> 
template <class MatT, int nk> inline
void
MatrixND<nr, nc, T>::addMatrixProduct(const MatrixND<nr, nk, T>& A, const MatT& B, double scale)
{
  if constexpr (nr*nc < BlasSize) {
    Repeat<nr> ([&](auto i_) {
      constexpr static int i = i_.value;
      Repeat<nc> ([&](auto j_) {
        constexpr static int j = j_.value;
        for (int k=0; k < nk; k++)
          (*this)(i,j) += scale*A(i,k)*B(k,j);
      });
    });
  }
  else
  {
    static_assert(std::is_same_v<T,double>, "Only double storage is supported");
  
    int m = nr,
        n = nc,
        k = nk;
    double one = 1.0;
    DGEMM("N", "N", &m, &n, &k, &scale, 
                                const_cast<double*>(&A(0,0)), &m,
                                const_cast<double*>(&B(0,0)), &k,
                                &one,   &(*this)(0,0),        &m);
  }
}


#if 1
template <index_t nr, index_t nc, typename T> 
template <class MatT, int nk> inline
void
MatrixND<nr, nc, T>::addMatrixProduct(double scale_this, 
                                      const MatrixND<nr, nk, T>& A, 
                                      const MatT& B, double scale)
{
  static_assert(std::is_same_v<T,double>, "Only double storage is supported");

  int m = nr,
      n = nc,
      k = nk;
  DGEMM("N", "N", &m, &n, &k, &scale, 
                              const_cast<double*>(&A(0,0)), &m,
                              const_cast<double*>(&B(0,0)), &k,
                              &scale_this,   &(*this)(0,0), &m);
}
#endif


// A  = a*B'*C + b*A
// A  = a*B'*C
// A += a*B'*C
template <index_t nr, index_t nc, typename T> 
template <int nk>
inline constexpr void
MatrixND<nr, nc, T>::setMatrixTransposeProduct(const MatrixND<nk, nr, T>& B,
                                               const MatrixND<nk, nc, T>& C) noexcept
{
  if constexpr (nr*nc >= BlasSize) {
    this->zero();
    double thisFact = 0.0;
    double otherFact = 1.0;
    int m = nr,
        n = nc,
        k = nk;
    DGEMM("T", "N", &m, &n, &k, &otherFact, 
                                const_cast<double*>(&B(0,0)), &k,
                                const_cast<double*>(&C(0,0)), &k,
                                &thisFact,   &(*this)(0,0),   &m);
    return;
  }

  else {
#if 0
    double *RESTRICT aij =this->data();
    for (int j=0; j<nc; j++) {
      for (int i=0; i<nr; i++) {
        const double *RESTRICT bki  = &(&B(0,0))[i*nk];
        const double *RESTRICT cjk  = &(&C(0,0))[j*nk];
        double sum = 0.0;
        for (int k=0; k<nk; k++)
          sum += (*bki++) * (*cjk++);

        *aij++ = sum;
      }
    }
  }
#else 
    double *RESTRICT aij =this->data();
    Repeat<nc> ([&](auto j_) {
      constexpr static int j = j_.value;
      Repeat<nr> ([&](auto i_) {
        constexpr static int i = i_.value;
        const double *RESTRICT bki  = &(&B(0,0))[i*nk];
        const double *RESTRICT cjk  = &(&C(0,0))[j*nk];
        double sum = 0.0;
        for (int k=0; k<nk; k++)
          sum += (*bki++) * (*cjk++);

        *aij++ = sum;
      });
    });
  }
#endif
}

// A  = a*B'*C + b*A
// A  = a*B'*C
// A += a*B'*C
template <index_t nr, index_t nc, typename T> 
template <int nk>
constexpr void
MatrixND<nr, nc, T>::addMatrixTransposeProduct(const MatrixND<nk, nr, T>& B,
                                               const MatrixND<nk, nc, T>& C) noexcept
{
  if constexpr (nr*nc >= BlasSize) {
    double thisFact = 1.0;
    double otherFact = 1.0;
    int m = nr,
        n = nc,
        k = nk;
    DGEMM("T", "N", &m, &n, &k, &otherFact, 
                                const_cast<double*>(&B(0,0)), &k,
                                const_cast<double*>(&C(0,0)), &k,
                                &thisFact,   &(*this)(0,0),   &m);
    return;
  }

  else {
    double *RESTRICT aijPtr =this->data();
    Repeat<nc> ([&](auto j_) {
      constexpr static int j = j_.value;
      Repeat<nr> ([&](auto i_) {
        constexpr static int i = i_.value;
        const double *RESTRICT bkiPtr  = &(&B(0,0))[i*nk];
        const double *RESTRICT cjkPtr  = &(&C(0,0))[j*nk];
        double sum = 0.0;
        for (int k=0; k<nk; k++) {
          sum += *bkiPtr++ * *cjkPtr++;
        }
        *aijPtr++ += sum;
      });
    });
  }
}

// A  = a*B'*C + b*A
// A  = a*B'*C
// A += a*B'*C
template <index_t nr, index_t nc, typename T> 
template <int nk>
constexpr void
MatrixND<nr, nc, T>::addMatrixTransposeProduct(double scale,
                                               const MatrixND<nk, nr, T>& B,
                                               const MatrixND<nk, nc, T>& C) noexcept
{
  if constexpr (nr*nk >= BlasSize) {
    double otherFact = 1.0;
    int m = nr,
        n = nc,
        k = nk;
    DGEMM("T", "N", &m, &n, &k, &otherFact, 
                                const_cast<double*>(&B(0,0)), &k,
                                const_cast<double*>(&C(0,0)), &k,
                                &scale,   &(*this)(0,0),      &m);
    return;
  }
  else {
    double *RESTRICT aijPtr =this->data();
    for (int j=0; j<nc; j++) {
      for (int i=0; i<nr; i++) {
        const double *RESTRICT bkiPtr  = &(&B(0,0))[i*nk];
        const double *RESTRICT cjkPtr  = &(&C(0,0))[j*nk];
        double sum = 0.0;
        for (int k=0; k<nk; k++) {
          sum += *bkiPtr++ * *cjkPtr++;
        }
        *aijPtr = *aijPtr * scale + sum;
        aijPtr++;
      }
    }
  }
}

// A  = a*B'*C + b*A
// A  = a*B'*C
// A += a*B'*C
template <index_t nr, index_t nc, typename T> 
template <class MatT, int nk> inline
void
MatrixND<nr, nc, T>::addMatrixTransposeProduct(double thisFact,
                                               const MatrixND<nk, nr, T>& B,
                                               const MatT& C,
                                               double otherFact)
{
  if constexpr (nr*nk >= BlasSize) {
    int m = nr,
        n = nc,
        k = nk;
    DGEMM("T", "N", &m, &n, &k, &otherFact, 
                                const_cast<double*>(&B(0,0)), &k,
                                const_cast<double*>(&C(0,0)), &k,
                                &thisFact,   &(*this)(0,0),   &m);
    return;
  }

  else {
#ifdef MATRIX_BRANCHING
    if (thisFact == 1.0) {
      double *RESTRICT aijPtr =this->data();
      for (int j=0; j<nc; j++) {
        for (int i=0; i<nr; i++) {
          const double *RESTRICT bkiPtr  = &(&B(0,0))[i*nk];
          const double *RESTRICT cjkPtr  = &(&C(0,0))[j*nk];
          double sum = 0.0;
          for (int k=0; k<nk; k++) {
            sum += *bkiPtr++ * *cjkPtr++;
          }
          *aijPtr++ += sum * otherFact;
        }
      } 
    }
    else if (thisFact == 0.0) {
      double *aijPtr =this->data();
      for (int j=0; j<nc; j++) {
        for (int i=0; i<nr; i++) {
          const double *RESTRICT bkiPtr  = &(&B(0,0))[i*nk];
          const double *RESTRICT cjkPtr  = &(&C(0,0))[j*nk];
          double sum = 0.0;
          for (int k=0; k<nk; k++) {
            sum += *bkiPtr++ * *cjkPtr++;
          }
          *aijPtr++ = sum * otherFact;
        }
      } 
    }
    else 
#endif
    {
      double *RESTRICT aijPtr =this->data();
      for (int j=0; j<nc; j++) {
        for (int i=0; i<nr; i++) {
          const double *RESTRICT bkiPtr  = &(&B(0,0))[i*nk];
          const double *RESTRICT cjkPtr  = &(&C(0,0))[j*nk];
          double sum = 0.0;
          for (int k=0; k<nk; k++) {
            sum += *bkiPtr++ * *cjkPtr++;
          }
          *aijPtr = *aijPtr * thisFact + sum * otherFact;
          aijPtr++;
        }
      }
    }
  }
}

// // A'BA
// template <int nr, int nc, class scalar_t> 
// template <int ncB> inline
// int
// MatrixND<nr,nc,scalar_t>::addMatrixTripleProduct(
//                                const  MatrixND<ncB, nr, scalar_t> &T, 
//                                const  MatrixND<ncB, ncB, scalar_t> &B)
// {
// // #ifdef MATRIX_BRANCHING
//   // if (otherFact == 0.0 && thisFact == 1.0)
//   //   return 0;
// // #endif 

//   MatrixND<ncB, nc> BT;
//   BT.zero();
//   BT.addMatrixProduct(B, T);
//   this->addMatrixTransposeProduct(T, BT); // , 1.0);
//   return 0;
// }

// A'BA
template <int nr, int nc, class scalar_t> 
template <int ncB> inline
int
MatrixND<nr,nc,scalar_t>::addMatrixTripleProduct( 
                               double thisFact,
                               const  MatrixND<ncB, nr, scalar_t> &T, 
                               const  MatrixND<ncB, ncB, scalar_t> &B, 
                               double otherFact)
{
// #ifdef MATRIX_BRANCHING
  // if (otherFact == 0.0 && thisFact == 1.0)
  //   return 0;
// #endif 

  MatrixND<ncB, nc> BT{};
  BT.addMatrixProduct(B, T, otherFact);
  this->addMatrixTransposeProduct(thisFact, T, BT , 1.0);
  return 0;
}


template <int nr, int nc, class scalar_t> 
template <int nk, int nl> inline int 
MatrixND<nr,nc,scalar_t>::addMatrixTripleProduct(double thisFact, 
                           const MatrixND<nk,nr> &A, 
                           const MatrixND<nk,nl> &B, 
                           const MatrixND<nl,nc> &C, double otherFact)
{
// #ifdef MATRIX_BRANCHING
//   if (otherFact == 0.0 && thisFact == 1.0)
//     return 0;
// #endif
  MatrixND<nk, nc> BC {};
  BC.addMatrixProduct(B, C, otherFact);
  this->addMatrixTransposeProduct(thisFact, A, BC , 1.0);
  return 0;
}




//
// 3D
//

template<int NR, int NC, typename T>
template<class VecT>
constexpr MatrixND<NR,NC,T>& 
MatrixND<NR,NC,T>::addSpin(const VecT& v) noexcept
{
  static_assert(NR == 3 && NC == 3, "addSpin requires a 3x3 matrix");

   const double v0 = v[0],
                v1 = v[1],
                v2 = v[2];

  (*this)(0, 0) += 0.0;   (*this)(0, 1) += -v2;     (*this)(0, 2) +=  v1;
  (*this)(1, 0) +=  v2;   (*this)(1, 1) +=  0.0;    (*this)(1, 2) += -v0;
  (*this)(2, 0) += -v1;   (*this)(2, 1) +=  v0;     (*this)(2, 2) += 0.0;

  return *this;
}


template<int NR, int NC, typename T>
template<class VecT>
constexpr MatrixND<NR,NC,T>&
MatrixND<NR,NC,T>::addSpin(const VecT& v, double mult) noexcept
{
   const double v0 = mult*v[0],
                v1 = mult*v[1],
                v2 = mult*v[2];

  (*this)(0, 0) += 0.0;   (*this)(0, 1) += -v2;     (*this)(0, 2) +=  v1;
  (*this)(1, 0) +=  v2;   (*this)(1, 1) += 0.00;    (*this)(1, 2) += -v0;
  (*this)(2, 0) += -v1;   (*this)(2, 1) +=  v0;     (*this)(2, 2) += 0.0;
  return *this;
}


template<int NR, int NC, typename T>
template <class VecT> 
constexpr
MatrixND<NR,NC,T>& 
MatrixND<NR,NC,T>::addSpinSquare(const VecT& v, const double scale) noexcept
{
  static_assert(NR == 3 && NC == 3, "addSpinSquare requires a 3x3 matrix");

  const double v1 = v[0],
               v2 = v[1],
               v3 = v[2];

  (*this)(0,0) += scale*( -v2*v2 - v3*v3 );
  (*this)(1,1) += scale*( -v1*v1 - v3*v3 );
  (*this)(2,2) += scale*( -v1*v1 - v2*v2 );

  (*this)(0,1) += scale*(  v1*v2 );
  (*this)(1,0) += scale*(  v1*v2 );
  (*this)(2,0) += scale*(  v1*v3 );
  (*this)(0,2) += scale*(  v1*v3 );
  (*this)(1,2) += scale*(  v2*v3 );
  (*this)(2,1) += scale*(  v2*v3 );
  return *this;
}


template<int NR, int NC, typename T>
template<class VecT> 
constexpr void 
MatrixND<NR,NC,T>::addSpinProduct(const VecT& a, const VectorND<NR,T>& b, const double scale) noexcept
{
  // a^b^ = boa - a.b 1
  // where 'o' denotes the tensor product and '.' the dot product
  //
  static_assert(NR == 3 && NC == 3, "Matrix must be 3x3");
  this->addTensorProduct(b, a, scale);
  this->addDiagonal(-b.dot(a)*scale);
}

template<int NR, int NC, typename T>
template<class VecT>
constexpr void 
MatrixND<NR,NC,T>::addMatrixSpinProduct(const MatrixND<NR,NC,T>& A, const VecT& b, const double scale) noexcept
{
  // this += s*A*[b^]
  // where b^ is the skew-symmetric representation of the three-vector b, s is a scalar,
  // and A a 3x3 matrix.
  //
  (*this)(0, 0) += scale*( A(0,1)*b[2] - A(0,2)*b[1]);
  (*this)(0, 1) += scale*(-A(0,0)*b[2] + A(0,2)*b[0]);
  (*this)(0, 2) += scale*( A(0,0)*b[1] - A(0,1)*b[0]);
  (*this)(1, 0) += scale*( A(1,1)*b[2] - A(1,2)*b[1]);
  (*this)(1, 1) += scale*(-A(1,0)*b[2] + A(1,2)*b[0]);
  (*this)(1, 2) += scale*( A(1,0)*b[1] - A(1,1)*b[0]);
  (*this)(2, 0) += scale*( A(2,1)*b[2] - A(2,2)*b[1]);
  (*this)(2, 1) += scale*(-A(2,0)*b[2] + A(2,2)*b[0]);
  (*this)(2, 2) += scale*( A(2,0)*b[1] - A(2,1)*b[0]);
}

template<int NR, int NC, typename T>
template<class MatT>
constexpr void 
MatrixND<NR,NC,T>::addSpinMatrixProduct(const VectorND<NR,T>& a, const MatT& B, const double scale) noexcept
{
  static_assert(NR == 3 && NC == 3, "Matrix must be 3x3");
  // this += s*[a^]*B
  // where a^ is the skew-symmetric representation of the three-vector a, s is a scalar,
  // and B a 3x3 matrix.
  //
  (*this)(0, 0) += scale*( -B(1,0)*a[2] + B(2,0)*a[1]);
  (*this)(0, 1) += scale*( -B(1,1)*a[2] + B(2,1)*a[1]);
  (*this)(0, 2) += scale*( -B(1,2)*a[2] + B(2,2)*a[1]);
  (*this)(1, 0) += scale*(  B(0,0)*a[2] - B(2,0)*a[0]);
  (*this)(1, 1) += scale*(  B(0,1)*a[2] - B(2,1)*a[0]);
  (*this)(1, 2) += scale*(  B(0,2)*a[2] - B(2,2)*a[0]);
  (*this)(2, 0) += scale*( -B(0,0)*a[1] + B(1,0)*a[0]);
  (*this)(2, 1) += scale*( -B(0,1)*a[1] + B(1,1)*a[0]);
  (*this)(2, 2) += scale*( -B(0,2)*a[1] + B(1,2)*a[0]);
}
} // namespace OpenSees
