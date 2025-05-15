//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Matrix exponential
//
// References:
//   N. J. Higham, The scaling and squaring method for the matrix
//      exponential revisited. SIAM J. Matrix Anal. Appl., 26(4), (2005),
//      pp. 1179-1193.
//   A. H. Al-Mohy and N. J. Higham, A new scaling and squaring algorithm
//      for the matrix exponential, SIAM J. Matrix Anal. Appl., 31(3),
//      (2009), pp. 970-989.
//
// See also
//
//   Golub, G. H. and C. F. Van Loan, Matrix Computation, p. 384, Johns Hopkins University Press, 1983.
// 
//   Moler, C. B. and C. F. Van Loan, 
//      "Nineteen Dubious Ways to Compute the Exponential of a Matrix," 
//      SIAM Review 20, 1978, pp. 801–836. 
//      Reprinted and updated as 
//      "Nineteen Dubious Ways to Compute the Exponential of a Matrix, Twenty-Five Years Later,”
//      SIAM Review 45, 2003, pp. 3–49.
//
//===----------------------------------------------------------------------===//
//
// Claudio M. Perez
//
// Adapted from:
//   https://eigen.tuxfamily.org/dox/unsupported/group__MatrixFunctions__Module.html
//
#pragma once
#include <cmath>
#include <complex>

namespace OpenSees {
namespace Internal {

//  Compute the (3,3)-Pade approximant to the exponential.
//
//  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pade
//  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
//
template <typename MatA, typename MatU, typename MatV>
void
matrix_exp_pade3(const MatA& A, MatU& U, MatV& V) {
  const double b[] = {120.L, 60.L, 12.L, 1.L};
  const MatA A2 = A * A;
  MatA tmp = b[3] * A2;
  tmp.addDiagonal(b[1]);
  U = A * tmp;
  V = b[2] * A2;
  V.addDiagonal(b[0]);
}

//  Compute the (5,5)-Pade approximant to the exponential.
//
//  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pade
//  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
// 
template <typename MatA, typename MatU, typename MatV>
void
matrix_exp_pade5(const MatA& A, MatU& U, MatV& V) {
  const double b[] = {30240.L, 15120.L, 3360.L, 420.L, 30.L, 1.L};
  const MatA A2 = A * A;
  const MatA A4 = A2 * A2;
  MatA tmp = b[5] * A4 + b[3] * A2;
  tmp.addDiagonal(b[1]);
  U = A * tmp;
  V = b[4] * A4 + b[2] * A2;
  V.addDiagonal(b[0]);
}

//  Compute the (7,7)-Pade approximant to the exponential.
//
//  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pade
//  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
// 
template <typename MatA, typename MatU, typename MatV>
void
matrix_exp_pade7(const MatA& A, MatU& U, MatV& V) {
  const double b[] = {17297280.L, 8648640.L, 1995840.L, 277200.L, 25200.L, 1512.L, 56.L, 1.L};
  const MatA A2 = A  * A;
  const MatA A4 = A2 * A2;
  const MatA A6 = A4 * A2;
  MatA tmp = b[7] * A6 + b[5] * A4 + b[3] * A2;
  tmp.addDiagonal(b[1]);
  U = A * tmp;
  V = b[6] * A6 + b[4] * A4 + b[2] * A2;
  V.addDiagonal(b[0]);
}

//  Compute the (9,9)-Pade approximant to the exponential.
//
//  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pade
//  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
// 
template <typename MatA, typename MatU, typename MatV>
void matrix_exp_pade9(const MatA& A, MatU& U, MatV& V) {
  const double b[] = {17643225600.L, 8821612800.L, 2075673600.L, 302702400.L, 30270240.L,
                          2162160.L,     110880.L,     3960.L,       90.L,        1.L};
  const MatA A2 = A * A;
  const MatA A4 = A2 * A2;
  const MatA A6 = A4 * A2;
  const MatA A8 = A6 * A2;
  MatA tmp = b[9] * A8 + b[7] * A6 + b[5] * A4 + b[3] * A2;
  tmp.addDiagonal(b[1]);

  U = A * tmp;
  V = b[8] * A8 + b[6] * A6 + b[4] * A4 + b[2] * A2;
  V.addDiagonal(b[0]);
}

//  Compute the (13,13)-Pade approximant to the exponential.
//
//  After exit, \f$ (V+U)(V-U)^{-1} \f$ is the Pade
//  approximant of \f$ \exp(A) \f$ around \f$ A = 0 \f$.
// 
template <typename MatA, typename MatU, typename MatV>
void matrix_exp_pade13(const MatA& A, MatU& U, MatV& V) {
  const double b[] = {64764752532480000.L,
                      32382376266240000.L,
                      7771770303897600.L,
                      1187353796428800.L,
                      129060195264000.L,
                      10559470521600.L,
                      670442572800.L,
                      33522128640.L,
                      1323241920.L,
                      40840800.L,
                      960960.L,
                      16380.L,
                      182.L,
                      1.L};
  const MatA A2 = A * A;
  const MatA A4 = A2 * A2;
  const MatA A6 = A4 * A2;
  V = b[13] * A6 + b[11] * A4 + b[9] * A2;  // used for temporary storage

  MatA tmp = A6 * V;
  tmp += b[7] * A6 + b[5] * A4 + b[3] * A2;
  tmp.addDiagonal(b[1]);

  U   = A * tmp;
  tmp = b[12] * A6 + b[10] * A4 + b[8] * A2;
  V  = A6 * tmp;
  V += b[6] * A6 + b[4] * A4 + b[2] * A2;
  V.addDiagonal(b[0]);
}


template <typename MatrixType>
static void
matrix_exp_computeUV(const MatrixType& arg, MatrixType& U, MatrixType& V, int& squarings)
{

  // TODO : This can surly be optimized
//const double l1norm = arg.cwiseAbs().colwise().sum().maxCoeff();
  double l1norm=0;
  arg.map([&l1norm](double x)->void {
      l1norm = std::abs(l1norm) > std::abs(x) ? l1norm : x;
  });


  squarings = 0;
  if (l1norm < 1.495585217958292e-002) {
    matrix_exp_pade3(arg, U, V);

  } else if (l1norm < 2.539398330063230e-001) {
    matrix_exp_pade5(arg, U, V);

  } else if (l1norm < 9.504178996162932e-001) {
    // m = 7
    matrix_exp_pade7(arg, U, V);

  } else if (l1norm < 2.097847961257068e+000) {
    // m = 9
    matrix_exp_pade9(arg, U, V);

  } else { 
    // m = 13
    const double maxnorm = 5.371920351148152;
    std::frexp(l1norm / maxnorm, &squarings);
    if (squarings < 0)
      squarings = 0;

    auto func = [=](double x) ->double {
      return std::ldexp(x, -squarings);
    };

    MatrixType A = arg;
    A.map(func, A);
    matrix_exp_pade13(A, U, V);
  }
}

}  // namespace Internal

template <typename MatrixType>
static inline MatrixType
ExpGLn(const MatrixType& arg)
{
  // Pade approximant is (U+V) / (-U+V)
  MatrixType U, V;
  int squarings;
  Internal::matrix_exp_computeUV<MatrixType>(arg, U, V, squarings);
  MatrixType numer =  U + V;
  MatrixType denom =  V - U;

  MatrixType result;
  Matrix rhs(result);
  Matrix lhs(numer);
  denom.solve(lhs, rhs);

  // Undo scaling by repeated squaring
  for (int i = 0; i < squarings; i++)
    result = result*result;

  return result;
}

} // end namespace OpenSees

#include <Rotations.hpp>

static inline int
test_expm()
{
  using namespace OpenSees;
  Vector3D x = {0, 0, 2.2};

  MatrixND<3,3> X = Hat(x);
  MatrixND<3,3> B = ExpGLn(X);
  double cs = std::cos(x[2]),
         sn = std::sin(x[2]);

  MatrixND<3,3> A = {{{cs, -sn, 0},
                      {sn,  cs, 0},
                      { 0,   0, 1}}};

  opserr << Matrix(A) ;
  opserr << Matrix(B) ;
  return 0;
}

