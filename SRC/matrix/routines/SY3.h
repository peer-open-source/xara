
/* Eigen-decomposition for symmetric 3x3 real matrices.
   Public domain, copied from the public domain Java library JAMA. */

#ifndef _eigSY3_h
#define _eigSY3_h
#ifdef __cplusplus
extern "C" {
#endif // __cplusplus
/* Symmetric matrix A => eigenvectors in columns of V, corresponding
   eigenvalues in d. */

int cmx_eigSY3(double A[3][3], double V[3][3], double d[3]);

// int cmx_eigSY3(double** A /*[3][3]*/, double** V /*[3][3]*/, double* d /*[3]*/);
int cmx_eigSY3v2(double A[3][3], double EE[3][3], double V[3][3], double d[3]);

#ifdef __cplusplus
}
#endif // __cplusplus
#endif
