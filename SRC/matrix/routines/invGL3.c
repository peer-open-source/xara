//
// adapted from https://caps.gsfc.nasa.gov/simpson/software/m33inv_f90.txt
//
#include <math.h>

#if 0
int 
cmx_inv3(double *a, double *ainv, int *ok_flag__)
{
/* **************************************************************************************** */
/*  m33inv  -  compute the inverse of a 3x3 matrix. */

/*  a       = input 3x3 matrix to be inverted */
/*  ainv    = output 3x3 inverse of matrix a */
/*  ok_flag = (output) .true. if the input matrix could be inverted, and */
/*            .false. if the input matrix is singular. */
/* **************************************************************************************** */


    /* Parameter adjustments */
    ainv -= 4;
    a -= 4;

    const double det = a[4]*a[8]*a[12] 
                     - a[4]*a[11]*a[9] 
                     - a[7]*a[5]*a[12] 
                     + a[7]*a[11]*a[6] 
                     + a[10]*a[5]*a[9] 
                     - a[10]*a[8]*a[6];

    const double eps=1e-10;
    if (fabs(det) <= eps) {
        *ok_flag__ = -1;
        // return 0;
    }

    double cofactor[9];
    cofactor[0] =    a[8]*a[12] - a[11]*a[9];
    cofactor[3] = -(a[5]*a[12] - a[11]*a[6]);
    cofactor[6] =    a[5]*a[9] - a[8]*a[6];
    cofactor[1] = -(a[7]*a[12] - a[10]*a[9]);
    cofactor[4] =   a[4]*a[12] - a[10]*a[6];
    cofactor[7] = -(a[4]*a[9] - a[7]*a[6]);
    cofactor[2] =   a[7]*a[11] - a[10]*a[8];
    cofactor[5] = -(a[4]*a[11] - a[10]*a[5]);
    cofactor[8] =   a[4]*a[8] - a[7]*a[5];

    for (int i = 1; i <= 3; ++i)
        for (int j = 1; j <= 3; ++j)
            ainv[j + i * 3] = cofactor[i + j * 3 - 4] / det;

    *ok_flag__ = 0;
    return 0;
}
#endif

#include <math.h>
#include <float.h>

#ifndef CMX_INV3_PIVOT_TOL
#define CMX_INV3_PIVOT_TOL (32.0 * DBL_EPSILON)
#endif

int
cmx_inv3(double *a, double *ainv, int *ok_flag__)
{
    /*
      Input and output are 0-based packed column-major 3x3 matrices.

          a[row + 3*col]
          ainv[row + 3*col]

      The input matrix a is not overwritten.
      The output ainv is written only on success.

      ok_flag__:
          0   success
         -1   singular or numerically singular
    */

    double w00 = a[0], w10 = a[1], w20 = a[2];
    double w01 = a[3], w11 = a[4], w21 = a[5];
    double w02 = a[6], w12 = a[7], w22 = a[8];

    double scale = fabs(w00);
    double v;

    v = fabs(w10); if (v > scale) scale = v;
    v = fabs(w20); if (v > scale) scale = v;
    v = fabs(w01); if (v > scale) scale = v;
    v = fabs(w11); if (v > scale) scale = v;
    v = fabs(w21); if (v > scale) scale = v;
    v = fabs(w02); if (v > scale) scale = v;
    v = fabs(w12); if (v > scale) scale = v;
    v = fabs(w22); if (v > scale) scale = v;

    if (!(scale > 0.0)) {
        *ok_flag__ = -1;
        return 0;
    }

    const double tol = CMX_INV3_PIVOT_TOL * scale;

    /*
      Step 0: pivot in column 0.
    */

    int p0 = 0;
    double amax = fabs(w00);

    v = fabs(w10);
    if (v > amax) {
        amax = v;
        p0 = 1;
    }

    v = fabs(w20);
    if (v > amax) {
        amax = v;
        p0 = 2;
    }

    if (!(amax > tol)) {
        *ok_flag__ = -1;
        return 0;
    }

    if (p0 == 1) {
        double t;

        t = w00; w00 = w10; w10 = t;
        t = w01; w01 = w11; w11 = t;
        t = w02; w02 = w12; w12 = t;
    } else if (p0 == 2) {
        double t;

        t = w00; w00 = w20; w20 = t;
        t = w01; w01 = w21; w21 = t;
        t = w02; w02 = w22; w22 = t;
    }

    const double inv0 = 1.0 / w00;

    w10 *= inv0;
    w20 *= inv0;

    w11 -= w10 * w01;
    w12 -= w10 * w02;

    w21 -= w20 * w01;
    w22 -= w20 * w02;

    /*
      Step 1: pivot in column 1.
    */

    int p1 = 1;
    amax = fabs(w11);

    v = fabs(w21);
    if (v > amax) {
        amax = v;
        p1 = 2;
    }

    if (!(amax > tol)) {
        *ok_flag__ = -1;
        return 0;
    }

    if (p1 == 2) {
        double t;

        t = w10; w10 = w20; w20 = t;
        t = w11; w11 = w21; w21 = t;
        t = w12; w12 = w22; w22 = t;
    }

    const double inv1 = 1.0 / w11;

    w21 *= inv1;
    w22 -= w21 * w12;

    /*
      Step 2.
    */

    if (!(fabs(w22) > tol)) {
        *ok_flag__ = -1;
        return 0;
    }

    const double inv2 = 1.0 / w22;

    /*
      Solve A X = I using the LU factors.

      First apply the same row interchanges to I.
      The variables rXY mean row X, column Y.
    */

    double r00 = 1.0, r10 = 0.0, r20 = 0.0;
    double r01 = 0.0, r11 = 1.0, r21 = 0.0;
    double r02 = 0.0, r12 = 0.0, r22 = 1.0;

    if (p0 == 1) {
        double t;

        t = r00; r00 = r10; r10 = t;
        t = r01; r01 = r11; r11 = t;
        t = r02; r02 = r12; r12 = t;
    } else if (p0 == 2) {
        double t;

        t = r00; r00 = r20; r20 = t;
        t = r01; r01 = r21; r21 = t;
        t = r02; r02 = r22; r22 = t;
    }

    if (p1 == 2) {
        double t;

        t = r10; r10 = r20; r20 = t;
        t = r11; r11 = r21; r21 = t;
        t = r12; r12 = r22; r22 = t;
    }

    /*
      Forward solve: L Y = P I.

      L =
          1    0    0
          w10  1    0
          w20  w21  1
    */

    r10 -= w10 * r00;
    r11 -= w10 * r01;
    r12 -= w10 * r02;

    r20 -= w20 * r00;
    r21 -= w20 * r01;
    r22 -= w20 * r02;

    r20 -= w21 * r10;
    r21 -= w21 * r11;
    r22 -= w21 * r12;

    /*
      Back solve: U X = Y.

      U =
          w00  w01  w02
          0    w11  w12
          0    0    w22
    */

    r20 *= inv2;
    r21 *= inv2;
    r22 *= inv2;

    r10 = (r10 - w12 * r20) * inv1;
    r11 = (r11 - w12 * r21) * inv1;
    r12 = (r12 - w12 * r22) * inv1;

    r00 = (r00 - w01 * r10 - w02 * r20) * inv0;
    r01 = (r01 - w01 * r11 - w02 * r21) * inv0;
    r02 = (r02 - w01 * r12 - w02 * r22) * inv0;

    ainv[0] = r00;
    ainv[1] = r10;
    ainv[2] = r20;

    ainv[3] = r01;
    ainv[4] = r11;
    ainv[5] = r21;

    ainv[6] = r02;
    ainv[7] = r12;
    ainv[8] = r22;

    *ok_flag__ = 0;
    return 0;
}
