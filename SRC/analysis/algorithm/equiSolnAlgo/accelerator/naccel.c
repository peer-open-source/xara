//
// Routine used by MillerAccelerator
//

// naccel.f -- translated by f2c (version 20200916).
// You must link the resulting object file with libf2c:
//      on Microsoft Windows system, link with libf2c.lib;
//      on Linux or Unix systems, link with .../path/to/libf2c.a -lm
//      or, if you install libf2c.a in a standard place, with -lf2c -lm
//      -- in that order, at the end of the command line, as in
//              cc *.o -lf2c -lm
//      Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,
//
//              http://www.netlib.org/f2c/libf2c.zip

#include <math.h>
#define LINK_MAX 10


int naccel_(int *n_, int *itr, int *mvec, 
        double *tol, double *u, double *f)
{
    // System generated locals
    int u_offset;

    // Local variables
    static double c__[LINK_MAX],
                  h__[(LINK_MAX+1)*(LINK_MAX+1)]; // was [11][11]
    static int k;
    static double t;
    static int tmp, head, nvec, link[LINK_MAX+1], 
               last, next, jptr, kptr, km1ptr;

// ***********************************************************************

//   NACCEL -- Newton iteration accelerator.

//   Argument  I/O/M/T  Description
//   --------  -------  -----------
//     N          I     Vector size.

//     ITR        I     Newton iteration count.  Initialization is done
//                      for ITR=1 and ITR is ignored thereafter.

//     MVEC       I     Maximum number of vectors to use in the
//                      acceleration algorithm.  May change from call
//                      to call but should be greater than 0 and no
//                      more than the internal parameter MAX (=10).

//     TOL        I     Tolerance for dropping vectors.  We drop the
//                      pair (z_k,w_k) if the sine of the angle between
//                      w_k and span{w_1, ..., w_(k-1)} is less than TOL.

//     U          T     Work space of length at least N*(2*MVEC+2).
//                      It should not be modified between calls.

//     F          M     Vector of length N.  On entry, it is the value
//                      of the function f at the current iterate.  It
//                      is overwritten with the accelerated correction.
//                      The unaccelerated correction would simply be f
//                      itself.

// ***********************************************************************
//     ==================================================================
//      First call: save f and the (unaccelerated) correction.
//     ==================================================================
    const int n = *n_;
    // Parameter adjustments
    --f;
    u_offset = 1 + n*3;
    u -= u_offset;

    // Function Body
    if (*itr == 1) {
        head = 1;
        for (int j = 1; j <= n; ++j) {
            u[j + ((head << 1) + 1) * n] = f[j];
            u[j + ((head << 1) + 2) * n] = f[j];
        }
        link[0] = 0;
        nvec = 1;
//      Free storage linked list.
        next = 2;
        for (k = 2; k <= LINK_MAX; ++k)
          link[k - 1] = k + 1;

        link[LINK_MAX] = 0;
        return 0;
    }

// ==================================================================
//  Compute w_1.
// ==================================================================

    for (int j = 1; j <= n; ++j) 
      u[j + ((head << 1) + 2) * n] -= f[j];

    t = 0.;

    for (int j = 1; j <= n; ++j) {
        // Compute 2nd power
        const double d__1 = u[j + ((head << 1) + 2) * n];
        t += d__1 * d__1;
    }
    t = 1. / sqrt(t);

//  Normalize w_1 and apply same factor to z_1.
    for (int j = 1; j <= n; ++j) {
        u[j + ((head << 1) + 1) * n] = t * u[j + ((head << 1) + 1)*n];
        u[j + ((head << 1) + 2) * n] = t * u[j + ((head << 1) + 2)*n];
    }

//  Update H.
    kptr = link[head - 1];
    for (k = 2; k <= nvec; ++k) {
        h__[k*11 - 11] = 0.;
        for (int j = 1; j <= n; ++j)
            h__[k * 11 - 11] += u[j + ((head << 1) + 2) * n] * u[j + ((kptr << 1) + 2) * n];

        kptr = link[kptr - 1];
    }

//  ==================================================================
//   Compute the Choleski factorization of H.
//  ==================================================================
    k = 2;
    h__[0] = 1.0;
L200:
    if (k > fmin(nvec,*mvec))
        goto L250;

    { // cmp

      for (int j = 1; j <= k - 1; ++j) {
          h__[k + j * 11 - 12] = h__[j + k * 11 - 12];
          for (int i = 1; i <= j - 1; ++i)
              h__[k + j * 11 - 12] -= h__[k + i * 11 - 12] * h__[j + i * 11 - 12];

          h__[k + j*11 - 12] /= h__[j + j*11 - 12];
  // L210:
      }
      h__[k + k*11 - 12] = 1.;
      for (int j = 1; j <= k - 1; ++j) {
          // Compute 2nd power
          const double d__1 = h__[k + j*11 - 12];
          h__[k + k*11 - 12] -= d__1 * d__1;
      }

      const double tol2 = (*tol) * (*tol);
      if (h__[k + k*11 - 12] < tol2) {
  //      -----------------------------------------------
  //       w_k is nearly in span{w_1, ..., w_(k-1)}
  //      -----------------------------------------------
  //      Remove w_k from linked list.
          km1ptr = head;
          for (int j = 2; j <= k - 1; ++j)
            km1ptr = link[km1ptr - 1];

          kptr = link[km1ptr - 1];
          link[km1ptr - 1] = link[kptr - 1];
          --nvec;
  //      Update free storage list.
          link[kptr - 1] = next;
          next = kptr;
  //      Update H.
          for (int j = k; j <= nvec; ++j) {
            for (int i = 1; i <= k - 1; ++i)
                h__[i + j*11 - 12] = h__[i + (j + 1)*11 - 12];

            for (int i = k; i <= j - 1; ++i) 
                h__[i + j*11 - 12] = h__[i + 1 + (j + 1)*11 - 12];
          }

          goto L200;

      } else {
          h__[k + k*11 - 12] = sqrt(h__[k + k*11 - 12]);
          ++k;
          goto L200;
      }
    }

// ---------------------------------------------------------------------------

//  ------------------------------
//   Retain at most MVEC vectors.
//  ------------------------------
L250:
    if (nvec > *mvec) {
//      truncate the linked list.
        last = head;
        for (int j = 2; j <= *mvec; ++j)
            last = link[last - 1];

        tmp  = link[last - 1];
        link[last - 1] = 0;
        last = tmp;
//      Update free storage list.
        for (int j = *mvec + 2; j <= nvec; ++j)
            last = link[last - 1];

        link[last - 1] = next;
        next = tmp;
        nvec = *mvec;
    }
//     ==================================================================
//      Compute the projection of f onto {w_1, ... , w_nvec}.
//     ==================================================================
    jptr = head;
    for (int j = 1; j <= nvec; ++j) {
        c__[j - 1] = 0.0;
        for (int i = 1; i <= n; ++i) // L301:
            c__[j - 1] += f[i] * u[i + ((jptr << 1) + 2)*n];

        jptr = link[jptr - 1];
        for (int i = 1; i <= j - 1; ++i) // L310:
            c__[j - 1] -= h__[j + i * 11 - 12] * c__[i - 1];

        c__[j - 1] /= h__[j + j*11 - 12];
// L300:
    }
    for (int j = nvec; j >= 1; --j) {
        for (int i = j + 1; i <= nvec; ++i) // L321:
            c__[j - 1] -= h__[i + j * 11 - 12] * c__[i - 1];

        c__[j - 1] /= h__[j + j*11 - 12];
// L320:
    }

//  ==================================================================
//   Compute the accelerated correction.
//  ==================================================================
//  Save f for the next call.
    for (int j = 1; j <= n; ++j) // L410:
        u[j + ((next << 1) + 2)*n] = f[j];

    kptr = head;
    for (k = 1; k <= nvec; ++k) {
        for (int j = 1; j <= n; ++j) { /* L421 */
            f[j] = f[j] - c__[k - 1] * u[j + ((kptr << 1) + 2)*n] + 
                    c__[k - 1] * u[j + ((kptr << 1) + 1)*n];
        }
        kptr = link[kptr - 1];
// L420:
    }

//  Save the correction for the next call.
    for (int j = 1; j <= n; ++j) // L430:
        u[j + ((next << 1) + 1) * n] = f[j];

//  ==================================================================
//   Shift the vectors to the right.
//  ==================================================================
    tmp = next;
    next = link[tmp - 1];
    link[tmp - 1] = head;
    head = tmp;
//  Update H.
    for (int j = nvec; j >= 1; --j) {
        for (int i = 1; i <= j - 1; ++i) // L500
            h__[i + 1 + (j + 1)*11 - 12] = h__[i + j*11 - 12];
    }
    ++nvec;
    return 0;
} // naccel_

