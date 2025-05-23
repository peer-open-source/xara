
/* 
** Eigen decomposition code for symmetric 3x3 matrices, copied from the public
** domain Java Matrix library JAMA. 
**
** http://barnesc.blogspot.com/2007/02/eigenvectors-of-3x3-symmetric-matrix.html
*/

#include <math.h>

#ifdef MAX
#undef MAX
#endif
#define MAX(a, b) ((a)>(b)?(a):(b))

#define n 3

static inline double 
hypot2(double x, double y) {
    return sqrt(x*x + y*y);
}

// Symmetric Householder reduction to tridiagonal form.

static void 
tred2(double V[n][n], double d[n], double e[n]) {

    //  This is derived from the Algol procedures tred2 by
    //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
    //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
    //  Fortran subroutine in EISPACK.

    int i, j, k;
    for (int j = 0; j < n; j++) {
      d[j] = V[n-1][j];
    }

    // Householder reduction to tridiagonal form.

    for (int i = n-1; i > 0; i--) {

        // Scale to avoid under/overflow.

        double scale = 0.0;
        double h = 0.0;
        for (k = 0; k < i; k++) {
            scale = scale + fabs(d[k]);
        }
        if (scale == 0.0) {
            e[i] = d[i-1];
            for (j = 0; j < i; j++) {
                d[j] = V[i-1][j];
                V[i][j] = 0.0;
                V[j][i] = 0.0;
            }
        }
        else {
            // Generate Householder vector.

            for (k = 0; k < i; k++) {
                d[k] /= scale;
                h += d[k] * d[k];
            }
            double f = d[i-1];
            double g = sqrt(h);
            if (f > 0) {
                g = -g;
            }
            e[i] = scale * g;
            h = h - f * g;
            d[i-1] = f - g;
            for (j = 0; j < i; j++) {
                e[j] = 0.0;
            }

            // Apply similarity transformation to remaining columns.

            for (j = 0; j < i; j++) {
                f = d[j];
                V[j][i] = f;
                g = e[j] + V[j][j] * f;
                for (k = j+1; k <= i-1; k++) {
                    g += V[k][j] * d[k];
                    e[k] += V[k][j] * f;
                }
                e[j] = g;
            }
            f = 0.0;
            for (j = 0; j < i; j++) {
                e[j] /= h;
                f += e[j] * d[j];
            }
            double hh = f / (h + h);
            for (j = 0; j < i; j++) {
                e[j] -= hh * d[j];
            }
            for (j = 0; j < i; j++) {
                f = d[j];
                g = e[j];
                for (k = j; k <= i-1; k++) {
                    V[k][j] -= (f * e[k] + g * d[k]);
                }
                d[j] = V[i-1][j];
                V[i][j] = 0.0;
            }
        }
        d[i] = h;
    }

    // Accumulate transformations.

    for (int i = 0; i < n-1; i++) {
        V[n-1][i] = V[i][i];
        V[i][i] = 1.0;
        double h = d[i+1];
        if (h != 0.0) {
            for (k = 0; k <= i; k++) {
                d[k] = V[k][i+1] / h;
            }
            for (j = 0; j <= i; j++) {
                double g = 0.0;
                for (k = 0; k <= i; k++) {
                    g += V[k][i+1] * V[k][j];
                }
                for (k = 0; k <= i; k++) {
                    V[k][j] -= g * d[k];
                }
            }
        }
        for (k = 0; k <= i; k++) {
            V[k][i+1] = 0.0;
        }
    }
    for (int j = 0; j < n; j++) {
        d[j] = V[n-1][j];
        V[n-1][j] = 0.0;
    }
    V[n-1][n-1] = 1.0;
    e[0] = 0.0;
} 

// Symmetric tridiagonal QL algorithm.

static void
tql2(double V[n][n], double d[n], double e[n]) {

    //  This is derived from the Algol procedures tql2, by
    //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
    //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
    //  Fortran subroutine in EISPACK.
    
    int i, j, k, l;
    for (i = 1; i < n; i++) {
        e[i-1] = e[i];
    }
    e[n-1] = 0.0;

    double f = 0.0;
    double tst1 = 0.0;
    double eps = pow(2.0,-52.0);
    for (l = 0; l < n; l++) {

        // Find small subdiagonal element

        tst1 = MAX(tst1,fabs(d[l]) + fabs(e[l]));
        int m = l;
        while (m < n) {
            if (fabs(e[m]) <= eps*tst1) {
                break;
            }
            m++;
        }

        // If m == l, d[l] is an eigenvalue,
        // otherwise, iterate.

        if (m > l) {
            int iter = 0;
            do {
              iter = iter + 1;  // (Could check iteration count here.)

              // Compute implicit shift

              double g = d[l];
              double p = (d[l+1] - g) / (2.0 * e[l]);
              double r = hypot2(p,1.0);
              if (p < 0) {
                  r = -r;
              }
              d[l] = e[l] / (p + r);
              d[l+1] = e[l] * (p + r);
              double dl1 = d[l+1];
              double h = g - d[l];
              for (i = l+2; i < n; i++) {
                  d[i] -= h;
              }
              f = f + h;

              // Implicit QL transformation.

              p = d[m];
              double c = 1.0;
              double c2 = c;
              double c3 = c;
              double el1 = e[l+1];
              double s = 0.0;
              double s2 = 0.0;
              for (i = m-1; i >= l; i--) {
                  c3 = c2;
                  c2 = c;
                  s2 = s;
                  g = c * e[i];
                  h = c * p;
                  r = hypot2(p,e[i]);
                  e[i+1] = s * r;
                  s = e[i] / r;
                  c = p / r;
                  p = c * d[i] - s * g;
                  d[i+1] = h + s * (c * g + s * d[i]);

                  // Accumulate transformation.

                  for (k = 0; k < n; k++) {
                      h = V[k][i+1];
                      V[k][i+1] = s * V[k][i] + c * h;
                      V[k][i] = c * V[k][i] - s * h;
                  }
              }
              p = -s * s2 * c3 * el1 * e[l] / dl1;
              e[l] = s * p;
              d[l] = c * p;

              // Check for convergence.

            } while (fabs(e[l]) > eps*tst1);
        }
        d[l] = d[l] + f;
        e[l] = 0.0;
    }
  
    // Sort eigenvalues and corresponding vectors.

    for (i = 0; i < n-1; i++) {
        int k = i;
        double p = d[i];
        for (j = i+1; j < n; j++) {
            if (d[j] < p) {
                k = j;
                p = d[j];
            }
        }
        if (k != i) {
            d[k] = d[i];
            d[i] = p;
            for (j = 0; j < n; j++) {
                p = V[j][i];
                V[j][i] = V[j][k];
                V[j][k] = p;
            }
        }
    }
}

void
cmx_eigSY3(double A[n][n], double V[n][n], double d[n]) {

  double e[n];
  for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
          V[i][j] = A[i][j];
      }
  }
  // Tridiagonalize
  tred2(V, d, e);
  // Diagonalize
  tql2(V, d, e);
}

void
cmx_eig3v2(double A[n][n], double EE[n][n], double V[n][n], double d[n])
{
  double e[n];
  double Vn[n][n], Dn[n][n], prod[n][n], U[n]; // Tesser
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      V[i][j] = A[i][j];
    }
  }
  // Tridiagonalize
  tred2(V, d, e);
  // Diagonalize
  tql2(V, d, e);

  // Tesser start
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      V[j][i] = V[j][i] /
                sqrt(V[0][i]*V[0][i] + V[1][i]*V[1][i] + V[2][i]*V[2][i]);
  }

  double supp, supp1; // Tesser
  for (int i = 0; i < n; i++) {
    supp = 0.0;
    for (int j = 0; j < n; j++) {
      supp = supp + V[j][i]*EE[j][0];
    }
    if (supp > -supp) {
      for (int j = 0; j < n; j++)
        Vn[j][i] = V[j][i];
    } else {
      for (int j = 0; j < n; j++)
        Vn[j][i] = -V[j][i];
    }
  }

  supp  = 0.0;
  supp1 = 0.0;
  for (int j = 0; j < n; j++) {
    supp  = supp + EE[j][0]*Vn[j][1];
    supp1 = supp1 + EE[j][0]*Vn[j][0];
  }
  if (supp > supp1) {
    for (int j = 0; j < n; j++) {
      U[j]     = Vn[j][0];
      Vn[j][0] = Vn[j][1];
      Vn[j][1] = U[j];
    }
  }

  supp  = 0.0;
  supp1 = 0.0;
  for (int j = 0; j < n; j++) {
    supp  = supp + EE[j][0]*Vn[j][2];
    supp1 = supp1 + EE[j][0]*Vn[j][0];
  }
  if (supp > supp1) {
    for (int j = 0; j < n; j++) {
      U[j]     = Vn[j][0];
      Vn[j][0] = Vn[j][2];
      Vn[j][2] = U[j];
    }
  }

  for (int i = 1; i < n; i++) {
    supp = 0.0;
    for (int j = 0; j < n; j++) {
      supp = supp + Vn[j][i]*EE[j][1];
    }
    if (supp > -supp) {
      for (int j = 0; j < n; j++)
        Vn[j][i] = Vn[j][i];
    } else {
      for (int j = 0; j < n; j++)
        Vn[j][i] = -Vn[j][i];
    }
  }
  supp  = 0.0;
  supp1 = 0.0;
  for (int j = 0; j < n; j++) {
    supp  = supp + EE[j][1]*Vn[j][2];
    supp1 = supp1 + EE[j][1]*Vn[j][1];
  }
  if (supp > supp1) {
    for (int j = 0; j < n; j++) {
      U[j]     = Vn[j][1];
      Vn[j][1] = Vn[j][2];
      Vn[j][2] = U[j];
    }
  }
  Vn[0][2] = Vn[1][0]*Vn[2][1] - Vn[1][1]*Vn[2][0];
  Vn[1][2] = Vn[0][1]*Vn[2][0] - Vn[0][0]*Vn[2][1];
  Vn[2][2] = Vn[0][0]*Vn[1][1] - Vn[0][1]*Vn[1][0];

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      V[i][j] = Vn[i][j];
    }
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      prod[i][j] = 0.0;
      for (int k = 0; k < n; k++)
        prod[i][j] = prod[i][j] + A[i][k]*Vn[k][j];
    }
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      Dn[i][j] = 0.0;
      for (int k = 0; k < n; k++)
        Dn[i][j] = Dn[i][j] + Vn[k][i]*prod[k][j];
    }
  }
  for (int i = 0; i < n; i++)
    d[i] = Dn[i][i];

  //Tesser end
}
