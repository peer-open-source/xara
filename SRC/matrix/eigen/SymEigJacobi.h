#include <VectorND.h>
#include <MatrixND.h>

namespace OpenSees {
template <int n>
int
EigSY3v0(const MatrixND<n,n>&A, VectorND<n>& d, MatrixND<n,n>& v) 
{
  //.... compute eigenvalues and vectors for a 3 x 3 symmetric matrix
  //     using classical cyclic Jacobi plane-rotation iteration
  //
  //.... INPUTS:
  //        A(3,3) - matrix with initial values (only upper half used)
  //
  //.... OUTPUTS
  //        v(3,3) - matrix of eigenvectors (by column)
  //        d(3)   - eigenvalues associated with columns of v
  //        rot    - number of rotations to diagonalize
  //
  //---------------------------------------------------------------eig3==

  //.... Storage done as follows:
  //
  //       | v(1,1) v(1,2) v(1,3) |     |  d(1)  a(1)  a(3)  |
  //       | v(2,1) v(2,2) v(2,3) |  =  |  a(1)  d(2)  a(2)  |
  //       | v(3,1) v(3,2) v(3,3) |     |  a(3)  a(2)  d(3)  |
  //
  //        Transformations performed on d(i) and a(i) and v(i,j) become
  //        the eigenvectors.  
  //
  //---------------------------------------------------------------eig3==
  //
  // Adapted from LovelyEig 
  //
  double  g, h, aij, c, s, tau;

  static constexpr double tol = 1.0e-08;

  v = A;

  //.... move array into one-d arrays
  double a[n];
  a[0] = v(0, 1);
  a[1] = v(1, 2);
  a[2] = v(2, 0);

  double b[n], z[n]{};
  for (int i = 0; i < n; i++) {
    d[i] = v(i, i);
    b[i] = v(i, i);

    for (int j = 0; j < 3; j++)
      v(i, j) = 0.0;

    v(i, i) = 1.0;
  }

  int rot = 0;
  int its = 0;

  double sm = fabs(a[0]) + fabs(a[1]) + fabs(a[2]);

  while (sm > tol) {
    //.... set convergence test and threshold

    double thresh;
    if (its < 3)
      thresh = 0.011 * sm;
    else
      thresh = 0.0;

    //.... perform sweeps for rotations
    for (int i = 0; i < n; i++) {

      int j = (i + 1) % 3;
      int k = (j + 1) % 3;

      aij = a[i];

      g = 100.0 * fabs(aij);

      if (fabs(d(i)) + g != fabs(d(i)) ||
          fabs(d(j)) + g != fabs(d(j))) {

        if (fabs(aij) > thresh) {

          a[i] = 0.0;
          h = d[j] - d[i];

          double t;
          if (std::fabs(h) + g == std::fabs(h))
            t = aij / h;
          else {
            //t = 2.0 * sign(h/aij) / ( fabs(h/aij) + sqrt(4.0+(h*h/aij/aij)));
            double hDIVaij = h / aij;
            if (hDIVaij > 0.0)
              t = 2.0 / (hDIVaij + std::sqrt(4.0 + (hDIVaij * hDIVaij)));
            else
              t = -2.0 / (-hDIVaij + std::sqrt(4.0 + (hDIVaij * hDIVaij)));
          }

          //.... set rotation parameters

          c = 1.0 / std::sqrt(1.0 + t * t);
          s = t * c;
          tau = s / (1.0 + c);

          //.... rotate diagonal terms

          h = t * aij;
          z[i] = z[i] - h;
          z[j] = z[j] + h;
          d[i] = d[i] - h;
          d[j] = d[j] + h;

          //.... rotate off-diagonal terms

          h = a[j];
          g = a[k];
          a[j] = h + s * (g - h * tau);
          a[k] = g - s * (h + g * tau);

          //.... rotate eigenvectors

          for (int k = 0; k < 3; k++) {
            g = v(k, i);
            h = v(k, j);
            v(k, i) = g - s * (h + g * tau);
            v(k, j) = h + s * (g - h * tau);
          }

          rot = rot + 1;

        } // end if fabs > thresh 
      }
      else
        a[i] = 0.0;
    }

    //.... update the diagonal terms
    for (int i = 0; i < 3; i++) {
      b[i] = b[i] + z[i];
      d[i] = b[i];
      z[i] = 0.0;
    }

    its += 1;

    sm = std::fabs(a[0]) + std::fabs(a[1]) + std::fabs(a[2]);
  }

  // sort in descending order (unrolled bubble sort)
  auto sortij = [&d, &v](int i, int j) {
      if (d[i] < d[j]) {
          std::swap(d(i), d(j));
          for (int k = 0; k < 3; ++k)
              std::swap(v(k, i), v(k, j));
      }
  };
  sortij(0, 1);
  sortij(1, 2);
  sortij(0, 1);

  // done
  return 0;
}
} // namespace OpenSees