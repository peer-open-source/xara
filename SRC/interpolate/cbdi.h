//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
#ifndef MatrixUtil_h
#define MatrixUtil_h
#include <MatrixND.h>
using OpenSees::MatrixND;
class Matrix;


static void getCBDIinfluenceMatrix(                             int nIntegrPts, const double *pts,    double L, Matrix &ls);
static void getCBDIinfluenceMatrix(int npts, const double *pts, int nIntegrPts, const double *ipts,   double L, Matrix &ls);


void vandermonde_inverse(int numSections, const double xi[], Matrix& Ginv);


template <int nip, typename T>
static inline void
vandermonde_inverse(int numSections, const double xi[], T& Ginv)
{
  MatrixND<nip,nip> G;
  // setup Vandermonde matrix
  for (int i = 0; i < nip; i++) {
    G(i, 0) = 1;
    for (int j = 1; j < nip; j++)
      G(i, j) = pow(xi[i], j);
  }

//   Matrix I(numSections, numSections);
//   for (int i = 0; i < numSections; i++)
//     I(i, i) = 1.0;

  G.invert(Ginv);
}


static inline void 
getCBDIinfluenceMatrix(int nPts,             const double *pts, 
                       int nIntegrPts, const double *integrPts, 
                       double L, Matrix &ls)
{
   // setup Vandermode and CBDI influence matrices
   double xi;
   Matrix G(nIntegrPts, nIntegrPts); 
   Matrix Ginv(nIntegrPts, nIntegrPts);
   Matrix l(nPts, nIntegrPts);

   // Loop over columns
   for (int j = 1; j <= nIntegrPts; j++) {
     int j0 = j - 1;
     for (int i = 0; i < nIntegrPts; i++) {
       xi = integrPts[i];
       G(i,j0) =  pow(xi,j-1);
     }
     for (int i = 0; i < nPts; i++) {
       xi = pts[i];
       l(i,j0) = (pow(xi,j+1)-xi)/(j*(j+1));
     }
   }

   G.Invert(Ginv);
      
   // ls = l * Ginv * (L*L);
   ls.addMatrixProduct(0.0, l, Ginv, L*L);
}


static inline void
getCBDIinfluenceMatrix(int nIntegrPts, const double *pts, double L, Matrix &ls)
{
  return getCBDIinfluenceMatrix(nIntegrPts, pts, nIntegrPts, pts, L, ls);
}


template <int nPts, int nIntegrPts>
static inline void 
getCBDIinfluenceMatrix(const double *pts, 
                       const double *integrPts, 
                       double L, Matrix &ls)
{
   // setup Vandermode and CBDI influence matrices
   MatrixND<nIntegrPts, nIntegrPts> G{}; 
   MatrixND<nIntegrPts, nIntegrPts> Ginv;
   MatrixND<nPts,nIntegrPts> l{};


   for (int j = 1; j <= nIntegrPts; j++) {
     int j0 = j - 1;
     for (int i = 0; i < nIntegrPts; i++) {
       double xi = integrPts[i];
       G(i,j0) =  pow(xi,j-1);
     }
     for (int i = 0; i < nPts; i++) {
       double xi = pts[i];
       l(i,j0) = (pow(xi,j+1)-xi)/(j*(j+1));
     }
   }

   G.invert(Ginv);
      
   // ls = l * Ginv * (L*L);
   ls.addMatrixProduct(0.0, l, Ginv, L*L);
}

template <int nip>
static inline void
getCBDIinfluenceMatrix(const double *pts, double L, Matrix &ls)
{
  return getCBDIinfluenceMatrix<nip, nip>(pts, pts, L, ls);
}

#endif
