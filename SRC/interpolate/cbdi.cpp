//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
#include <math.h>                                                                        
#include <stdlib.h>
#include <Vector.h>
#include <Matrix.h>
#include "cbdi.h"

void
vandermonde(int numSections, const double xi[], Matrix& G)
{
  for (int i = 0; i < numSections; i++) {
    G(i, 0) = 1;
    for (int j = 1; j < numSections; j++)
      G(i, j) = pow(xi[i], j);
  }

  return;
}

void
vandermonde_inverse(int numSections, const double xi[], Matrix& Ginv)
{
  Matrix G(numSections, numSections);
  vandermonde(numSections, xi, G);

  Matrix I(numSections, numSections);
  for (int i = 0; i < numSections; i++)
    I(i, i) = 1.0;

  G.Solve(I, Ginv);
}

void 
getCBDIinfluenceMatrix(int nIntegrPts, const Matrix &xi_pt, double L, Matrix &ls)
{
   // setup Vandermode and CBDI influence matrices
   Matrix G(nIntegrPts, nIntegrPts); 
   Matrix Ginv(nIntegrPts, nIntegrPts);
   Matrix l(nIntegrPts, nIntegrPts);

   for (int i = 1; i <= nIntegrPts; i++)
      for (int j = 1; j <= nIntegrPts; j++) {
         int i0 = i - 1;
         int j0 = j - 1;
         double xi = xi_pt(i0,0);
         G(i0,j0) =  pow(xi,j-1);
         l(i0,j0) = (pow(xi,j+1)-xi)/(j*(j+1));
      }

   G.Invert(Ginv);

   // ls = l * Ginv * (L*L);
   ls.addMatrixProduct(0.0, l, Ginv, L*L);
}

void getCBDIinfluenceMatrix(int nIntegrPts, const double *pts, double L, Matrix &ls)
{
   // setup Vandermode and CBDI influence matrices
   Matrix G(nIntegrPts, nIntegrPts); 
   Matrix Ginv(nIntegrPts, nIntegrPts);
   Matrix l(nIntegrPts, nIntegrPts);

   for (int i = 0; i < nIntegrPts; i++) {
     double xi = pts[i];
     for (int j = 1; j <= nIntegrPts; j++) {
       int j0 = j - 1;
       G(i,j0) =  pow(xi,j-1);
       l(i,j0) = (pow(xi,j+1)-xi)/(j*(j+1));
     }
   }

   G.Invert(Ginv);

   // ls = l * Ginv * (L*L);
   ls.addMatrixProduct(0.0, l, Ginv, L*L);
}

void 
getCBDIinfluenceMatrix(int nPts, const double *pts, int nIntegrPts, const double *integrPts, double L, Matrix &ls)
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
