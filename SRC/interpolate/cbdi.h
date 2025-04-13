//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
#ifndef MatrixUtil_h
#define MatrixUtil_h

class Matrix;
// void   invertMatrix(int n, const Matrix &a, Matrix &b);
void   getCBDIinfluenceMatrix(int nIntegrPts, const Matrix &xi_pt, double L, Matrix &ls);
void   getCBDIinfluenceMatrix(int nIntegrPts, const double *pts, double L, Matrix &ls);
void   getCBDIinfluenceMatrix(int npts, const double *pts, int nIntegrPts, const double *ipts, double L, Matrix &ls);


void vandermonde(int numSections, const double xi[], Matrix& G);
void vandermonde_inverse(int numSections, const double xi[], Matrix& Ginv);

#endif
