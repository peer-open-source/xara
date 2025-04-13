/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** ****************************************************************** */
//
// This file contains the definition for NURBS derivatives
// Written originally by Vinh Phu Nguyen, nvinhphu@gmail.com
//
#ifndef NurbsDers_h
#define NurbsDers_h

class Vector;
class Matrix;

int      FindSpan(int n, int p, double u, Vector& U);
void     BasisFuns( int i, double u, int p, Vector& U, Vector& N);
void     dersBasisFuns(int i, double u, int p, int order, Vector& knot, Matrix& ders);
double   OneBasisFun(int p, int m, Vector U, int i, double u);
void     dersOneBasisFuns(int p, int m, Vector U, int i, double u, int n, double* ders);

#endif
