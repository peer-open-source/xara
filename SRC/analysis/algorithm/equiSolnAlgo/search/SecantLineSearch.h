/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
//
// Description: This file contains the class definition for SecantLineSearch.
// This performs the search for U(i+1) = U(i) + eta * deltaU(i) by using the 
// secant method to find the best solution.
//
//                eta(j+1) = eta(j) -  s(j) * (eta(j-1)-eta(j))
//                                     ------------------------
//                                           s(j-1) - s(j)
//
// where     s(j) = U(i+1,j) ^ R(U(i+1, j))
//
//  and      U(i+1,j) = U(i) + eta(j)*deltaU(i)
// 
// What: "@(#)NewtonLineSearch.h, revA"
//
// Written: fmk 
// Created: 11/01
//
#pragma once
#include <LineSearch.h>
#include <Vector.h>

class SecantLineSearch: public LineSearch
{
  public:
    SecantLineSearch(double tolerance,// = 0.8, 
                     int    maxIter,//   = 10, 
                     double minEta,//    = 0.1, 
                     double maxEta,//    = 10.0, 
                     int    printFlag = 1);

    ~SecantLineSearch();

    int newStep(const Vector &Go) override;
    int search(double s0, 
               double s1, 
               const Vector& dU,
               Vector& G,
               Vector& dXs,
               IncrementalResidual &theIntegrator) override;

    void Print(OPS_Stream &s, int flag) override;
    
  private:
    double tolerance;
    int    maxIter;
    double minEta;
    double maxEta;
    int    printFlag;
};
