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
// Written: fmk 
// Created: 11/01

// Description: This file contains the class definition for LinearInterpolatedSearch.
// This performs the search by using a form of linear interpolation to find the best solution.
// Solution procedure follows the one in Crissfields book.
// (M.A. Crissfield, Nonlinear Finite Element Analysis of Solid and Structures, Wiley. 97).
// NOTE: it is not quite linear interpolation/false-position/regula-falsi as eta(0) = 0.0
// does not change. Uses eta(i) = eta(i-1)*s0
//                                -----------
//                                s0 - s(i-1)  to compute eta(i)
//
# pragma once
#include <LineSearch.h>
#include <Vector.h>
class OPS_Stream;

class InitialInterpolatedLineSearch: public LineSearch
{
  public:
  InitialInterpolatedLineSearch(
          double tolerance,  // = 0.8, 
          int    maxIter  ,  // = 10,
          double minEta   ,  // = 0.1,
          double maxEta   ,  // = 10.0, 
          int    printFlag); // = 1

    ~InitialInterpolatedLineSearch();

    int newStep(const Vector &) override;
    int search(double s0, 
               double s1, 
               const Vector& dU,
               Vector& G,
               Vector& dUs,
               IncrementalResidual &) override;

    void Print(OPS_Stream &s, int flag) override;
    
  private:
    double tolerance;
    int    maxIter;
    double minEta;
    double maxEta;
    int    printFlag;
};
