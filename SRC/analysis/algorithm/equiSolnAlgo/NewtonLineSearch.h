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
// Description: This file contains the class definition for 
// NewtonLineSearch. NewtonLineSearch is a class which performs a Newton-Raphson 
// with line search solution algorithm in solving the equations as outline in
// Crissfields book [1].
//
// Written: fmk 
// Created: 11/96 
// Modified: Ed "C++" Love 10/00 to perform the line search
//
#pragma once
#include <EquiSolnAlgo.h>
#include <LineSearch.h>
#include <Vector.h>

class NewtonLineSearch: public EquiSolnAlgo
{
  public:
    NewtonLineSearch(LineSearch *theLineSearch,
                     IncrementalIntegrator::TangentFlagType prediction_tangent,
                     IncrementalIntegrator::TangentFlagType correction_tangent);
    ~NewtonLineSearch();

    int solveCurrentStep();

    void Print(OPS_Stream &, int flag) const final;    
    
  private:
    LineSearch *theLineSearch;
    IncrementalIntegrator::TangentFlagType
      correction_tangent, 
      prediction_tangent;
    Vector Go, Gn, dX, dXs;
};
