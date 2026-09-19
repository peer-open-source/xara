/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
//
// Description: This file contains the class definition for 
// LineSearch. LineSearch is an abstract base class, 
// i.e. no objects of it's type can be created.  Its subclasses seek
// to find a better solution to R(U)=0 than the solution Ui-1 + delta Ui
// would give, typically Ui = Ui-1 + factor * delta Ui.
//
#pragma once
#include <Logging.h>
#include <Vector.h>

class SolutionAlgorithm;
class IncrementalResidual;
class ConvergenceTest;
class LinearSOE;
class OPS_Stream;

class LineSearch
{
  public:
    LineSearch(int classTag);
    virtual ~LineSearch();

    int apply(IncrementalResidual &theIntegrator, 
              LinearSOE &theSOE, 
              ConvergenceTest &theTest,
              Vector &dU, 
              Vector &Go);
    //
    virtual int newStep(const Vector &) =0;

    virtual int search(double so,
                       double su,
                       const Vector& Xo,
                       Vector&       Gus,
                       Vector&       Xs,
                       IncrementalResidual &) =0;

    virtual void Print(OPS_Stream &, int flag) =0;

    void printTrial(int count, double trialEta, double trialPhi) {}

private:
  Vector Gn;
  Vector dXs;
};

