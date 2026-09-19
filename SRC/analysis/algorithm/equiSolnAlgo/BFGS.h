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
// Description: This file contains the class definition for BFGS.
//
// See also:
// - FEAP: feap/program/iterat.f
//
// [1] 
//
// Written: Ed Love
// Created: 06/01
//
#pragma once
#include <vector>
#include <EquiSolnAlgo.h>
#include <Vector.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h> 
#include <LinearAction.h>
class LineSearch;

class BFGS: public EquiSolnAlgo
{
  public:

    BFGS(IncrementalIntegrator::TangentFlagType tangent, int n, LineSearch* search=nullptr);

    ~BFGS();

    int solveCurrentStep() final;

    void Print(OPS_Stream &, int flag) const final;    
    

    
  private:
    int BFGSUpdate(IncrementalIntegrator *,
                    LinearSOE *theSOE,
                    Vector &du, 
                    const Vector &b, 
                    int count);


    struct ApplyBFGS //: public LinearAction
    {
      public:
        ApplyBFGS(int n, BFGS *theAlgo) : 
          rdotz(n+3), sdotr(n+3), temp(0)
        {
        }

        int link(IncrementalIntegrator* integrator, LinearSOE* soe) {
          theIntegrator = integrator;
          theSOE = soe;
          temp.resize(soe->getNumEqn());
          return 0;
        }

      public:
        IncrementalIntegrator* theIntegrator;
        LinearSOE* theSOE;
        std::vector<double> rdotz;
        std::vector<double> sdotr;
        Vector temp; // temporary vector 

        Vector **s;  // displacement increments
        Vector **z;
    } action;


    IncrementalIntegrator::TangentFlagType tangent;
    int numberLoops;
    LineSearch* search;

    Vector **s;  // displacement increments
    Vector **z;
    Vector Go;  // residuals
    Vector Gn;

    Vector du; // displacement increment
    Vector b;  // current right-hand side
};
