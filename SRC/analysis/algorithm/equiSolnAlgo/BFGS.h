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
// Written: Ed Love
// Created: 06/01
//
#ifndef BFGS_h
#define BFGS_h

#include <EquiSolnAlgo.h>
#include <Vector.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h> 

class BFGS: public EquiSolnAlgo
{
  public:

    BFGS(IncrementalIntegrator::TangentFlagType tangent = CURRENT_TANGENT, int n = 10);
    ~BFGS();

    int solveCurrentStep() final;
    
    virtual int sendSelf(int commitTag, Channel &) final;
    virtual int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &) final;

    void Print(OPS_Stream &, int flag) const final;    
    
  protected:
    
  private:
    void BFGSUpdate(IncrementalIntegrator *theIntegrator,
                    LinearSOE *theSOE,
                    Vector &du, 
                    Vector &b, 
                    int count);

    IncrementalIntegrator::TangentFlagType tangent;
    int numberLoops;

    Vector **s;  // displacement increments
    Vector **z;  
    Vector *residOld;  // residuals
    Vector *residNew;

    Vector *du; //displacement increment

    Vector *b;  //current right-hand side

    Vector *temp; //temporary vector 

    double *rdotz;

    double *sdotr;
  
};

#endif


