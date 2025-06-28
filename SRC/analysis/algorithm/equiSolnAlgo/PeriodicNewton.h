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
#ifndef PeriodicNewton_h
#define PeriodicNewton_h
// Description: This file contains the class definition for 
// PeriodicNewton. PeriodicNewton is a class which performs a Periodic 
// Newton-Raphson  solution algorithm in solving the equations.
// No member functions are declared as virtual as 
// it is not expected that this class will be subclassed.
//
// Written: MHS
// Created: Oct 2002
//
#include <EquiSolnAlgo.h>
#include <Vector.h>

class ConvergenceTest;

class PeriodicNewton: public EquiSolnAlgo
{
  public:
    PeriodicNewton(int tangent = CURRENT_TANGENT, int maxCount = 3);
    PeriodicNewton(ConvergenceTest &theTest, int tangent = CURRENT_TANGENT, int maxCount = 3);
    ~PeriodicNewton();

    int solveCurrentStep(void);    

    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &, int flag) final;    
    
  protected:
    
  private:
    int tangent;
    int maxCount;
};

#endif


