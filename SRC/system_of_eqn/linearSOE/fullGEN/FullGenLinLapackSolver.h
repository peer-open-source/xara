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
// FullGenLinLapackSolver. It solves the FullGenLinSOE object by calling
// Lapack routines.
//
// Written: fmk 
// Created: Tue Sep 26 16:27:47: 1996
// Revision: A
//
#ifndef FullGenLinLapackSolver_h
#define FullGenLinLapackSolver_h

#include "FullGenLinSolver.h"

class FullGenLinLapackSolver : public FullGenLinSolver
{
  public:
    FullGenLinLapackSolver();    
    ~FullGenLinLapackSolver();

    int solve();
    int solve(const Vector& B, Vector& X) override;
    int setSize();
    
    virtual int setLinearSOE(FullGenLinSOE &theSOE);

    virtual double getDeterminant() override;

  private:
    int *iPiv;
    int sizeIpiv;
    double det;
    FullGenLinSOE *theSOE;
    void setDeterminant();
};

#endif

