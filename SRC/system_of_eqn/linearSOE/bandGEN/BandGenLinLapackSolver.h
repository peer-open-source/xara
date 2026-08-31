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
// Created: Tue Sep 26 16:27:47: 1996
// Revision: A
//
// Description: This file contains the class definition for 
// BandGenLinLapackSolver. It solves the BandGenLinSOE object by calling
// Lapack routines.
//
// What: "@(#) BandGenLinLapackSolver.h, revA"
#pragma once

#include <BandGenLinSolver.h>

class BandGenLinLapackSolver : public BandGenLinSolver
{
  public:
    BandGenLinLapackSolver(bool doDet=true);
    ~BandGenLinLapackSolver();

    int solve();
    int solve(const Vector& B, Vector& X) override;
    int setSize();

    int sendSelf(int commitTag, Channel &);
    int recvSelf(int commitTag, Channel &,  FEM_ObjectBroker &theBroker);

    virtual double getDeterminant() override;
    
  private:
    int *iPiv;
    int iPivSize;
    double det;
    bool doDet;
    void setDeterminant();
};
