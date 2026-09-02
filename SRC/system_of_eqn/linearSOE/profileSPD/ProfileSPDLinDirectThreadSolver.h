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
// File: ~/system_of_eqn/linearSOE/profileSPD/ProfileSPDLinDirectThreadSolver.h
//
// Written: fmk 
// Created: February 1997
// Revision: A
//
// Description: This file contains the class definition for 
// ProfileSPDLinDirectThreadSolver. ProfileSPDLinDirectThreadSolver is a subclass 
// of LinearSOESOlver. It solves a ProfileSPDLinSOE object using
// the LDL^t factorization.

// What: "@(#) ProfileSPDLinDirectThreadSolver.h, revA"

#ifndef ProfileSPDLinDirectThreadSolver_h
#define ProfileSPDLinDirectThreadSolver_h

#include <ProfileSPDLinSolver.h>
class  ProfileSPDLinSOE;
struct ProfileTCB;

class ProfileSPDLinDirectThreadSolver : public ProfileSPDLinSolver
{
  public:
    ProfileSPDLinDirectThreadSolver(int numProcessors, int blockSize, double tol);    
    virtual ~ProfileSPDLinDirectThreadSolver();

    int solve() override;        
    int setSize() override;
    bool requireDeterminant() override { return false; }

    virtual int setProfileSOE(ProfileSPDLinSOE &theSOE);

  protected:
    int NP;
    int running;
    
    double minDiagTol;
    int blockSize;
    int maxColHeight;
    int size;
    int *RowTop;
    double **topRowPtr, *invD;
    
  private:
    struct ProfileTCB* tcb;

};

#endif

