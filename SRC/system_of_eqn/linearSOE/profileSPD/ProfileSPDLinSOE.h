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
// Revision: A
//
// Description: This file contains the class definition for ProfileSPDLinSOE
// ProfileSPDLinSOE is a subclass of LinearSOE. It uses the LAPACK Upper storage
// scheme to store the components of the A matrix.

// What: "@(#) ProfileSPDLinSOE.h, revA"

#ifndef ProfileSPDLinSOE_h
#define ProfileSPDLinSOE_h

#include <LinearSOE.h>
#include <Vector.h>
class ProfileSPDLinSolver;

class ProfileSPDLinSOE : public LinearSOE
{
  public:
    ProfileSPDLinSOE(ProfileSPDLinSolver &theSolver);
    ProfileSPDLinSOE(ProfileSPDLinSolver &theSolver, int classTag);
    ProfileSPDLinSOE(int classTag);
    ProfileSPDLinSOE(int N, int *iLoc, ProfileSPDLinSolver &theSolver);

    virtual ~ProfileSPDLinSOE();

    int getNumEqn() const override;
    int setSize(Graph &theGraph) override;
    int addA(const Matrix &, const ID &, double fact = 1.0) override;

    virtual int addB(const Vector &, const ID &, double fact = 1.0);    
    virtual int setB(const Vector &, double fact = 1.0);
    
    void zeroA() override;
    void zeroB() override;

    virtual void setX(int loc, double value);
    virtual void setX(const Vector &x);
    
    virtual const Vector &getX();
    virtual const Vector &getB();

    // int formAp(const Vector &p, Vector &Ap) override;

    virtual int setProfileSPDSolver(ProfileSPDLinSolver &newSolver);

    friend class ProfileSPDLinSolver;    
    friend class ProfileSPDLinDirectSolver;
    friend class ProfileSPDLinDirectBlockSolver;
    friend class ProfileSPDLinDirectThreadSolver;    
    friend class ProfileSPDLinDirectSkypackSolver;    
    friend class ProfileSPDLinSubstrSolver;
    friend class ProfileSPDLinSubstrThreadSolver;
    
  protected:
    int size, profileSize;    
    double *A;
    Vector B, X;
    int *iDiagLoc;
    int Asize, Bsize;
    bool isAfactored, isAcondensed;
    int numInt;
    
  private:
};


#endif



