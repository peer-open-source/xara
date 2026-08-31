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
#ifndef DistributedProfileSPDLinSOE_h
#define DistributedProfileSPDLinSOE_h

// Written: fmk 
// Description: This file contains the class definition for DistributedProfileSPDLinSOE
// DistributedProfileSPDLinSOE is a subclass of LinearSOE. It uses the LAPACK storage
// scheme to store the components of the A matrix, which is a banded 
// unsymmetric matrix.
//
// What: "@(#) DistributedProfileSPDLinSOE.h, revA"

#include <ProfileSPDLinSOE.h>
#include <Vector.h>

class DistributedProfileSPDLinSolver;

class DistributedProfileSPDLinSOE : public ProfileSPDLinSOE
{
  public:
    DistributedProfileSPDLinSOE(ProfileSPDLinSolver &theSolver);
    DistributedProfileSPDLinSOE();
    
    ~DistributedProfileSPDLinSOE();

    // these methods need to be rewritten
    int addA(const Matrix &, const ID &, double fact = 1.0);
    int addB(const Vector &, const ID &, double fact = 1.0);    
    int setB(const Vector &, double fact = 1.0);            
    void zeroB();
    int setSize(Graph &);
    int solve();
    const Vector &getB();

    int sendSelf(int commitTag, Channel &);
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);

    friend class ProfileSPDLinSolver;    
    friend class ProfileSPDLinDirectSolver;
    friend class ProfileSPDLinDirectBlockSolver;
    friend class ProfileSPDLinDirectThreadSolver;    
    friend class ProfileSPDLinDirectSkypackSolver;    
    friend class ProfileSPDLinSubstrSolver;
    friend class ProfileSPDLinSubstrThreadSolver;

    int setProcessID(int processTag) override;
    int setChannels(int numChannels, Channel **theChannels) override;

  private:
    int processID;
    int numChannels;
    Channel **theChannels;
    ID **localCol;
    ID *sizeLocal;

    double *workArea;
    int sizeWork;
    Vector *myVectB;
    double *myB;
};


#endif

