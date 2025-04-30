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
// Description: This file contains the implementation of BandSPDLinLapackSolver.
//
// Written: fmk
// Created: 11/96
//
#include <assert.h>
#include <math.h>
#include <BandSPDLinLapackSolver.h>
#include <BandSPDLinSOE.h>
#include <blasdecl.h>


BandSPDLinLapackSolver::BandSPDLinLapackSolver()
:BandSPDLinSolver(SOLVER_TAGS_BandSPDLinLapackSolver)
{

}

BandSPDLinLapackSolver::~BandSPDLinLapackSolver()
{

}


int
BandSPDLinLapackSolver::solve(void)
{
  assert(theSOE != nullptr);

    int n = theSOE->size;
    int kd = theSOE->half_band -1;
    int ldA = kd +1;
    int nrhs = 1;
    int ldB = n;
    int info;
    double *Aptr = theSOE->A;
    double *Xptr = theSOE->X;
    double *Bptr = theSOE->B;

    // first copy B into X
    for (int i=0; i<n; i++)
      *(Xptr++) = *(Bptr++);
    Xptr = theSOE->X;

    // now solve AX = Y


    char tflag[] = "U";
    if (theSOE->factored == false) {
      // factor and solve
      DPBSV(tflag, &n,&kd,&nrhs,Aptr,&ldA,Xptr,&ldB,&info);

    } else {
      // solve only using factored matrix
      // unsigned int sizeC = 1;
      // DPBTRS("U", sizeC, &n,&kd,&nrhs,Aptr,&ldA,Xptr,&ldB,&info);

        DPBTRS(tflag, &n,&kd,&nrhs,Aptr,&ldA,Xptr,&ldB,&info);
    }


    // check if successful
    if (info != 0) {
      if (info > 0) {
	return -info+1;
      } else {
	return info;
      }
    }

    theSOE->factored = true;
    return 0;
}


int
BandSPDLinLapackSolver::setSize()
{
  // nothing to do
  return 0;
}

int
BandSPDLinLapackSolver::sendSelf(int cTag,
				 Channel &theChannel)
{
  // nothing to do
  return 0;
}

int
BandSPDLinLapackSolver::recvSelf(int tag,
				 Channel &theChannel,
				 FEM_ObjectBroker &theBroker)
{
  // nothing to do
  return 0;
}


