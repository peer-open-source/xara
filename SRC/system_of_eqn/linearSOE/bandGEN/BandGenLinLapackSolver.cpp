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
#include <math.h>
#include <assert.h>
#include <blasdecl.h>
#include <BandGenLinLapackSolver.h>
#include <BandGenLinSOE.h>


BandGenLinLapackSolver::BandGenLinLapackSolver(bool doDet_)
:BandGenLinSolver(SOLVER_TAGS_BandGenLinLapackSolver),
 iPiv(0), iPivSize(0), doDet(doDet_)
{

}

BandGenLinLapackSolver::~BandGenLinLapackSolver()
{
  if (iPiv != 0)
     delete [] iPiv;
}

static inline double *
index(double* A, int kl, int ku, int i, int j)
{
  int ldA = 2*kl + ku + 1;
  double *coliiptr = A + j*ldA + kl + ku;
  int diff = j - i;
  if (diff > 0 && diff <= ku) {
      return coliiptr - diff;
  } else {
    diff *= -1;
    if (diff <= kl) {
      return coliiptr + diff;
    }
    else
      return nullptr;
  }
}

void
BandGenLinLapackSolver::setDeterminant()
{
  // A[(ku+1+i-j,j)]
  int kl = theSOE->numSubD;
  int ku = theSOE->numSuperD;

  det = 1.0;
  for (int i=0; i < theSOE->size; i++) {
    det *= *index(theSOE->A, kl, ku, i, i);
  }

  for (int i=0; i < iPivSize; i++)
    if (i+1 != iPiv[i]) {
      det = -det;
    }
}

double
BandGenLinLapackSolver::getDeterminant()
{
  return det;
}


int
BandGenLinLapackSolver::solve(void)
{
    assert(theSOE != nullptr);
    // if (theSOE == nullptr) {
    //     opserr << "WARNING BandGenLinLapackSolver::solve(void)- ";
    //     opserr << " No LinearSOE object has been set\n";
    //     return -1;
    // }

    int n = theSOE->size;
    // check iPiv is large enough
    assert(!(iPivSize < n));
    // if (iPivSize < n) {
    //     opserr << "WARNING BandGenLinLapackSolver::solve(void)- ";
    //     opserr << " iPiv not large enough - has setSize() been called?\n";
    //     return -1;
    // }	

    int kl = theSOE->numSubD;
    int ku = theSOE->numSuperD;
    int ldA = 2*kl + ku +1;
    int nrhs = 1;
    int ldB = n;
    int info;
    double *Aptr = theSOE->A;
    double *Xptr = theSOE->X;
    double *Bptr = theSOE->B;
    int    *iPIV = iPiv;

    // first copy B into X
    for (int i=0; i<n; i++) {
	*(Xptr++) = *(Bptr++);
    }
    Xptr = theSOE->X;

    // now solve AX = B

    char type[] = "N";
    if (theSOE->factored == false)
      // factor and solve
      DGBSV(&n,&kl,&ku,&nrhs,Aptr,&ldA,iPIV,Xptr,&ldB,&info);

    else  {
      // solve only using factored matrix
      // unsigned int sizeC = 1;
      //DGBTRS("N",&sizeC,&n,&kl,&ku,&nrhs,Aptr,&ldA,iPIV,Xptr,&ldB,&info);
      DGBTRS(type,&n,&kl,&ku,&nrhs,Aptr,&ldA,iPIV,Xptr,&ldB,&info);
    }

    // check if successful
    if (info != 0) {
      if (info > 0) {
	// opserr << "WARNING BandGenLinLapackSolver::solve() -";
	// opserr << "factorization failed, matrix singular U(i,i) = 0, i= " << info-1 << endln;
	return -info+1;
      } else {
	// opserr << "WARNING BandGenLinLapackSolver::solve() - OpenSees code error\n";
	return info;
      }
    }

    theSOE->factored = true;
    if (doDet)
      this->setDeterminant();
    return 0;
}



int
BandGenLinLapackSolver::setSize()
{
    // if iPiv not big enough, free it and get one large enough
    if (iPivSize < theSOE->size) {
      if (iPiv != nullptr)
          delete [] iPiv;

      iPiv = new int[theSOE->size];
      iPivSize = theSOE->size;
    }
    return 0;
}

int
BandGenLinLapackSolver::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}

int
BandGenLinLapackSolver::recvSelf(int commitTag,
				 Channel &theChannel,
				 FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}
