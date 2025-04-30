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
//
// Description: This file contains the class definition for 
// FullGenLinLapackSolver. It solves the FullGenLinSOE object by calling
// Lapack routines.
//
#include <math.h>
#include <assert.h>
#include <blasdecl.h>
#include <FullGenLinLapackSolver.h>
#include <FullGenLinSOE.h>
#include <Matrix.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>


FullGenLinLapackSolver::FullGenLinLapackSolver()
:FullGenLinSolver(SOLVER_TAGS_FullGenLinLapackSolver),iPiv(0),sizeIpiv(0),det(1.0)
{
    
}

FullGenLinLapackSolver::~FullGenLinLapackSolver()
{
    if (iPiv != nullptr)
     delete [] iPiv;
}

void
FullGenLinLapackSolver::setDeterminant()
{

  det = 1.0;
  const Matrix* Amat = theSOE->getA();

  for (int i=0; i < theSOE->size; i++)
    det *= (*Amat)(i, i);

  for (int i=0; i < sizeIpiv; i++)
    if (i+1 != iPiv[i]) {
      det = -det;
    }
}

double
FullGenLinLapackSolver::getDeterminant()
{
  return det;
}


int
FullGenLinLapackSolver::solve(void)
{
    assert(theSOE != nullptr);
    
    int n = theSOE->size;
    
    // check for quick return
    if (n == 0)
      return 0;
    
    // check iPiv is large enough
    assert(!(sizeIpiv < n));
    //  opserr << " iPiv not large enough - has setSize() been called?\n";
     
    int ldA = n;
    int nrhs = 1;
    int ldB = n;
    int info;
    double *Aptr = theSOE->A;
    double *Xptr = theSOE->X;
    double *Bptr = theSOE->B;
    int *iPIV = iPiv;
    
    // first copy B into X
    for (int i=0; i<n; i++)
     *(Xptr++) = *(Bptr++);
    Xptr = theSOE->X;

    //
    // now solve AX = Y
    //
    char tran[] = "N";
    if (theSOE->factored == false)  
     // factor and solve 
      DGESV(&n,&nrhs,Aptr,&ldA,iPIV,Xptr,&ldB,&info);
    else {
     // solve only using factored matrix      
      DGETRS(tran, &n,&nrhs,Aptr,&ldA,iPIV,Xptr,&ldB,&info);      
    }

    
    // check if successful
    if (info != 0) {
      if (info > 0) {
     // opserr << "WARNING FullGenLinLapackSolver::solve() -";
     // opserr << "factorization failed, matrix singular U(i,i) = 0, i= " << info-1 << endln;
     return -info+1;
      } else {
     // opserr << "WARNING FullGenLinLapackSolver::solve() - OpenSees code error\n";
     return info;
      }      
    }

    theSOE->factored = true;

    // we must call setDeterminant while A is factored
    this->setDeterminant();
    return 0;
}


int
FullGenLinLapackSolver::setSize()
{
    int n = theSOE->size;
    if (n > 0) {
     if (sizeIpiv < n) {
         if (iPiv != nullptr)
          delete [] iPiv;
         iPiv = new int[n];          
         sizeIpiv = n;
     }
    } else if (n == 0)
     return 0;
     
    return 0;
}

int
FullGenLinLapackSolver::sendSelf(int commitTag,
                     Channel &theChannel)
{
    // nothing to do
    return 0;
}

int
FullGenLinLapackSolver::recvSelf(int commitTag,
                     Channel &theChannel, 
                     FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}



