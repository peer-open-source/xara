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
// UmfpackGenLinSolver. It solves the UmfpackGenLinSOEobject by calling
// UMFPACK5.7.1 routines.
//
// Written: fmk 
// Created: 11/98
//
#include <UmfpackGenLinSOE.h>
#include <UmfpackGenLinSolver.h>
#include <math.h>
#include <assert.h>
#include <Constants.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>


UmfpackGenLinSolver::UmfpackGenLinSolver(bool doDet_)
    :LinearSOESolver(SOLVER_TAGS_UmfpackGenLinSolver), 
     Symbolic(nullptr), theSOE(nullptr),
     det(0.0), doDet(doDet_)
{
}


UmfpackGenLinSolver::~UmfpackGenLinSolver()
{
    if (Symbolic != nullptr) {
	umfpack_di_free_symbolic(&Symbolic);
    }
}

double
UmfpackGenLinSolver::getDeterminant()
{
  return doDet? det : OpenSees::Constants::nan;
}

int
UmfpackGenLinSolver::solve()
{
    int n = theSOE->X.Size();
    int nnz = (int)theSOE->Ai.size();
    if (n == 0 || nnz==0) return 0;
    
    int* Ap = &(theSOE->Ap[0]);
    int* Ai = &(theSOE->Ai[0]);
    double* Ax = &(theSOE->Ax[0]);
    double* X = &(theSOE->X(0));
    double* B = &(theSOE->B(0));

    // check if symbolic is done
    assert(Symbolic != 0);
    // if (Symbolic == 0) {
    //     opserr<<"WARNING: setSize has not been called -- Umfpackgenlinsolver::solve\n";
    //     return -1;
    // }
    
    //  perform the numerical factorization
    // numerical analysis
    void* Numeric = nullptr;
    int status = umfpack_di_numeric(Ap,Ai,Ax,Symbolic,&Numeric,Control,Info);

    // check error
    if (status!=UMFPACK_OK) {
      // TODO
      // opserr<<"WARNING: numeric analysis returns "<<status<<" -- Umfpackgenlinsolver::solve\n";
      return -1;
    }

    // solve
    status = umfpack_di_solve(UMFPACK_A,Ap,Ai,Ax,X,B,Numeric,Control,Info);
    
    if (doDet == true)
      umfpack_di_get_determinant(&det, nullptr, Numeric, Info);

    // delete Numeric
    if (Numeric != nullptr) {
      umfpack_di_free_numeric(&Numeric);
    }

    // check error
    if (status != UMFPACK_OK) {
      // opserr<<"WARNING: solving returns "<<status<<" -- Umfpackgenlinsolver::solve\n";
      return -1;
    }

    return 0;
}


int
UmfpackGenLinSolver::setSize()
{
    // set default control parameters
    umfpack_di_defaults(Control);
    Control[UMFPACK_PIVOT_TOLERANCE] = 1.0;
    Control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_SYMMETRIC;

    int n = theSOE->X.Size();
    int nnz = (int)theSOE->Ai.size();
    if (n == 0 || nnz==0) return 0;
    
    int* Ap = &(theSOE->Ap[0]);
    int* Ai = &(theSOE->Ai[0]);
    double* Ax = &(theSOE->Ax[0]);

    // symbolic analysis
    if (Symbolic != nullptr) {
	umfpack_di_free_symbolic(&Symbolic);
    }

    //  perform a column pre-ordering to reduce fill-in
    //  and a symbolic factorization.
    int status = umfpack_di_symbolic(n,n,Ap,Ai,Ax,&Symbolic,Control,Info);

    // check error
    if (status!=UMFPACK_OK) {
	// opserr<<"WARNING: symbolic analysis returns "<<status<<" -- Umfpackgenlinsolver::setsize\n";
	Symbolic = 0;
	return -1;
    }
    return 0;
}

int
UmfpackGenLinSolver::setLinearSOE(UmfpackGenLinSOE &theLinearSOE)
{
    theSOE = &theLinearSOE;
    return 0;
}

int
UmfpackGenLinSolver::sendSelf(int cTag, Channel &theChannel)
{
    // nothing to do
    return 0;
}

int
UmfpackGenLinSolver::recvSelf(int ctag,
			      Channel &theChannel, 
			      FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}

