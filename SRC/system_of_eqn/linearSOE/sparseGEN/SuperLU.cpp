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
// Created: Tue 11/96
// Revision: A
//
// Description: This file contains the implementation of SuperLU
//
// What: "@(#) SuperLU.h, revA"

#include <SuperLU.h>
#include <SparseGenColLinSOE.h>
#include <math.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <DataFileStream.h>
#include <iostream>
#include <string>
using std::nothrow;


SuperLU::SuperLU(int perm, 
		 double drop_tolerance, 
		 int panel, 
		 int relx, 
		 char symm)
:SparseGenColLinSolver(SOLVER_TAGS_SuperLU),
 perm_r(0),perm_c(0), etree(0), sizePerm(0),
 relax(relx), permSpec(perm), panelSize(panel), 
 drop_tol(drop_tolerance), symmetric(symm)
{
  // set_default_options(&options);
  options.Fact = DOFACT;
  options.Equil = YES;
  options.ColPerm = COLAMD;
  options.DiagPivotThresh = 1.0;
  options.Trans = NOTRANS;
  options.IterRefine = NOREFINE;
  options.SymmetricMode = NO;
  options.PivotGrowth = NO;
  options.ConditionNumber = NO;
  options.PrintStat = NO;

  if (symmetric == 'Y')
    options.SymmetricMode = YES;

  L.ncol = 0;
  U.ncol = 0;
  A.ncol = 0;
  B.ncol = 0;
  AC.ncol = 0;
}


SuperLU::~SuperLU()
{
  if (perm_r != 0)
    delete [] perm_r;
  if (perm_c != 0)
    delete [] perm_c;
  if (etree != 0) {
    delete [] etree;
    StatFree(&stat);
  }

  if (L.ncol != 0)
    Destroy_SuperNode_Matrix(&L);
  if (U.ncol != 0)
    Destroy_CompCol_Matrix(&U);
  if (AC.ncol != 0) {
    NCPformat *ACstore = (NCPformat *)AC.Store;
    SUPERLU_FREE(ACstore->colbeg);
    SUPERLU_FREE(ACstore->colend);
    SUPERLU_FREE(ACstore);
  }
  if (A.ncol != 0) {
    SUPERLU_FREE(A.Store);
  }
  if (B.ncol != 0) {
    SUPERLU_FREE(B.Store);
  }
}

/*
extern "C" void  dgstrf (char *refact, SuperMatrix* A, double pivThresh, 
			 double dropTol, int relax, 
			 int panelSize, int *etree,
			 void *work, int lwork, int *perm_r, int *perm_c, 
			 SuperMatrix *L, SuperMatrix *U, int *info);

extern "C" void  dgstrs (char *trans, SuperMatrix *L, SuperMatrix *U, 
			 int *perm_r, int *perm_c, 
			 SuperMatrix *B, int *info);    

extern "C" void    StatInit(int, int);

extern "C" void    dCreate_CompCol_Matrix(SuperMatrix *, int, int, int, double *,
					  int *, int *, Stype_t, Dtype_t, Mtype_t);


extern "C" void    get_perm_c(int, SuperMatrix *, int *);

extern "C" void    sp_preorder (char*, SuperMatrix*, int*, int*, SuperMatrix*);

extern "C" void    dCreate_Dense_Matrix(SuperMatrix *, int, int, double *, int,
	       	                     Stype_t, Dtype_t, Mtype_t);

*/

int
SuperLU::solve()
{
    if (theSOE == nullptr) {
      opserr << "WARNING SuperLU::solve(void)- ";
      opserr << " No LinearSOE object has been set\n";
      return -1;
    }
    
    int n = theSOE->size;

    // check for quick return
    if (n == 0)
      return 0;

    if (sizePerm == 0) {
      opserr << "WARNING SuperLU::solve(void)- ";
      opserr << " size for row and col permutations 0 - has setSize() been called?\n";
      return -1;
    }

    /*
    DataFileStream dataStream("K.txt");
    dataStream.open();
    dataStream << n << " " << theSOE->nnz-n << endln;

    // output diagonal entries
    for (int i=0; i<n; i++) 
      for (int j=theSOE->colStartA[i]; j<theSOE->colStartA[i+1]; j++)
	if (theSOE->rowA[j] == i)
	  dataStream << theSOE->A[j] << endln;

    // rowA - diagonals
    for (int i=0; i<=n; i++) 
      dataStream << theSOE->colStartA[i]-i+1 << endln; 

    for (int i=0; i<n; i++) 
      for (int j=theSOE->colStartA[i]; j<theSOE->colStartA[i+1]; j++)
	if (theSOE->rowA[j] != i)
	  dataStream << theSOE->rowA[j]+1 << endln;

    for (int i=0; i<n; i++) 
      for (int j=theSOE->colStartA[i]; j<theSOE->colStartA[i+1]; j++)
	if (theSOE->rowA[j] != i)
	  dataStream << theSOE->A[j] << endln;
    dataStream.close();
    */

    // first copy B into X
    double *Xptr = &theSOE->X[0];
    double *Bptr = &theSOE->B[0];
    for (int i=0; i<n; i++)
      *(Xptr++) = *(Bptr++);

        GlobalLU_t Glu; /* Not needed on return. */

        if (theSOE->factored == false) {
      // factor the matrix
      int info;

      if (L.ncol != 0 && symmetric == 'N') {
        Destroy_SuperNode_Matrix(&L);
        Destroy_CompCol_Matrix(&U);	  
      }

      dgstrf(&options, &AC, relax, panelSize,
            etree, NULL, 0, perm_c, perm_r, &L, &U, &Glu, &stat, &info);


      if (info != 0) {	
        opserr << "WARNING SuperLU::solve(void)- ";
        opserr << " Error " << info << " returned in factorization dgstrf()\n";
        return -info;
      }

      if (symmetric == 'Y')
        options.Fact= SamePattern_SameRowPerm;
      else
        options.Fact = SamePattern;
      
      theSOE->factored = true;
    }	

    // do forward and backward substitution
    trans_t trans = NOTRANS;
    int info;
    dgstrs (trans, &L, &U, perm_c, perm_r, &B, &stat, &info);    

    if (info != 0) {	
       opserr << "WARNING SuperLU::solve(void)- ";
       opserr << " Error " << info << " returned in substitution dgstrs()\n";
       return -info;
    }

    return 0;
}




int
SuperLU::setSize()
{
    int n = theSOE->size;
    if (n > 0) {

      // create space for the permutation vectors 
      // and the elimination tree
      if (sizePerm < n) {

	if (perm_r != 0)
	  delete [] perm_r;
	perm_r = new (nothrow) int[n];		

	if (perm_c != 0)
	  delete [] perm_c;
	perm_c = new (nothrow) int[n];		

	if (etree != 0)
	  delete [] etree;
	etree = new (nothrow) int[n];		

	if (perm_r == 0 || perm_c == 0 || etree == 0) {
	  opserr << "WARNING SuperLU::setSize()";
	  opserr << " - ran out of memory\n";
	  sizePerm = 0;
	  return -1;
	}		
	sizePerm = n;
      }

      // initialisation
      StatInit(&stat);

      // create the SuperMatrix A	
      dCreate_CompCol_Matrix(&A, n, n, theSOE->nnz, theSOE->A, 
			     theSOE->rowA, theSOE->colStartA, 
			     SLU_NC, SLU_D, SLU_GE);

      // obtain and apply column permutation to give SuperMatrix AC
      get_perm_c(permSpec, &A, perm_c);

      sp_preorder(&options, &A, perm_c, etree, &AC);

      // create the rhs SuperMatrix B 
      dCreate_Dense_Matrix(&B, n, 1, &theSOE->X[0], n, SLU_DN, SLU_D, SLU_GE);
	
      // set the refact variable to 'N' after first factorization with new size 
      // can set to 'Y'.
      options.Fact = DOFACT;

      if (symmetric == 'Y')
        options.SymmetricMode=YES;

    } else if (n == 0)
      return 0;
    else {
      opserr << "WARNING SuperLU::setSize()";
      opserr << " - order of system <  0\n";
      return -1;	
    }
	
    return 0;
}

int
SuperLU::sendSelf(int cTag, Channel &theChannel)
{
    // nothing to do
    return 0;
}

int
SuperLU::recvSelf(int ctag,
		  Channel &theChannel, 
		  FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}















