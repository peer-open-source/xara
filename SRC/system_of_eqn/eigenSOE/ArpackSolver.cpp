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
// This is the solver that works on the ArpackSOE. It uses the LinearSOE
// in the SOE to perform the solve() operation if required.
// It uses the ARPACK library to perform the eigenvalue analysis.
// ARPACK is an eigen analysis package which was developed by 
// R.B.Lehoucq, D.C.Sorensen and C.Yang at Rice University. ARPACK is a
// collection of FORTRAN77 subroutines designed to solve large scale eigen
// problems. ARPACK is capable of solving large scale non-Hermitian standard 
// and generalized eigen problems. When the matrix <B>K</B> is symmetric, 
// the method is a variant of the Lanczos process called Implicitly Restarted
// Lanczos Method (IRLM).
//
// It is based on previous work of Jun Peng(Stanford)
//
// Written: fmk
// Created: 05.09
//
#include <math.h>
#include <stdio.h>
#include <DataFileStream.h>
#include <ArpackSolver.h>
#include <ArpackSOE.h>
#include <LinearSOE.h>
#include <AnalysisModel.h>
#include <DOF_GrpIter.h>
#include <DOF_Group.h>
#include <FE_EleIter.h>
#include <FE_Element.h>
#include <Integrator.h>
#include <Channel.h>
#include <arpack.hpp>

#include <string.h>
#if !defined(_WIN32)
#  define DSAUPD dsaupd_
#  define DSEUPD dseupd_
#endif

#define ARPACK_ICB

#if !defined(ARPACK_ICB)
extern "C" void DSAUPD(int *ido, char* bmat, 
                       int *n, char *which,
                       int *nev, 
                       double *tol, double *resid, int *ncv, double *v, 
                       int *ldv,
                       int *iparam, int *ipntr, double *workd, double *workl,
                       int *lworkl, int *info);

extern "C" void DSEUPD(bool *rvec, char *howmny,
                       int *select, double *d, double *z,
                       int *ldz, double *sigma, char *bmat,
                       int         *n, char *which,
                       int *nev, double *tol, double *resid, int *ncv, 
                       double *v,
                       int *ldv, int *iparam, int *ipntr, double *workd, 
                       double *workl, int *lworkl, int *info);
#endif 

static double *workArea = nullptr;
static int sizeWork = 0;

ArpackSolver::ArpackSolver()
:EigenSolver(EigenSOLVER_TAGS_ArpackSolver),
 theSOE(0), numModesMax(0), numMode(0), size(0),
 eigenvalues(0), eigenvectors(0), 
 v(0), workl(0), workd(0), resid(0), select(0)
{
  // do nothing here.    
}


ArpackSolver::~ArpackSolver()
{
  if (eigenvalues != 0)
    delete [] eigenvalues;
  if (eigenvectors !=0)
    delete [] eigenvectors;
  if (v != 0)
    delete [] v;
  if (workl != 0)
    delete [] workl;
  if (workd != 0)
    delete [] workd;
  if (resid != 0)
    delete [] resid;
  if (select != 0)
    delete [] select;
  if (workArea != 0)
    delete [] workArea;
  workArea = 0;
  sizeWork = 0;
}


int
ArpackSolver::solve(int numModes, bool generalized, bool findSmallest)
{
  if (generalized == false) {
    opserr << "ArpackSolver::solve() - at moment only solves generalized problem\n";
    return -1;
  }

  theSOE = theArpackSOE->theSOE;

  if (theSOE == nullptr) {
    opserr << "ArpackSolver::setSize() - no LinearSOE set\n";
    return -1;
  }
  
  // set up the space for ARPACK functions.
  // this is done each time method is called!! .. this needs to be cleaned up
  
  int n = size;
  int nev = numModes;
  int ncv = getNCV(n, nev);
  int ldv = n;
  int lworkl = ncv*ncv + 8*ncv;

  int processID = theArpackSOE->processID;

  if (numModes > numModesMax) {
    
    if (v     != nullptr) delete [] v;
    if (workl != nullptr) delete [] workl;
    if (workd != nullptr) delete [] workd;
    if (eigenvalues != nullptr) delete [] eigenvalues;
    if (eigenvectors != nullptr) delete [] eigenvectors;
    if (resid  != nullptr) delete [] resid;
    if (select != nullptr) delete [] select;
    
    v = new double[ldv * ncv];
    workl = new double[lworkl + 1];
    workd = new double[3 * n + 1];
    eigenvalues = new double[nev];
    eigenvectors = new double[n * nev];
    resid = new double[n];
    select = new int[ncv];

    for (int i=0; i<lworkl+1; i++)
      workl[i] = 0;
    for (int i=0; i<3*n+1; i++)
      workd[i] = 0;
    for (int i=0; i<ldv*ncv; i++)
      v[i] = 0;
    
    numModesMax = numModes;
  }

  int ido = 0;
  int ierr = 0;
  double tol = 0.0;
  int info = 0;
  int maxitr = 1000;
  int mode = 3;
  bool rvec = true;
  iparam[0] = 1;
  iparam[2] = maxitr;
  iparam[6] = mode;

#if !defined(ARPACK_ICB)
char which[3];
  if (findSmallest == true)
    strcpy(which, "LM");
  else
    strcpy(which, "SM");

  char bmat  = 'G';
  char howmy = 'A';

#else
  arpack::which which = findSmallest
                      ? arpack::which::largest_magnitude
                      : arpack::which::smallest_magnitude; 

  arpack::bmat bmat = arpack::bmat::generalized;
  arpack::howmny howmy = arpack::howmny::ritz_vectors;
#endif

  while (1) { 
      
#if  !defined(ARPACK_ICB)
    DSAUPD(&ido, &bmat, &n, which, &nev, &tol, resid, 
           &ncv, v, &ldv,
           iparam, ipntr, workd, workl, &lworkl, &info);
#else
    arpack::saupd(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv,
                  iparam, ipntr, workd, workl, lworkl, info);
#endif

    if (theArpackSOE->checkSameInt(ido) != 1) {
      opserr << "ArpackSolver::solve - ido values not the same .. ido, processID: "
             << ido << " " << processID << endln;
      return -1;
    }

    if (ido == -1) {
      // IDO = -1: compute  Y = OP * X  where
      //           IPNTR(1) is the pointer into WORKD for X,
      //           IPNTR(2) is the pointer into WORKD for Y.
      //           This is for the initialization phase to force the
      //           starting vector into the range of OP.
    
      myMv(n, &workd[ipntr[0]-1], &workd[ipntr[1]-1]); 
    
      theVector.setData(&workd[ipntr[1] - 1], size);
     
      if (processID > 0)
        theSOE->zeroB();
      else
        theSOE->setB(theVector);

      ierr = theSOE->solve();
      const Vector &X = theSOE->getX();
      theVector = X;

      continue; 
    }

    else if (ido == 1) {
      // compute  Y = OP * X  where
      // IPNTR(1) is the pointer into WORKD for X,
      // IPNTR(2) is the pointer into WORKD for Y.
      // In mode 3, the vector B * X is already
      // available in WORKD(ipntr(3)).  It does not
      // need to be recomputed in forming OP * X.

      // double ratio = 1.0;
      myCopy(n, &workd[ipntr[2]-1], &workd[ipntr[1]-1]);
      theVector.setData(&workd[ipntr[1] - 1], size);
      if (processID > 0)
        theSOE->zeroB();
      else
        theSOE->setB(theVector);

      theSOE->solve();
   
      const Vector &X = theSOE->getX();
      theVector = X;
      // theVector.setData(&workd[ipntr[1] - 1], size);

      continue;
    }

    else if (ido == 2) {
      // IDO =  2: compute  Y = M * X  where
      //           IPNTR(1) is the pointer into WORKD for X,
      //           IPNTR(2) is the pointer into WORKD for Y.
      myMv(n, &workd[ipntr[0]-1], &workd[ipntr[1]-1]);
      continue;
    }
    break;
  }
  
  if (info < 0) {
    opserr << "Arpack returned with flag " << info << "\n  ";
    switch(info) {
    case -1: 
      opserr << "N must be positive.\n";
      break;
    case -2: 
      opserr << "NEV must be positive.\n";
      break;
    case -3: 
      opserr << "NCV must be greater than NEV and less than or equal to N.\n";
      break;
    case -4:
      opserr << "The maximum number of Arnoldi update iterations allowed\n";
      break;
    case -5: 
      opserr << "WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.\n";
      break;
    case -6: 
      opserr << "BMAT must be one of 'I' or 'G'.\n";
      break;
    case -7: 
      opserr << "Length of private work array WORKL is not sufficient.\n";
      break;
    case -8: 
      opserr << "Error return from trid. eigenvalue calculation\n";
      opserr << "Informatinal error from LAPACK routine dsteqr.\n";
      break;
    case -9: 
      opserr << "Starting vector is zero.\n";
      break;
    case -10: 
      opserr << "IPARAM(7) must be 1,2,3,4,5.\n";
      break;
    case -11: 
      opserr << "IPARAM(7) = 1 and BMAT = 'G' are incompatible.\n";
      break;
    case -12: 
      opserr << "IPARAM(1) must be equal to 0 or 1.\n";
      break;
    case -13:
      opserr << "NEV and WHICH = 'BE' are incompatible.\n";
      break;
    case -9999:
      opserr << "Could not build an Arnoldi factorization.\n";
      opserr << "IPARAM(5) - the size of the current Arnoldi factorization is " << iparam[4] << ".\n";
      opserr << "Check that enough workspace and array storage have been allocated.\n";
      break;
    default:
      opserr << "unrecognised return value\n";
    }

    if (eigenvalues != nullptr) 
      delete [] eigenvalues;
    eigenvalues = nullptr;
    if (eigenvectors != nullptr)
      delete [] eigenvectors;
    eigenvectors = nullptr;
    
    return info;

  } else {
    if (info == 1) {
      opserr << "ArpackSolver::Maximum number of iteration reached.\n";
    }
    else if (info == 3) {
      opserr << "ArpackSolver::No Shifts could be applied during implicit,";
      opserr << "Arnoldi update, try increasing NCV." << endln;
    }
    
    double sigma = shift;
    if (iparam[4] > 0) {
      rvec = true;
      n = size;    
      ldv = n;

#if !defined(ARPACK_ICB)
      DSEUPD(&rvec, &howmy, (int*)(select), 
             eigenvalues, eigenvectors, &ldv, &sigma, &bmat, &n, which,
             &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd,
             workl, &lworkl, &info);
#else
      arpack::seupd(rvec, howmy, select,
                    eigenvalues, eigenvectors, 
                    ldv, sigma, bmat, n, which,
                    nev, tol, resid, ncv, v, ldv, iparam, ipntr, workd,
                    workl, lworkl, info);
#endif
      if (info != 0) {
        opserr << "ArpackSolver::Error with dseupd_" << info;
        switch(info) {
          
        case -1: 
          opserr << " N must be positive.\n";
          break;
        case -2: 
          opserr << " NEV must be positive.\n";
          break;
        case -3: 
          opserr << " NCV must be greater than NEV and less than or equal to N.\n";
          break;
        case -5: 
          opserr << " WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'.\n";
          break;
        case -6: 
          opserr << " BMAT must be one of 'I' or 'G'.\n";
          break;
        case -7: 
          opserr << " Length of private work WORKL array is not sufficient.\n";
          break;
        case -8: 
          opserr << " Error return from trid. eigenvalue calculation";
          opserr << "Information error from LAPACK routine dsteqr.\n";
          break;
        case -9: 
          opserr << " Starting vector is zero.\n";
          break;
        case -10: 
          opserr << " IPARAM(7) must be 1,2,3,4,5.\n";
          break;
        case -11: 
          opserr << " IPARAM(7) = 1 and BMAT = 'G' are incompatibl\n";
          break;
        case -12: 
          opserr << " NEV and WHICH = 'BE' are incompatible.\n";
          break;
        case -14: 
          opserr << " DSAUPD did not find any eigenvalues to sufficient accuracy.\n";
          break;
        case -15: 
          opserr << " HOWMNY must be one of 'A' or 'S' if RVEC = .true.\n";
          break;
        case -16: 
          opserr << " HOWMNY = 'S' not yet implemented\n";
          break;
        default:
          ;
        }

        return info;
        
      }
    }
  }
  
  numMode = numModes;

  return 0;
}


int 
ArpackSolver::getNCV(int n, int nev)
{
  int result;
  if (2*nev > nev+8) {
    result = nev+8;
  } else {
    result = 2*nev;
  }
  
  if (result >= n) {
    result = n;
  }
  
  return result;
}


void
ArpackSolver::myMv(int n, double *v, double *result)
{
  Vector x(v, n);
  Vector y(result,n);
    
  bool mDiagonal = theArpackSOE->mDiagonal;

  if (mDiagonal == true) {

    int Msize = theArpackSOE->Msize;
    double *M = theArpackSOE->M;

    if (n <= Msize) {
      for (int i=0; i<n; i++)
        result[i] = M[i]*v[i];
    } else {
      opserr << "ArpackSolver: n > Msize!\n";
      return;
    }

  } else {

    y.Zero();

    AnalysisModel *theAnalysisModel = theArpackSOE->theModel;
    
    // loop over the FE_Elements
    FE_Element *elePtr;
    FE_EleIter &theEles = theAnalysisModel->getFEs();    
    while((elePtr = theEles()) != 0) {
      const Vector &b = elePtr->getM_Force(x, 1.0);
      y.Assemble(b, elePtr->getID(), 1.0);
    }

    // loop over the DOF_Groups
    DOF_Group *dofPtr;
    DOF_GrpIter &theDofs = theAnalysisModel->getDOFs();
    while ((dofPtr = theDofs()) != 0) {
      const Vector &a = dofPtr->getM_Force(x,1.0);      
      y.Assemble(a, dofPtr->getID(), 1.0);
    }
  }

  // if parallel we have to merge the results
  int processID = theArpackSOE->processID;
  if (processID != -1) {
    Channel **theChannels = theArpackSOE->theChannels;
    int numChannels = theArpackSOE->numChannels;
    if (processID != 0) {
      theChannels[0]->sendVector(0, 0, y);
      theChannels[0]->recvVector(0, 0, y);
    } else {
      Vector other(workArea, n);
      // recv contribution from remote & add
      for (int i=0; i<numChannels; i++) {
        theChannels[i]->recvVector(0,0,other);
        y += other;
      }
      // send result back
      for (int i=0; i<numChannels; i++) {
        theChannels[i]->sendVector(0,0,y);
      }
    }
  }
}
    
void
ArpackSolver::myCopy(int n, double *v, double *result)
{
  for (int i=0; i<n; i++) {
    result[i] = v[i];
  }
}


int
ArpackSolver::setEigenSOE(ArpackSOE &theArpSOE)
{
  theArpackSOE = &theArpSOE;
  shift = theArpackSOE->getShift();
  return 0;
}


const Vector &
ArpackSolver::getEigenvector(int mode)
{
  if (mode <= 0 || mode > numMode) {
    theVector.Zero();
  }
  
  int index = (mode - 1) * size;
  
  theVector.setData(&eigenvectors[index], size);

  return theVector;;  
}


double
ArpackSolver::getEigenvalue(int mode)
{
  if (mode <= 0 || mode > numMode) {
    opserr << "ArpackSOE::getEigenvalue() - mode is out of range(1 - nev)";
    return -1;
  }
  
  if (eigenvalues != 0)
    return eigenvalues[mode-1];
  else {
    opserr << "ArpackSOE::getEigenvalue() - eigenvalues not yet determined";
    return -2;
  }      
}


int
ArpackSolver::setSize()
{
  size = theArpackSOE->Msize;

  if (sizeWork < size)
    if (workArea != 0)
      delete [] workArea;

  workArea = new double[size];
  
  return 0;
}


int    
ArpackSolver::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}

int
ArpackSolver::recvSelf(int commitTag, Channel &theChannel, 
                       FEM_ObjectBroker &theBroker)
{
  return 0;
}


