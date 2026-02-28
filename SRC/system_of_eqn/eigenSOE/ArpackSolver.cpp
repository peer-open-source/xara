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
**                                                                     **
** ****************************************************************** */
//
// This is the solver that works on the ArpackSOE. It uses the LinearSOE
// in the SOE to perform the solve() operation if required.
// It uses the ARPACK library to perform the eigenvalue analysis.
// ARPACK is an eigen analysis package which was developed by 
// R.B.Lehoucq, D.C.Sorensen and C.Yang at Rice University. 
// ARPACK is a collection of FORTRAN77 subroutines designed to solve large scale eigen
// problems. ARPACK is capable of solving large scale non-Hermitian standard 
// and generalized eigen problems. When the matrix <B>K</B> is symmetric, 
// the method is a variant of the Lanczos process called Implicitly Restarted
// Lanczos Method (IRLM).
//
// It is based on previous work of Jun Peng(Stanford)
//
// From SciPy (Symmetric):
//
// # The following modes are supported:
// #  mode = 1:
// #    Solve the standard eigenvalue problem:
// #      A*x = lambda*x :
// #       A - symmetric
// #    Arguments should be
// #       matvec      = left multiplication by A
// #       M_matvec    = None [not used]
// #       Minv_matvec = None [not used]
// #
// #  mode = 2:
// #    Solve the general eigenvalue problem:
// #      A*x = lambda*M*x
// #       A - symmetric
// #       M - symmetric positive definite
// #    Arguments should be
// #       matvec      = left multiplication by A
// #       M_matvec    = left multiplication by M
// #       Minv_matvec = left multiplication by M^-1
// #
// #  mode = 3:
// #    Solve the general eigenvalue problem in shift-invert mode:
// #      A*x = lambda*M*x
// #       A - symmetric
// #       M - symmetric positive semi-definite
// #    Arguments should be
// #       matvec      = None [not used]
// #       M_matvec    = left multiplication by M
// #                     or None, if M is the identity
// #       Minv_matvec = left multiplication by [A-sigma*M]^-1
// #
// #  mode = 4:
// #    Solve the general eigenvalue problem in Buckling mode:
// #      A*x = lambda*AG*x
// #       A  - symmetric positive semi-definite
// #       AG - symmetric indefinite
// #    Arguments should be
// #       matvec      = left multiplication by A
// #       M_matvec    = None [not used]
// #       Minv_matvec = left multiplication by [A-sigma*AG]^-1
// #
// #  mode = 5:
// #    Solve the general eigenvalue problem in Cayley-transformed mode:
// #      A*x = lambda*M*x
// #       A - symmetric
// #       M - symmetric positive semi-definite
// #    Arguments should be
// #       matvec      = left multiplication by A
// #       M_matvec    = left multiplication by M
// #                     or None, if M is the identity
// #       Minv_matvec = left multiplication by [A-sigma*M]^-1
//
// For unsymmetric:
// # The following modes are supported:
// #  mode = 1:
// #    Solve the standard eigenvalue problem:
// #      A*x = lambda*x
// #       A - square matrix
// #    Arguments should be
// #       matvec      = left multiplication by A
// #       M_matvec    = None [not used]
// #       Minv_matvec = None [not used]
// #
// #  mode = 2:
// #    Solve the generalized eigenvalue problem:
// #      A*x = lambda*M*x
// #       A - square matrix
// #       M - symmetric, positive semi-definite
// #    Arguments should be
// #       matvec      = left multiplication by A
// #       M_matvec    = left multiplication by M
// #       Minv_matvec = left multiplication by M^-1
// #
// #  mode = 3,4:
// #    Solve the general eigenvalue problem in shift-invert mode:
// #      A*x = lambda*M*x
// #       A - square matrix
// #       M - symmetric, positive semi-definite
// #    Arguments should be
// #       matvec      = None [not used]
// #       M_matvec    = left multiplication by M
// #                     or None, if M is the identity
// #       Minv_matvec = left multiplication by [A-sigma*M]^-1
// #    if A is real and mode==3, use the real part of Minv_matvec
// #    if A is real and mode==4, use the imag part of Minv_matvec
// #    if A is complex and mode==3,
// #       use real and imag parts of Minv_matvec
//
// Written: fmk
// Created: 05.09
//
#include <cmath>
#include <vector>
#include <iostream>
#include <random>
#include <algorithm>
#include <cstring> // memcpy
#include <limits>
#include <ArpackSolver.h>
#include <ArpackSOE.h>
#include <LinearSOE.h>
#include <Channel.h>
#include <Logging.h>

#define ARPACK_ICB

#if defined(ARPACK_ICB)
# include <arpack.hpp>
#else
#endif 


namespace {

typedef struct { int code; const char *msg; } ErrorEntry;

static const ErrorEntry DnaupdErrors[] = {
  { 0,     "Normal exit." },
  { 1,     "Maximum number of iterations taken. All possible eigenvalues of OP has been found. IPARAM(5) returns the number of wanted converged Ritz values." },
  { 2,     "No longer an informational error. Deprecated starting with release 2 of ARPACK." },
  { 3,     "No shifts could be applied during a cycle of the Implicitly restarted Arnoldi iteration. One possibility is to increase the size of NCV relative to NEV." },
  { -1,    "N must be positive." },
  { -2,    "NEV must be positive." },
  { -3,    "NCV-NEV >= 2 and less than or equal to N." },
  { -4,    "The maximum number of Arnoldi update iterations allowed must be greater than zero." },
  { -5,    "WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'" },
  { -6,    "BMAT must be one of 'I' or 'G'." },
  { -7,    "Length of private work array WORKL is not sufficient." },
  { -8,    "Error return from LAPACK eigenvalue calculation;" },
  { -9,    "Starting vector is zero." },
  { -10,   "IPARAM(7) must be 1,2,3,4." },
  { -11,   "IPARAM(7) = 1 and BMAT = 'G' are incompatible." },
  { -12,   "IPARAM(1) must be equal to 0 or 1." },
  { -13,   "NEV and WHICH = 'BE' are incompatible." },
  { -9999, "Could not build an Arnoldi factorization. IPARAM(5) returns the size of the current Arnoldi factorization. The user is advised to check that enough workspace and array storage has been allocated." }
};
// static const size_t DnaupdErrorsCount = sizeof(DnaupdErrors)/sizeof(DnaupdErrors[0]);


static const ErrorEntry DneupdErrors[] = {
  { 0,    "Normal exit." },
  { 1,    "The Schur form computed by LAPACK routine dlahqr could not be reordered by LAPACK routine dtrsen. Re-enter subroutine dneupd with IPARAM(5)NCV and increase the size of the arrays DR and DI to have dimension at least dimension NCV and allocate at least NCV columns for Z. NOTE: Not necessary if Z and V share the same space. Please notify the authors if this erroroccurs." },
  { -1,   "N must be positive." },
  { -2,   "NEV must be positive." },
  { -3,   "NCV-NEV >= 2 and less than or equal to N." },
  { -5,   "WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'" },
  { -6,   "BMAT must be one of 'I' or 'G'." },
  { -7,   "Length of private work WORKL array is not sufficient." },
  { -8,   "Error return from calculation of a real Schur form. Informational error from LAPACK routine dlahqr ." },
  { -9,   "Error return from calculation of eigenvectors. Informational error from LAPACK routine dtrevc." },
  { -10,  "IPARAM(7) must be 1,2,3,4." },
  { -11,  "IPARAM(7) = 1 and BMAT = 'G' are incompatible." },
  { -12,  "HOWMNY = 'S' not yet implemented" },
  { -13,  "HOWMNY must be one of 'A' or 'P' if RVEC = .true." },
  { -14,  "DNAUPD  did not find any eigenvalues to sufficient accuracy." },
  { -15,  "DNEUPD got a different count of the number of converged Ritz values than DNAUPD got.  This indicates the user probably made an error in passing data from DNAUPD to DNEUPD or that the data was modified before entering DNEUPD" }
};

static const ErrorEntry DsaupdErrors[] = {
  { 0,     "Normal exit." },
  { 1,     "Maximum number of iterations taken. All possible eigenvalues of OP has been found." },
  { 2,     "No longer an informational error. Deprecated starting with release 2 of ARPACK." },
  { 3,     "No shifts could be applied during a cycle of the Implicitly restarted Arnoldi iteration. One possibility is to increase the size of NCV relative to NEV." },
  { -1,    "N must be positive." },
  { -2,    "NEV must be positive." },
  { -3,    "NCV must be greater than NEV and less than or equal to N." },
  { -4,    "The maximum number of Arnoldi update iterations allowed must be greater than zero." },
  { -5,    "WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'." },
  { -6,    "BMAT must be one of 'I' or 'G'." },
  { -7,    "Length of private work array WORKL is not sufficient." },
  { -8,    "Error return from trid. eigenvalue calculation; Informational error from LAPACK routine dsteqr ." },
  { -9,    "Starting vector is zero." },
  { -10,   "IPARAM(7) must be 1,2,3,4,5." },
  { -11,   "IPARAM(7) = 1 and BMAT = 'G' are incompatible." },
  { -12,   "IPARAM(1) must be equal to 0 or 1." },
  { -13,   "NEV and WHICH = 'BE' are incompatible." },
  { -9999, "Could not build an Arnoldi factorization. IPARAM(5) returns the size of the current Arnoldi factorization. The user is advised to check that enough workspace and array storage has been allocated." }
};

static const ErrorEntry DseupdErrors[] = {
  { 0,    "Normal exit." },
  { -1,   "N must be positive." },
  { -2,   "NEV must be positive." },
  { -3,   "NCV must be greater than NEV and less than or equal to N." },
  { -5,   "WHICH must be one of 'LM', 'SM', 'LA', 'SA' or 'BE'." },
  { -6,   "BMAT must be one of 'I' or 'G'." },
  { -7,   "Length of private work WORKL array is not sufficient." },
  { -8,   "Error return from trid. eigenvalue calculation; Information error from LAPACK routine dsteqr." },
  { -9,   "Starting vector is zero." },
  { -10,  "IPARAM(7) must be 1,2,3,4,5." },
  { -11,  "IPARAM(7) = 1 and BMAT = 'G' are incompatible." },
  { -12,  "NEV and WHICH = 'BE' are incompatible." },
  { -14,  "DSAUPD  did not find any eigenvalues to sufficient accuracy." },
  { -15,  "HOWMNY must be one of 'A' or 'S' if RVEC = .true." },
  { -16,  "HOWMNY = 'S' not yet implemented" },
  { -17,  "DSEUPD  got a different count of the number of converged Ritz values than DSAUPD  got.  This indicates the user probably made an error in passing data from DSAUPD  to DSEUPD  or that the data was modified before entering  DSEUPD." }
};


static inline const char* 
LookupArpackError(const ErrorEntry *table, size_t count, int code) {
  for (size_t i = 0; i < count; ++i)
    if (table[i].code == code) return table[i].msg;
  return "Unknown error code.";
}

struct ArpackWorkspace {
  ArpackWorkspace(int n, int nev, int symm) :
    size(n),
    ldv(n),
    nev(nev),
    ncv(getNCV(n, nev, symm)),
    lworkl(symm==ArpackWorkspace::Symmetric ? 1*(  ncv*ncv + 8*ncv) 
                                            : 1*(3*ncv*ncv + 6*ncv)),
    v(new double[n * ncv]{}),
    workl(new double[lworkl+1]{}),
    workd(new double[3*n]{}),
    resid(new double[n]{}),
    select(new int[ncv]{})
  {
    assert(n >=0);
    assert(ncv <= n);
    if (symm == ArpackWorkspace::NonSymmetric) {
      assert(ncv > nev + 2);
    }
    else {
      assert(ncv > nev);
    }
  }

  ~ArpackWorkspace() {
    delete [] v;
    delete [] workl;
    delete [] workd;
    delete [] resid;
    delete [] select;
  }

  ArpackWorkspace(const ArpackWorkspace&) = delete;
  ArpackWorkspace& operator=(const ArpackWorkspace&) = delete;
  ArpackWorkspace(ArpackWorkspace&&) = delete;
  ArpackWorkspace& operator=(ArpackWorkspace&&) = delete;

  enum { Symmetric = 1, NonSymmetric = 2 };
  int getNCV(int n, int nev, int driver)
  {
    // compute the number of Arnoldi vectors to use
    // n is the system size, nev is the number of eigenvectors.
    // dsaupd: NCV must be greater than NEV and less than or equal to N.

    // For dnaupd: NCV must satisfy the two inequalities 
    //.       2 <= NCV-NEV and NCV <= N.
    // The only formal requrement is that NCV > NEV + 2.
    // However, it is recommended that NCV .ge. 2*NEV+1. See Chapter 8 of:
    // 
    // 2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
    //    Restarted Arnoldi Iteration", Rice University Technical Report
    //    TR95-13, Department of Computational and Applied Mathematics.
    //
    // Scipy uses min(max(2 * n + 1, 20), n)
  #if 1
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
  #elif 1
    return std::min(std::max(2*nev + 1, 20), n);
  #else
    // ensure headroom and a sensible floor
    int ncv = std::max({2*nev + 8, nev + 20, 20});
    ncv = std::min(ncv, n);
    // ARPACK requires ncv > nev+1; enforce with a small bump if saturated
    if (ncv <= nev + 1 && n > nev + 1) 
      ncv = std::min(n, nev + 2);
    return ncv;
  #endif
  }

  int size; // ldv
  int ldv;
  int nev;
  int ncv;
  int lworkl;
  double* v;
  double* workl;
  double* workd;
  double* resid;
  int*    select;
};

} // namespace


ArpackSolver::ArpackSolver()
:EigenSolver(EigenSOLVER_TAGS_ArpackSolver),
 theSOE(nullptr), numMode(0), size(0), iparam{0}, ipntr{0}, shift(0.0)
{

}


ArpackSolver::~ArpackSolver()
{
  solution.clear();
}


int
ArpackSolver::solve(int numModes, bool generalized, bool findSmallest)
{
  if (generalized == false) {
    return solveI(numModes, generalized, findSmallest);
    // opserr << "ArpackSolver::solve - at moment only solves generalized problem\n";
    return -1;
  }

  theSOE = theArpackSOE->theSOE;
  numMode = 0;

  if (theSOE == nullptr) {
    opserr << "ArpackSolver::setSize() - no LinearSOE set\n";
    return -1;
  }
  
  // set up the space for ARPACK functions.
  //
  // this is done each time method is called!! .. this needs to be cleaned up
  
  int n = theArpackSOE->getNumEqn(); // size;
  const int nev = numModes;
  if (n < nev || nev < 1) {
    opserr << "ArpackSolver::solve - no. of modes requested is invalid\n";
    return -1;
  }
  solution.reserve(n, numModes);

  int ido = 0;
  double tol = std::numeric_limits<double>::epsilon();
  int info = 0;
  int maxitr = 1000;
  int mode = 3; //
  iparam[0] = 1; // exact shifts
  iparam[1] = 0; // not used by ARPACK
  iparam[2] = maxitr;
  iparam[3] = 1; // NB; must be 1.
  iparam[6] = mode;

  int processID = theArpackSOE->processID;

#if !defined(ARPACK_ICB)
char which[3];
  if (findSmallest == true)
    strcpy(which, "LM");
  else
    strcpy(which, "SM");

  char bmat  = 'G';
  char howmy = 'A';

#else
  arpack::which  which = findSmallest
                      ? arpack::which::largest_magnitude
                      : arpack::which::smallest_magnitude; 

  arpack::bmat   bmat  = arpack::bmat::generalized; // 'G'
  arpack::howmny howmy = arpack::howmny::ritz_vectors; // 'A'
#endif

  ArpackWorkspace w(n, nev, ArpackWorkspace::Symmetric);
  int ncv = w.ncv;

  if (false) {
    std::vector<double> ones(n, 1.0), mdiag(n);
    theArpackSOE->opM(n, ones.data(), mdiag.data());  // mdiag[i] == M_ii

    double mmax = 0.0;
    for (double mi : mdiag) mmax = std::max(mmax, std::abs(mi));
    const double tau = (mmax > 0 ? 1e-12 * mmax : 0.0); // tune as needed

    std::vector<int> keep; keep.reserve(n);
    for (int i = 0; i < n; ++i)
      if (std::abs(mdiag[i]) > tau) keep.push_back(i);  // translational DOFs

#if 1
    info = 1; // user supplies resid
    std::mt19937_64 rng(1234567);
    std::normal_distribution<double> N(0.0, 1.0);

    std::vector<double> x(n), r(n);
    for (int i = 0; i < n; ++i)
      x[i] = N(rng);

    // Project into Range(M): r = M * x
    theArpackSOE->opM(n, x.data(), r.data());

    // Zero out entries known to be zero-mass (helps numerically)
    for (int i = 0; i < n; ++i)
      if (std::abs(mdiag[i]) <= tau) 
        r[i] = 0.0;

    std::memcpy(w.resid, r.data(), n*sizeof(double));
#endif
  }

  // 
  // I Arnoldi Iteration
  //
  while (true) {
    arpack::saupd(ido, bmat, n, which, nev, tol, w.resid, ncv, 
                  w.v, w.ldv,
                  iparam, ipntr, w.workd, w.workl, w.lworkl, info);

    if (theArpackSOE->checkSameInt(ido) != 1) {
      opserr << "ArpackSolver::solve - ido values not the same\n";
      return -1;
    }

    if (ido == -1) {
      // IDO = -1: compute  Y = OP * X  where
      //           IPNTR(1) is the pointer into WORKD for X,
      //           IPNTR(2) is the pointer into WORKD for Y.
      //           This is for the initialization phase to force the
      //           starting vector into the range of OP.
      
      // y <- x
      theArpackSOE->opM(n, &w.workd[ipntr[0]-1], &w.workd[ipntr[1]-1]); 
    
      theVector.setData(&w.workd[ipntr[1] - 1], size);
     
      if (processID > 0)
        theSOE->zeroB();
      else
        theSOE->setB(theVector);

      theSOE->solve();
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

      myCopy(n, &w.workd[ipntr[2]-1], &w.workd[ipntr[1]-1]);
      theVector.setData(&w.workd[ipntr[1] - 1], size);
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
      theArpackSOE->opM(n, &w.workd[ipntr[0]-1], &w.workd[ipntr[1]-1]);
      continue;
    }
    break;
  }

  //
  //
  //
  // iparam[4]: number of converged eigenvalues
  if (info < 0 || iparam[4] == 0) {
    opserr << "Arpack returned with flag " << info << "\n  ";
    opserr << LookupArpackError(DsaupdErrors, sizeof(DsaupdErrors) / sizeof(ErrorEntry), info);
    opserr << OpenSees::SignalMessageEnd;
    if (info == -9999) {
      return this->solveI(numModes, generalized, findSmallest);
    }
  
    return info;
  }
  //
  // II Compute Eigenvectors
  //
  else {
    if (info == 1) {
      opserr << "ArpackSolver::Maximum number of iteration reached.\n";
    }
    else if (info == 3) {
      opserr << "ArpackSolver::No Shifts could be applied during implicit,";
      opserr << "Arnoldi update, try increasing NCV." << "\n";
    }
    
    double sigma = shift;
    if (iparam[4] > 0) {
      this->numMode = numModes;
      bool rvec = true;
      solution.reserve(n, nev);
      assert(w.ldv >= 1);
      if (rvec)
        assert(w.ldv >= n);

      arpack::seupd(rvec, howmy, 
                    w.select,
                    solution.eigenvalues,  // D
                    solution.eigenvectors, // Z
                    w.ldv,
                    sigma, bmat, n, which,
                    nev, tol, w.resid, ncv, w.v, w.ldv, 
                    iparam, ipntr, 
                    w.workd,
                    w.workl, w.lworkl, info);

      if (info != 0) {
        opserr << "Failed to recover eigenvectors; " << info;
        opserr << LookupArpackError(DseupdErrors, sizeof(DseupdErrors) / sizeof(ErrorEntry), info) 
               << OpenSees::SignalMessageEnd;
        return info;
      }
    }
  }

  return 0;
}


#if 1
int
ArpackSolver::solveI(int numModes, bool generalized, bool findSmallest)
{
  // Asymmetric solver

  // Solve
  theSOE = theArpackSOE->theSOE;
  numMode = 0;
  int n = theArpackSOE->getNumEqn(); // size;

  if (n < numModes || numModes < 1) {
    opserr << "ArpackSolver::solve - no. of modes requested is invalid\n";
    return -1;
  }
  assert (theSOE != nullptr);

  solution.reserve(n, numModes);

  int ido = 0;
  double tol = std::numeric_limits<double>::epsilon();
  int  info = 0;
  int  maxitr = 1000;
  int  mode = 1;
  double sigma = shift;

  int nev = numModes;
  ArpackWorkspace work(n, nev, ArpackWorkspace::NonSymmetric);
  int ncv = work.ncv;

  arpack::bmat  bmat   = arpack::bmat::identity;
  arpack::which which  = arpack::which::largest_magnitude;
  iparam[0] = 1;        // exact shifts
  iparam[2] = maxitr;
  iparam[6] = mode;     // mode 1 (std. problem)


  //
  // I Arnoldi Iteration
  //
  while (true) { 
    arpack::naupd(ido, bmat, n, which, nev, tol,
                  work.resid, ncv, 
                  work.v, work.ldv, 
                  iparam, ipntr,
                  work.workd, work.workl, work.lworkl, 
                  info);

    assert(ipntr[0] >= 1 && ipntr[0] <= 3*n);
    assert(ipntr[1] >= 1 && ipntr[1] <= 3*n);

    if (ido == -1 || ido == 1) {
      // x is at workd[ipntr[0]-1]
      // y should be written to workd[ipntr[1]-1]
      double* x = &work.workd[ipntr[0]-1];
      double* y = &work.workd[ipntr[1]-1];

      // y = M * x
      if (generalized == true)
        theArpackSOE->opM(n, x, y);
      else
        // theArpackSOE->opM(n, x, y);
        this->myCopy(n, x, y);

      // solve (K - σ M) z = y, store solution (z) into y
      if (true) { // sigma != 0.0 && generalized == true) {
        theVector.setData(y, n);
        theSOE->setB(theVector);
        theSOE->solve();             // factor Aσ and solve
        const Vector& X = theSOE->getX();
        for (int i=0; i<n; i++)
          y[i] = X[i];
      }

      continue;
    }

    break;  // ido == 99
  }

  //
  // II Compute Eigenvectors
  //
  if (info != 0) {
    opserr << LookupArpackError(DnaupdErrors, sizeof(DnaupdErrors) / sizeof(ErrorEntry), info)
           << OpenSees::SignalMessageEnd;
    return info;
  }
  else {
    bool rvec = true;
    arpack::howmny howmny = arpack::howmny::ritz_vectors;
    double* di = new double[nev]{};
    double* workev = new double[3*ncv]{};
    double* z = new double[n*(nev+1)]{};

    arpack::neupd(rvec, 
                  howmny, 
                  work.select, 
                  solution.eigenvalues, di, z, 
                  work.ldv,
                  sigma, 0.0, workev,
                  bmat, n, which, nev, tol, 
                  work.resid, work.ncv, work.v, work.ldv,
                  iparam, ipntr, 
                  work.workd, work.workl, work.lworkl, info);
    
    if (info == 0) {
      numMode = iparam[4];
      for (int i=0; i<nev; i++) {
        if (std::abs(solution.eigenvalues[i]) < 1e-16)
          solution.eigenvalues[i] = 0.0;
        else
          solution.eigenvalues[i] = sigma + 1.0/solution.eigenvalues[i];
      }

      // populate eigenvectors (all real)
      {
        const int ldz = work.ldv;  // leading dimension of z from neupd
        for (int j = 0; j < nev; ++j) {
          const double* zj = &z[j * ldz];     // column j of z
          double* vj       = &solution.eigenvectors[j * n]; // column j of output
          std::memcpy(vj, zj, sizeof(double) * n);
        }
      }
    }
    delete [] z;
    delete [] di;
    delete [] workev;
  }

  if (info != 0) {
    opserr << LookupArpackError(DneupdErrors, sizeof(DneupdErrors) / sizeof(ErrorEntry), info)
           << OpenSees::SignalMessageEnd;

    solution.zero();
  }

  return info;
}
#endif




void
ArpackSolver::myCopy(int n, double *v, double *result)
{
  for (int i=0; i<n; i++)
    result[i] = v[i];
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
  this->getEigenvector(mode, theVector);

  return theVector;
}

int 
ArpackSolver::getEigenvector(int mode, Vector &eigenvector)
{
  if (mode <= 0 || mode > numMode) {
    return -1;
  }
  
  int index = (mode - 1) * size;
  
  eigenvector.setData(&solution.eigenvectors[index], size);

  return 0;  
}


double
ArpackSolver::getEigenvalue(int mode)
{
  if (mode <= 0 || mode > numMode) {
    opserr << "ArpackSOE::getEigenvalue() - mode is out of range(1 - nev)";
    return -1;
  }
  
  if (mode <= numMode)
    return solution.eigenvalues[mode-1];
  else {
    opserr << "ArpackSOE::getEigenvalue() - eigenvalues not yet determined";
    return -2;
  }      
}


int
ArpackSolver::setSize()
{
  size = theArpackSOE->getNumEqn(); // Msize;
  return 0;
}


int    
ArpackSolver::sendSelf(int commitTag, Channel &)
{
  return 0;
}

int
ArpackSolver::recvSelf(int commitTag, Channel &, FEM_ObjectBroker &)
{
  return 0;
}


