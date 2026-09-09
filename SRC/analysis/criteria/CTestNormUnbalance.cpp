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
//
#include <CTestNormUnbalance.h>
#include <Vector.h>
#include <Channel.h>
#include <EquiSolnAlgo.h>
#include <LinearSOE.h>
#include <Logging.h>
#include <iostream>
#include <fstream>


CTestNormUnbalance::CTestNormUnbalance(double theTol, int maxIter, int printIt, int normType, int maxincr, double max)
    : ConvergenceTest(CONVERGENCE_TEST_CTestNormUnbalance),
      tol(theTol), maxTol(max), maxNumIter(maxIter), currentIter(0), printFlag(printIt),
      nType(normType), norms(maxNumIter), maxIncr(maxincr), numIncr(0)
{
  if (maxIncr < 0) {
    maxIncr = maxNumIter;
  }
}


CTestNormUnbalance::~CTestNormUnbalance()
{

}


ConvergenceTest* 
CTestNormUnbalance::getCopy(int iterations)
{
  return new CTestNormUnbalance(this->tol, iterations, 
                                9,// this->printFlag, 
                                this->nType, this->maxIncr, this->maxTol);

}


void 
CTestNormUnbalance::setTolerance(double newTol)
{
  tol = newTol;
}


int
CTestNormUnbalance::start(LinearSOE& theSOE)
{
  // set iteration count = 1
  norms.Zero();
  currentIter = 1;
  numIncr = 0;

  if (printFlag & ConvergenceTest::PrintTest) {
    pstream << LOG_ITERATE << "Iter: " << pad(0);
    pstream << ", R : " << pad(theSOE.getB().pNorm(nType)) << "\n";
  }
  return 0;
}

int
CTestNormUnbalance::test(const Vector& b, const Vector& x)
{

  // check to ensure the algo does invoke start() - this is needed otherwise
  // may never get convergence later on in analysis!
  if (currentIter == 0) {
    opserr << "WARNING: CTestNormUnbalance::test() - start() was never invoked.\n";
    return -2;
  }

  // get the B vector & determine it's norm & save the value in norms vector
  // const Vector &b = theSOE.getB();
  double norm = b.pNorm(nType);
  if (currentIter <= maxNumIter)
    norms(currentIter-1) = norm;

  if (currentIter > 1) {
    // check if the norm is increasing
    if (norms(currentIter-2) < norm) {
      numIncr++;
    }
  }

  // print the data if required
  if (printFlag & ConvergenceTest::PrintTest) {
    pstream << LOG_ITERATE << "Iter: " << pad(currentIter);
    pstream << ", R : " << pad(norm);
    pstream << ", dX: " << x.pNorm(nType) << "\n";
  }
  if (printFlag & ConvergenceTest::PrintTest02) {
    pstream << LOG_ITERATE << "Iter: " << pad(currentIter);
    pstream << ", Norm: " << pad(norm) << "\n";
    pstream << "\t|dX|: " << x.pNorm(nType) 
            << ", |R|: " << pad(norm) << "\n";
    pstream << "\tdX: " << x << "\tR: " << b;
  }

  if (printFlag == 7) {
    // std::ofstream outDu;
    // std::ofstream outDp;

    // if (currentIter == 1) {
    //   outDu.open("dX.out",std::ios::out);
    //   outDp.open("dP.out", std::ios::out);
    // } else {
    //   outDu.open("dX.out",std::ios::app);
    //   outDp.open("dP.out", std::ios::app);
    // }
    // const Vector &Du = theSOE.getX();
    // const Vector &Dp = theSOE.getB();
    // for (int i=0; i<Du.Size(); i++) {
    //   outDu << Du[i] << " ";
    //   outDp << Dp[i] << " ";
    // }
    // outDu << "\n";
    // outDp << "\n";
    // outDu.close();
    // outDp.close();
  }

  //
  // check if the algorithm converged
  //

  // if converged - print & return ok
  if (norm <= tol) {

    // Print data
    if (printFlag & ConvergenceTest::PrintTest || printFlag & ConvergenceTest::PrintTest02)
      pstream << "\n";
    if (printFlag & ConvergenceTest::PrintSuccess || printFlag == 7) {
      pstream << LOG_SUCCESS << "Iter: " << pad(currentIter);
      pstream << ", Norm: " << pad(norm) << " (max: " << tol;
      pstream << ", Norm deltaX: " << x.pNorm(nType) << ")\n";
    }

    // return the number of times test has been called
    return currentIter;
  }

  // algo failed to converged after specified number of iterations - but RETURN OK
  else if ((printFlag & ConvergenceTest::AlwaysSucceed) && (currentIter >= maxNumIter||numIncr>=maxIncr)) {
    if (printFlag & ConvergenceTest::PrintFailure) {
      pstream << LOG_FAILURE
              << ", Norm: " << pad(norm) 
              << ", Norm deltaX: " << pad(x.pNorm(nType))
              << "\n";
    }
    return currentIter;
  }

  // algo failed to converged after specified number of iterations - return FAILURE -2
  else if (currentIter >= maxNumIter || numIncr >= maxIncr || norm > maxTol) { // the algorithm failed to converge
      if (printFlag & ConvergenceTest::PrintFailure) {
          pstream << LOG_FAILURE
                  //<< "criteria CTestNormUnbalance"
                  // << LOG_CONTINUE
                  << "Iter: "           << pad(currentIter)
                  << ", Norm: "         << pad(norm)
                  << ", Norm deltaX: "  << pad(x.pNorm(nType)) 
                  << "\n";
      }
      currentIter++;  // we increment in case analysis does not check for convergence
      return ConvergenceTest::Failure;
  }

  // algorithm not yet converged - increment counter and return -1
  else {
    currentIter++;
    return ConvergenceTest::Continue;
  }
}



int
CTestNormUnbalance::getNumTests()
{
  return currentIter;
}


int
CTestNormUnbalance::getMaxNumTests()
{
  return maxNumIter;
}


double
CTestNormUnbalance::getRatioNumToMax()
{
    double div = maxNumIter;
    return currentIter/div;
}


const Vector&
CTestNormUnbalance::getNorms()
{
  return norms;
}
