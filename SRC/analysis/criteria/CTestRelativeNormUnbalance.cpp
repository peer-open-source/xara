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
#include <CTestRelativeNormUnbalance.h>
#include <Vector.h>
#include <Channel.h>
#include <EquiSolnAlgo.h>
#include <LinearSOE.h>
#include <Logging.h>



CTestRelativeNormUnbalance::CTestRelativeNormUnbalance(double theTol, int maxIter, int printIt, int normType)
    : ConvergenceTest(CONVERGENCE_TEST_CTestRelativeNormUnbalance),
    tol(theTol), maxNumIter(maxIter), currentIter(0), printFlag(printIt),
    norms(maxNumIter+1), norm0(0.0), nType(normType),
    first_step(true)
{

}


CTestRelativeNormUnbalance::~CTestRelativeNormUnbalance()
{

}


ConvergenceTest* 
CTestRelativeNormUnbalance::getCopy(int iterations)
{
  return new CTestRelativeNormUnbalance(this->tol, iterations, 9, this->nType); //this->printFlag
}


void
CTestRelativeNormUnbalance::setTolerance(double newTol)
{
  tol = newTol;
}


int
CTestRelativeNormUnbalance::start(LinearSOE& theSOE)
{
  double norm_last = 0.0;
  if (currentIter != 0)
      norm_last = norms(currentIter-1)*norm0;

  norms.Zero();
  currentIter = 1;

  // determine the initial norm .. the the norm of the initial unbalance
  const Vector &b = theSOE.getB();
  double norm = b.pNorm(nType);

  // if (currentIter <= maxNumIter)
  //     norms(0) = norm;
  // if (first_step)
  norm0 += norm;// - norm_last;

  if (printFlag & ConvergenceTest::PrintTest) {
    pstream << LOG_ITERATE << "Iter: " << pad(0)
            << ", R : " << pad(b.pNorm(nType)) 
            << ", R0: " << pad(norm0) 
            << "\n";
  }
  // first_step = false;
  return 0;
}


int
CTestRelativeNormUnbalance::test(const Vector& b, const Vector& x)
{
  // check to ensure the algo does invoke start() - this is needed otherwise
  // may never get convergence later on in analysis!
  if (currentIter == 0) {
      opserr << "WARNING: CTestRelativeNormUnbalance::test - start() was never invoked.\n";
      return -2;
  }

  // get the B vector & determine it's norm & save the value in norms vector
  // const Vector &x = theSOE.getB();
  double norm = b.pNorm(nType);

  // determine the ratio
  if (norm0 != 0.0)
    norm /= norm0;

  if (currentIter <= maxNumIter)
    norms(currentIter-1) = norm;

  // print the data if required
  if (printFlag & ConvergenceTest::PrintTest) {
    pstream << LOG_ITERATE 
            << "Iter: "    << pad(currentIter)
            << ", |dR|/|dR0|: " << pad(norm) 
            << "\n";
  }
  if (printFlag & ConvergenceTest::PrintTest02) {
    pstream << LOG_ITERATE 
            << "Iter: "     << pad(currentIter)
            << ", |dR|/|dR0|: "  << pad(norm) 
            << "\n"
            << "\tNorm deltaX: " << pad(x.pNorm(nType)) 
            << ", Norm deltaR: " << pad(norm) 
            << "\n"
            << "\tdeltaX: "      << x 
            << "\tdeltaR: "      << b;
  }

  //
  // check if the algorithm converged
  //

  // if converged - print & return ok

  if (norm <= tol) { // the algorithm converged
    // do some printing first
    if (printFlag & ConvergenceTest::PrintTest || printFlag & ConvergenceTest::PrintTest02)
        pstream << "\n";
    if (printFlag & ConvergenceTest::PrintSuccess) {
        pstream << LOG_SUCCESS 
                << "Iter: "    << pad(currentIter)
                << ", |dR|/|dR0|: " << pad(norm) 
                << "\n";
    }

    // return the number of times test has been called
    return currentIter;
  }

  // algo failed to converged after specified number of iterations - but RETURN OK
  else if ((printFlag & ConvergenceTest::AlwaysSucceed) && currentIter >= maxNumIter) {
    if (printFlag & ConvergenceTest::PrintFailure) {
      pstream << LOG_FAILURE
              << ", dR/dR0: "       << pad(norm)
              << ", Norm dX: "  << pad(x.pNorm(nType)) 
              << "\n";
    }
    return currentIter;
  }

  // algo failed to converged after specified number of iterations - return FAILURE -2
  else if (currentIter >= maxNumIter) { // the algorithm failed to converge
    if (printFlag & ConvergenceTest::PrintFailure) {
        pstream << LOG_FAILURE
                //<< "criteria CTestRelativeNormUnbalance"
                // << LOG_CONTINUE
                << "Iter: "         << pad(currentIter)
                << ", |dR|/|dR0|: " << pad(norm) 
                << "\n";
    }
    // we increment in case analysis does not check for convergence
    currentIter++;
    return ConvergenceTest::Failure;
  }

  // algorithm not yet converged - increment counter and return -1
  else {
    currentIter++;
    return ConvergenceTest::Continue;
  }
}



int
CTestRelativeNormUnbalance::getNumTests()
{
  return currentIter;
}


int
CTestRelativeNormUnbalance::getMaxNumTests()
{
  return maxNumIter;
}


double
CTestRelativeNormUnbalance::getRatioNumToMax()
{
  return currentIter/double(maxNumIter);
}


const Vector&
CTestRelativeNormUnbalance::getNorms()
{
  return norms;
}


