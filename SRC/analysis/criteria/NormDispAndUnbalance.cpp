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
#include <NormDispAndUnbalance.h>
#include <Vector.h>
#include <EquiSolnAlgo.h>
#include <LinearSOE.h>
#include <Logging.h>

NormDispAndUnbalance::NormDispAndUnbalance(double theTolDisp, double theTolUnbalance, int maxIter, int printIt, int normType, int maxincr)
  : ConvergenceTest(CONVERGENCE_TEST_NormDispAndUnbalance),
    tolDisp(theTolDisp), tolUnbalance(theTolUnbalance),
    maxNumIter(maxIter), currentIter(0), printFlag(printIt),
    norms(2*maxIter), nType(normType), maxIncr(maxincr), numIncr(0)
{
    if(maxIncr < 0) maxIncr = maxNumIter;
}


NormDispAndUnbalance::~NormDispAndUnbalance()
{

}


ConvergenceTest*
NormDispAndUnbalance::getCopy(int iterations)
{
  return new NormDispAndUnbalance(this->tolDisp,
                                   this->tolUnbalance,
                                   iterations,
                                   this->printFlag,
                                   this->nType,
                                   this->maxIncr);
}


void 
NormDispAndUnbalance::setTolerance(double newTolDisp)
{
  tolDisp = newTolDisp;
}



int
NormDispAndUnbalance::start(LinearSOE& theSOE)
{
  // set iteration count = 1
  norms.Zero();
  currentIter = 1;
  numIncr = 0;
  return 0;
}


int
NormDispAndUnbalance::test(const Vector& b, const Vector& x)
{
  // check to ensure the algo does invoke start() - this is needed otherwise
  // may never get convergence later on in analysis!
  if (currentIter == 0) {
      opserr << "WARNING: NormDispAndUnbalance::test() - start() was never invoked.\n";
      return -2;
  }

  // get the X vector & determine it's norm & save the value in norms vector
  // const Vector &x = theSOE.getX();
  double normX = x.pNorm(nType);

  // const Vector &b = theSOE.getB();
  double normB = b.pNorm(nType);

  if ((currentIter>1 && norms(currentIter-2)<normX) || 
      (currentIter>1 && norms(maxNumIter+currentIter-2)<normB)) {
      numIncr++;
  }

  if (currentIter <= maxNumIter) {
      norms(currentIter-1) = normX;
      norms(maxNumIter+currentIter-1) = normB;
  }

  // print the data if required
  if (printFlag == ConvergenceTest::PrintTest) {
      pstream << LOG_ITERATE;
      pstream << "Iter: "    << pad(currentIter);
      pstream << ", NormX: " << pad(normX);
      pstream << ", NormB: " << pad(normB)  
              << ", NormIncr: " << numIncr << "\n";
  }
  if (printFlag == ConvergenceTest::PrintTest02) {
      pstream << LOG_ITERATE;
      pstream << "Iter: " << pad(currentIter);
      pstream << ", NormX: " << pad(normX);
      pstream << ", NormB: " << pad(normB)  
              << ", NormIncr: " << numIncr << "\n";
      pstream << "\tdeltaX: " << x << "\tdeltaR: " << b;
  }

  //
  // check if the algorithm converged
  //

  // if converged - print & return ok
  if (normX <= tolDisp && normB <= tolUnbalance) {

      // do some printing first
      if (printFlag == ConvergenceTest::PrintTest || printFlag == ConvergenceTest::PrintTest02)
          pstream << "\n";
      if (printFlag == ConvergenceTest::PrintSuccess) {
          pstream << "NormDispAndUnbalance::test() - iteration: " << pad(currentIter);
          pstream << ", NormX: " << pad(normX);
          pstream << ", NormB: " << pad(normB)  
                  << ", NormIncr: " << numIncr << "\n";
      }

      // return the number of times test has been called
      return currentIter;
  }

  // algo failed to converged after specified number of iterations - but RETURN OK
  else if ((printFlag == ConvergenceTest::AlwaysSucceed) && (currentIter >= maxNumIter || numIncr > maxIncr)) {
      if (printFlag & ConvergenceTest::PrintFailure) {
          pstream << "WARNING Failed to converge with criteria NormDispAndUnbalance but going on - ";
          pstream << ", NormX: " << pad(normX);
          pstream << ", NormB: " << pad(normB)  
                  << ", NormIncr: " << numIncr << "\n";
      }
      return currentIter;
  }

  // algo failed to converged after specified number of iterations - return FAILURE -2
  else if (currentIter >= maxNumIter || numIncr > maxIncr) { // failes to converge
      pstream << LOG_FAILURE;
      pstream << "Iter: " << pad(currentIter);
      pstream << ", NormX: " << pad(normX);
      pstream << ", NormB: " << pad(normB)  
              << ", NormIncr: " << numIncr << "\n";
      currentIter++;
      return ConvergenceTest::Failure;
  }

  // algorithm not yet converged - increment counter and return -1
  else {
      currentIter++;
      return ConvergenceTest::Continue;
  }
}


int NormDispAndUnbalance::getNumTests()
{
  return currentIter;
}


int NormDispAndUnbalance::getMaxNumTests()
{
    return maxNumIter;
}


double NormDispAndUnbalance::getRatioNumToMax()
{
    double div = maxNumIter;
    return currentIter/div;
}


const Vector& NormDispAndUnbalance::getNorms()
{
    return norms;
}




