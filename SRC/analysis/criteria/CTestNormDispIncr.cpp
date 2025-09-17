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
#include <CTestNormDispIncr.h>
#include <Vector.h>
#include <Channel.h>
#include <EquiSolnAlgo.h>
#include <LinearSOE.h>
#include <Logging.h>

CTestNormDispIncr::CTestNormDispIncr()
    : ConvergenceTest(CONVERGENCE_TEST_CTestNormDispIncr),
      tol(0), maxTol(OPS_MAXTOL), maxNumIter(0), currentIter(0), printFlag(0),
      norms(25), nType(2)
{

}


CTestNormDispIncr::CTestNormDispIncr(double theTol, int maxIter, int printIt, int normType, double max)
    : ConvergenceTest(CONVERGENCE_TEST_CTestNormDispIncr),
      tol(theTol), maxTol(max), maxNumIter(maxIter), currentIter(0), printFlag(printIt),
      nType(normType), norms(maxIter)
{

}


CTestNormDispIncr::~CTestNormDispIncr()
{

}


ConvergenceTest*
CTestNormDispIncr::getCopy(int iterations)
{
    return new CTestNormDispIncr(this->tol, iterations, 0, this->nType, this->maxTol) ;
}


void
CTestNormDispIncr::setTolerance(double newTol)
{
    tol = newTol;
}



int
CTestNormDispIncr::test(LinearSOE& theSOE)
{

    // check to ensure the algo does invoke start() - this is needed otherwise
    // may never get convergence later on in analysis!
    if (currentIter == 0) {
        opserr << "WARNING: CTestNormDispIncr::test - start() was never invoked.\n";
        return -2;
    }

    // get the X vector & determine it's norm & save the value in norms vector
    const Vector &x = theSOE.getX();
    double norm = x.pNorm(nType);
    if (currentIter <= maxNumIter)
        norms(currentIter-1) = norm;

    // print the data if required
    if (printFlag & ConvergenceTest::PrintTest) {
        opserr << LOG_ITERATE 
               << "Iter: "           << pad(currentIter)
               << ", Norm: "         << pad(norm) 
               << ", Norm deltaR: "  << pad(theSOE.getB().pNorm(nType))
               << endln;
    }
    else if (printFlag & ConvergenceTest::PrintTest02) {
        opserr << LOG_ITERATE 
               << "Iter: "     << currentIter
               << ", Norm: "         << pad(norm) 
               << endln;
        opserr << "\tNorm deltaX: "  << pad(norm) 
               << ", Norm deltaR: "  << pad(theSOE.getB().pNorm(nType))
               << endln
               << "\tdeltaX: "       << x
               << "\tdeltaR: "       << theSOE.getB();
    }

    //
    // check if the algorithm converged
    //

    // if converged - print & return ok
    if (norm <= tol) {

        // do some printing first
        if (printFlag & ConvergenceTest::PrintTest || printFlag & ConvergenceTest::PrintTest02)
            opserr << endln;

        if (printFlag & ConvergenceTest::PrintSuccess) {
            opserr << LOG_SUCCESS 
                   << "Iter: "          << pad(currentIter)
                   << ", Norm: "        << pad(norm)
                   << ", Norm deltaR: " << pad(theSOE.getB().pNorm(nType))
                   << endln;
        }

        // return the number of times test has been called
        return currentIter;
    }

    // failed to converged after specified number of iterations - but RETURN OK
    else if ((printFlag & ConvergenceTest::AlwaysSucceed) && currentIter >= maxNumIter) {
        if (printFlag & ConvergenceTest::PrintFailure) {
            opserr << LOG_FAILURE
                   << ", Norm: " << pad(norm)  // << " (max: " << tol;
                   << ", Norm deltaR: " << pad(theSOE.getB().pNorm(nType))
                   << LOG_CONTINUE
                   << "failed to converge but going on - "
                   << endln;
        }
        return currentIter;
    }

    // failed to converged after specified number of iterations - return FAILURE -2
    else if (currentIter >= maxNumIter || norm > maxTol) { // failes to converge
        if (printFlag & ConvergenceTest::PrintFailure) {
          opserr << LOG_FAILURE 
                 //<< "criteria CTestNormDispIncr" 
                 // << LOG_CONTINUE
                 << "Iter: "             << pad(currentIter)
                 << ", Norm: "           << pad(norm)
                 << ", Norm deltaR: "    << pad(theSOE.getB().pNorm(nType))
                 << endln;
        }
        currentIter++;
        return ConvergenceTest::Failure;
    }

    // algorithm not yet converged - increment counter and return -1
    else {
        currentIter++;
        return ConvergenceTest::Continue;
    }
}


int CTestNormDispIncr::start(LinearSOE& theSOE)
{
    // set iteration count = 1
    norms.Zero();
    currentIter = 1;
    return 0;
}


int CTestNormDispIncr::getNumTests()
{
    return currentIter;
}


int CTestNormDispIncr::getMaxNumTests(void)
{
    return maxNumIter;
}


double CTestNormDispIncr::getRatioNumToMax(void)
{
    double div = maxNumIter;
    return currentIter/div;
}


const Vector& CTestNormDispIncr::getNorms()
{
    return norms;
}


int CTestNormDispIncr::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  Vector x(5);
  x(0) = tol;
  x(1) = maxNumIter;
  x(2) = printFlag;
  x(3) = nType;
  x(4) = maxTol;
  res = theChannel.sendVector(this->getDbTag(), cTag, x);
  if (res < 0)
    opserr << "CTestNormDispIncr::sendSelf() - failed to send data\n";

  return res;
}

int
CTestNormDispIncr::recvSelf(int cTag, Channel &theChannel,
                          FEM_ObjectBroker &theBroker)
{
  int res = 0;
  Vector x(5);
  res = theChannel.recvVector(this->getDbTag(), cTag, x);


  if (res < 0) {
    opserr << "CTestNormDispIncr::sendSelf() - failed to send data\n";
    tol = 1.0e-8;
    maxNumIter = 25;
    printFlag = 0;
    nType = 2;
    norms.resize(maxNumIter);
  } else {
    tol = x(0);
    maxNumIter = (int)x(1);
    printFlag = (int)x(2);
    nType = (int)x(3);
    norms.resize(maxNumIter);
    maxTol = x(4);
  }
  return res;
}


