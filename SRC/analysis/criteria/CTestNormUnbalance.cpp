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

CTestNormUnbalance::CTestNormUnbalance()
    : ConvergenceTest(CONVERGENCE_TEST_CTestNormUnbalance),
      tol(0.0), maxTol(OPS_MAXTOL), maxNumIter(0), currentIter(0), printFlag(0),
      norms(1), nType(2), maxIncr(0), numIncr(0)
{

}


CTestNormUnbalance::CTestNormUnbalance(double theTol, int maxIter, int printIt, int normType, int maxincr, double max)
    : ConvergenceTest(CONVERGENCE_TEST_CTestNormUnbalance),
      tol(theTol), maxTol(max), maxNumIter(maxIter), currentIter(0), printFlag(printIt),
      nType(normType), norms(maxNumIter), maxIncr(maxincr), numIncr(0)
{
    if(maxIncr < 0) {
        maxIncr = maxNumIter;
    }
}


CTestNormUnbalance::~CTestNormUnbalance()
{

}


ConvergenceTest* CTestNormUnbalance::getCopy(int iterations)
{
    CTestNormUnbalance *theCopy ;
    theCopy = new CTestNormUnbalance(this->tol, iterations, this->printFlag, this->nType, this->maxIncr, this->maxTol);

    return theCopy ;
}


void CTestNormUnbalance::setTolerance(double newTol)
{
    tol = newTol;
}


int
CTestNormUnbalance::test(LinearSOE& theSOE)
{

    // check to ensure the algo does invoke start() - this is needed otherwise
    // may never get convergence later on in analysis!
    if (currentIter == 0) {
        opserr << "WARNING: CTestNormUnbalance::test() - start() was never invoked.\n";
        return -2;
    }

    // get the B vector & determine it's norm & save the value in norms vector
    const Vector &x = theSOE.getB();
    double norm = x.pNorm(nType);
    if (currentIter <= maxNumIter)
        norms(currentIter-1) = norm;

    if (currentIter > 1) {
        if (norms(currentIter-2) < norm) {
            numIncr++;
        }
    }

    // print the data if required
    if (printFlag & ConvergenceTest::PrintTest) {
        opserr << LOG_ITERATE << "Iter: " << pad(currentIter);
        opserr << ", Norm: " << pad(norm) << " (max: " << tol;
        opserr << ", Norm deltaX: " << theSOE.getX().pNorm(nType) << ")\n";
    }
    if (printFlag & ConvergenceTest::PrintTest02) {
        opserr << LOG_ITERATE << "Iter: " << pad(currentIter);
        opserr << ", Norm: " << pad(norm) << " (max: " << tol << ")\n";
        opserr << "\tNorm deltaX: " << theSOE.getX().pNorm(nType) << ", Norm deltaR: " << pad(norm) << "\n";
        opserr << "\tdeltaX: " << theSOE.getX() << "\tdeltaR: " << x;
    }

    if (printFlag == 7) {
      std::ofstream outDu;
      std::ofstream outDp;

      if (currentIter == 1) {
        outDu.open("dX.out",std::ios::out);
        outDp.open("dP.out", std::ios::out);
      } else {
        outDu.open("dX.out",std::ios::app);
        outDp.open("dP.out", std::ios::app);
      }
      const Vector &Du = theSOE.getX();
      const Vector &Dp = theSOE.getB();
      for (int i=0; i<Du.Size(); i++) {
        outDu << Du[i] << " ";
        outDp << Dp[i] << " ";
      }
      outDu << "\n";
      outDp << "\n";
      outDu.close();
      outDp.close();
    }

    //
    // check if the algorithm converged
    //

    // if converged - print & return ok
    if (norm <= tol) {

        // do some printing first
        if (printFlag & ConvergenceTest::PrintTest || printFlag & ConvergenceTest::PrintTest02)
            opserr << "\n";
        if (printFlag & ConvergenceTest::PrintSuccess || printFlag == 7) {
            opserr << LOG_SUCCESS << "Iter: " << pad(currentIter);
            opserr << ", Norm: " << pad(norm) << " (max: " << tol;
            opserr << ", Norm deltaX: " << theSOE.getX().pNorm(nType) << ")\n";
        }

        // return the number of times test has been called
        return currentIter;
    }

    // algo failed to converged after specified number of iterations - but RETURN OK
    else if ((printFlag & ConvergenceTest::AlwaysSucceed) && (currentIter >= maxNumIter||numIncr>=maxIncr)) {
        if (printFlag & ConvergenceTest::PrintFailure) {
            opserr << LOG_FAILURE
                   //<< "criteria CTestNormUnbalance but going on -";
                   << ", Norm: " << pad(norm) 
                   << ", Norm deltaX: " << pad(theSOE.getX().pNorm(nType))
                   << "\n";
        }
        return currentIter;
    }

    // algo failed to converged after specified number of iterations - return FAILURE -2
    else if (currentIter >= maxNumIter || numIncr >= maxIncr || norm > maxTol) { // the algorithm failed to converge
        if (printFlag & ConvergenceTest::PrintFailure) {
            opserr << LOG_FAILURE 
                   //<< "criteria CTestNormUnbalance"
                   // << LOG_CONTINUE
                   << "Iter: "           << pad(currentIter)
                   << ", Norm: "         << pad(norm)
                   << ", Norm deltaX: "  << pad(theSOE.getX().pNorm(nType)) 
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


int CTestNormUnbalance::start(LinearSOE& theSOE)
{
    // set iteration count = 1
    norms.Zero();
    currentIter = 1;
    numIncr = 0;
    return 0;
}


int CTestNormUnbalance::getNumTests()
{
    return currentIter;
}


int CTestNormUnbalance::getMaxNumTests(void)
{
    return maxNumIter;
}


double CTestNormUnbalance::getRatioNumToMax(void)
{
    double div = maxNumIter;
    return currentIter/div;
}


const Vector& CTestNormUnbalance::getNorms()
{
    return norms;
}


int CTestNormUnbalance::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    static Vector x(5);
    x(0) = tol;
    x(1) = maxNumIter;
    x(2) = printFlag;
    x(3) = nType;
    x(4) = maxTol;
    res = theChannel.sendVector(this->getDbTag(), cTag, x);
    if (res < 0)
        opserr << "CTestNormUnbalance::sendSelf() - failed to send data\n";

    return res;
}


int CTestNormUnbalance::recvSelf(int cTag, Channel &theChannel,
    FEM_ObjectBroker &theBroker)
{
    int res = 0;
    static Vector x(5);
    res = theChannel.recvVector(this->getDbTag(), cTag, x);

    if (res < 0) {
        opserr << "CTestNormUnbalance::sendSelf() - failed to send data\n";
        tol = 1.0e-8;
        maxNumIter = 25;
        printFlag = 0;
        nType = 2;
    }
    else {
        tol = x(0);
        maxNumIter = (int) x(1);
        printFlag = (int) x(2);
        nType = (int) x(3);
        norms.resize(maxNumIter);
        maxTol = x(4);
    }
    return res;
}
