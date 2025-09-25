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
// Date: 02/02
//
#include <CTestRelativeEnergyIncr.h>
#include <Vector.h>
#include <Channel.h>
#include <EquiSolnAlgo.h>
#include <LinearSOE.h>
#include <Logging.h>

CTestRelativeEnergyIncr::CTestRelativeEnergyIncr()
    : ConvergenceTest(CONVERGENCE_TEST_CTestRelativeEnergyIncr),
    tol(0), maxNumIter(0), currentIter(0), printFlag(0),
    norms(1), norm0(0.0), nType(2)
{

}


CTestRelativeEnergyIncr::CTestRelativeEnergyIncr(double theTol, int maxIter, int printIt, int normType)
    : ConvergenceTest(CONVERGENCE_TEST_CTestRelativeEnergyIncr),
    tol(theTol), maxNumIter(maxIter), currentIter(0),printFlag(printIt),
    norms(maxNumIter), norm0(0.0), nType(normType)
{

}


CTestRelativeEnergyIncr::~CTestRelativeEnergyIncr()
{

}


ConvergenceTest*
CTestRelativeEnergyIncr::getCopy(int iterations)
{
    return new CTestRelativeEnergyIncr(this->tol, iterations, this->printFlag, this->nType);
}


void CTestRelativeEnergyIncr::setTolerance(double newTol)
{
    tol = newTol;
}


int
CTestRelativeEnergyIncr::test(LinearSOE& theSOE)
{
    // check to ensure the algo does invoke start() - this is needed otherwise
    // may never get convergence later on in analysis!
    if (currentIter == 0) {
        opserr << "WARNING: CTestRelativeEnergyIncr::test() - start() was never invoked.\n";
        return -2;
    }


    // determine the energy & save value in norms vector
    const Vector &b = theSOE.getB();
    const Vector &x = theSOE.getX();
    double product = x ^ b;
    if (product < 0.0)
        product *= -0.5;
    else
        product *= 0.5;

    if (currentIter <= maxNumIter)
        norms(currentIter-1) = product;

    // if first pass through .. set norm0
    if (currentIter == 1) {
        norm0 = product;
    }

    // get ratio
    if (norm0 != 0.0)
        product /= norm0;

    // print the data if required
    if (printFlag & ConvergenceTest::PrintTest) {
        pstream << LOG_ITERATE
               << "Iter: "            << pad(currentIter)
               << ", dX*dR/dX1*dR1: " << pad(product)
               << "\n";
    }
    if (printFlag & ConvergenceTest::PrintTest02) {
        pstream << LOG_ITERATE
               << "Iter: "            << pad(currentIter)
               << ", dX*dR/dX1*dR1: " << pad(product)
               << "\n"
               << ", Norm deltaX: "   << pad(x.pNorm(nType))
               << ", Norm deltaR: "   << pad(b.pNorm(nType)) 
               << "\n"
               << "\tdeltaX: "        << x 
               << "\tdeltaR: "        << b;
    }

    //
    // check if the algorithm converged
    //

    // if converged - print & return ok
    if (product <= tol) {

        // do some printing first
        if (printFlag & ConvergenceTest::PrintTest || printFlag & ConvergenceTest::PrintTest02)
            pstream << "\n";
        if (printFlag & ConvergenceTest::PrintSuccess) {
            pstream << LOG_SUCCESS
                   << "Iter: "           << pad(currentIter)
                   << ", dX*dR/dX1*dR1: " << pad(product)
                   << "\n";
        }

        // return the number of times test has been called - SUCCESSFULL
        return currentIter;
    }

    // Failed to converged after specified number of iterations - but RETURN OK
    else if ((printFlag & ConvergenceTest::AlwaysSucceed) && currentIter >= maxNumIter) {
        if (printFlag & ConvergenceTest::PrintFailure) {
            pstream << LOG_FAILURE
                   //<< "criteria CTestRelativeEnergyIncr but goin on -"
                   << "Iter: "            << pad(currentIter)
                   << ", Norm dX: "  << pad(x.pNorm(nType))
                   << ", Norm dR: "  << pad(b.pNorm(nType))
                   << ", dX*dR/dX1*dR1: " << pad(product)
                   << "\n";
        }
        return currentIter;
    }

    // algo failed to converged after specified number of iterations - return FAILURE -2
    else if (currentIter >= maxNumIter) { // >= in case algorithm does not check
        if (printFlag & ConvergenceTest::PrintFailure) {
            pstream << LOG_FAILURE
                   // << LOG_CONTINUE
                   << "Iter: "           << pad(currentIter)
                   << ", dX*dR/dX1*dR1: " << pad(product)
                   << ", Norm deltaX: "  << pad(x.pNorm(nType))
                   // << LOG_CONTINUE
                   <<   "Norm deltaR: "  << pad(b.pNorm(nType))
                   << "\n";
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


int
CTestRelativeEnergyIncr::start(LinearSOE& theSOE)
{
    currentIter = 1;
    norms.Zero();
    norm0 = 0.0;

    return 0;
}


int CTestRelativeEnergyIncr::getNumTests()
{
    return currentIter;
}


int CTestRelativeEnergyIncr::getMaxNumTests()
{
    return maxNumIter;
}


double CTestRelativeEnergyIncr::getRatioNumToMax()
{
    double div = maxNumIter;
    return currentIter/div;
}


const Vector&
CTestRelativeEnergyIncr::getNorms()
{
    return norms;
}


int CTestRelativeEnergyIncr::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    static Vector x(4);
    x(0) = tol;
    x(1) = maxNumIter;
    x(2) = printFlag;
    x(3) = nType;
    res = theChannel.sendVector(this->getDbTag(), cTag, x);
    if (res < 0)
        opserr << "CTestRelativeEnergyIncr::sendSelf() - failed to send data\n";

    return res;
}


int
CTestRelativeEnergyIncr::recvSelf(int cTag, Channel &theChannel,
    FEM_ObjectBroker &theBroker)
{
    int res = 0;
    static Vector x(4);
    res = theChannel.recvVector(this->getDbTag(), cTag, x);

    if (res < 0) {
        opserr << "CTestRelativeEnergyIncr::sendSelf() - failed to send data\n";
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
    }

    return res;
}
