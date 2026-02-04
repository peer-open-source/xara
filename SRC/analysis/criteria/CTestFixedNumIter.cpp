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
// Purpose: This file contains the class implementation of CTestFixedNumIter.
//
// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 09/05
// Revision: A
//
#include <CTestFixedNumIter.h>
#include <Vector.h>
#include <Channel.h>
#include <EquiSolnAlgo.h>
#include <LinearSOE.h>
#include <Logging.h>


CTestFixedNumIter::CTestFixedNumIter()
    : ConvergenceTest(CONVERGENCE_TEST_CTestFixedNumIter),
    maxNumIter(0), currentIter(0), printFlag(0),
    norms(1), nType(2)
{

}


CTestFixedNumIter::CTestFixedNumIter(int maxIter, int printIt, int normType)
    : ConvergenceTest(CONVERGENCE_TEST_CTestFixedNumIter),
    maxNumIter(maxIter), currentIter(0), printFlag(printIt),
    norms(maxNumIter), nType(normType)
{

}


CTestFixedNumIter::~CTestFixedNumIter()
{

}


ConvergenceTest*
CTestFixedNumIter::getCopy(int iterations)
{
  return new CTestFixedNumIter(iterations, this->printFlag, this->nType) ;
}


int CTestFixedNumIter::test(LinearSOE& theSOE)
{

    // check to ensure the algo does invoke start() - this is needed otherwise
    // may never get convergence later on in analysis!
    if (currentIter == 0)  {
        opserr << "WARNING: CTestFixedNumIter::test() - start() was never invoked.\n";
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

    // print the data if required
    if (printFlag & ConvergenceTest::PrintTest)  {
        pstream << LOG_ITERATE << "Iter: " << pad(currentIter);
        pstream << ", EnergyIncr: " << product;
        pstream << " (Norm deltaX: " << x.pNorm(nType) << ", Norm dR: " << b.pNorm(nType) << ")\n";
    }

    if (printFlag & ConvergenceTest::PrintTest02)  {
        pstream << LOG_ITERATE << "Iter: " << pad(currentIter);
        pstream << ", EnergyIncr: " << product;
        pstream << " (Norm deltaX: " << x.pNorm(nType) << ", Norm dR: " << b.pNorm(nType) << ")\n";
        pstream << "\tdeltaX: " << x << "\tdR: " << b;
    }

    //
    // check if the algorithm converged
    //

    // if converged - print & return ok
    if (currentIter == maxNumIter)  {
        if (printFlag & ConvergenceTest::PrintTest || printFlag & ConvergenceTest::PrintTest02)
            pstream << "\n";

        if (printFlag & ConvergenceTest::PrintSuccess)  {
            pstream << LOG_SUCCESS << "Iter: " << pad(currentIter);
            pstream << " last EnergyIncr: " << product;
            pstream << " (Norm deltaX: " << x.pNorm(nType) << ", Norm dR: " << b.pNorm(nType) << ")\n";
        }

        // return the number of times test has been called
        return currentIter;
    }

    // algorithm not yet converged - increment counter and return -1
    else {
        currentIter++;
        return ConvergenceTest::Continue;
    }
}


int
CTestFixedNumIter::start(LinearSOE&)
{
    // set iteration count = 1
    currentIter = 1;
    norms.Zero();
    return 0;
}


int 
CTestFixedNumIter::getNumTests()
{
  return currentIter;
}


int 
CTestFixedNumIter::getMaxNumTests()
{
  return maxNumIter;
}


double
CTestFixedNumIter::getRatioNumToMax()
{
    double div = maxNumIter;
    return currentIter/div;
}


const Vector&
CTestFixedNumIter::getNorms()
{
    return norms;
}


int CTestFixedNumIter::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;
    Vector x(3);
    x(0) = maxNumIter;
    x(1) = printFlag;
    x(2) = nType;
    res = theChannel.sendVector(this->getDbTag(), cTag, x);
    if (res < 0)
        opserr << "CTestFixedNumIter::sendSelf() - failed to send data\n";

    return res;
}


int CTestFixedNumIter::recvSelf(int cTag, Channel &theChannel,
    FEM_ObjectBroker &theBroker)
{
    int res = 0;
    Vector x(3);
    res = theChannel.recvVector(this->getDbTag(), cTag, x);

    if (res < 0) {
        opserr << "CTestFixedNumIter::sendSelf() - failed to send data\n";
        maxNumIter = 25;
        printFlag = 0;
        nType = 2;
    } else  {
        maxNumIter = (int) x(0);
        printFlag = (int) x(1);
        nType = (int) x(2);
        norms.resize(maxNumIter);
    }
    return res;
}
