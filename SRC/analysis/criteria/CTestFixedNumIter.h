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
// Purpose: This file contains the class definition for CTestFixedNumIter.
// A CTestFixedNumIter object performs a fixed number of iterations without
// testing for convergence. This test is useful for hybrid simulation where
// the residual error is corrected for.
//
// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 09/05
// Revision: A
//
#ifndef CTestFixedNumIter_h
#define CTestFixedNumIter_h

#include <ConvergenceTest.h>
class LinearSOE;


class CTestFixedNumIter: public ConvergenceTest
{
public:
    // constructors
    CTestFixedNumIter();
    CTestFixedNumIter(int maxNumIter, int printFlag, int normType=2);

    ~CTestFixedNumIter();

    const char* getClassType() const override { return "FixedNumIter"; }

    ConvergenceTest *getCopy(int iterations);

    int test(LinearSOE&) override;
    int start(LinearSOE&) override;

    int getNumTests() override;
    int getMaxNumTests() override;
    double getRatioNumToMax() override;
    const Vector &getNorms() override;

private:
    int maxNumIter;     // max number of iterations
    int currentIter;    // number of times test() has been invokes since last start()
    int printFlag;      // a flag indicating if to print on test
    int nType;          // type of norm to use (1-norm, 2-norm, p-norm, max-norm)

    Vector norms;       // vector to hold the norms
};

#endif
