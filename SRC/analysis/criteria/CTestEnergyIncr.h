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
// Purpose: This file contains the class definition for CTestEnergyIncr.
// A CTestEnergyIncr object tests for convergence using the energy increment,
// which is 0.5 times the absolute value of the product of the rhs and
// the solution vector of the LinearSOE.
//
// Written: fmk
// Date: 09/98
// Modified: 05/05 ahs
//
#ifndef CTestEnergyIncr_h
#define CTestEnergyIncr_h

#include <ConvergenceTest.h>
class LinearSOE;


class CTestEnergyIncr: public ConvergenceTest
{
public:
    CTestEnergyIncr();
    CTestEnergyIncr(double tol, int maxNumIter, int printFlag, int normType =2, double maxTol = OPS_MAXTOL);

    ~CTestEnergyIncr();

    const char* getClassType() const override { return "EnergyIncr"; }

    ConvergenceTest *getCopy(int iterations) override;

    void setTolerance(double newTol);

    int start(LinearSOE&) override;
    int test(LinearSOE&) override;

    int getNumTests() override;
    int getMaxNumTests() override;
    double getRatioNumToMax() override;
    const Vector &getNorms() override;

    int sendSelf(int commitTag, Channel &) override;
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &) override;

private:
    double tol;         // the tol on the energy used to test for convergence
    double maxTol;

    int maxNumIter;     // max number of iterations
    int currentIter;    // number of times test() has been invokes since last start()
    int printFlag;      // a flag indicating if to print on test
    int nType;          // type of norm to use (1-norm, 2-norm, p-norm, max-norm)

    Vector norms;       // vector to hold the norms
};

#endif
