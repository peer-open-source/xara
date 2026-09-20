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
// Description: This file contains the class definition for LoadControl.
// LoadControl is an algorithmic class for performing a static analysis
// using a load control integration scheme.
// 
// Written: fmk 
// Created: 07/98
// Revision: A
//
#pragma once
//
#include <StaticIntegrator.h>
#include <classTags.h>
#include <control/PredictorControl.h>
#include <Vector.h>

class LinearSOE;
class AnalysisModel;
class FE_Element;
class Vector;
class EquiSolnAlgo;
class ReliabilityDomain;

class LoadControl : public StaticIntegrator
{
  public:
    LoadControl(double deltaLambda, int numIncr, 
                double minLambda, double maxlambda, 
                double exponent = 1.0);

    ~LoadControl();

    virtual int newStep() override;
    int update(const Vector &deltaU) override;
    int commit() override;
    int revertToLastStep() override;
    int domainChanged() override;

    void Print(OPS_Stream &, int flag) override;

    //
    virtual int updateGradient(const Vector &v, int gradNum, int numGrads);
    virtual int commitSensitivity(int gradNum, int numGrads);
    int computeSensitivities() override;
    bool computeSensitivityAtEachIteration() override;

    
  private:
    int formSensitivityRHS(int gradNum);

    double deltaLambda;                      // dlambda at step (i-1)
    double expon;                            // exponent for J(i-1)/Jd
    double specNumIncrStep, numIncrLastStep; // Jd & J(i-1) 
    double dLambdaMin, dLambdaMax;           // min & max values for dlambda at step (i)
    PredictorControl pcontrol;
    
    // Vector R; // residual vector
    // Vector Pu; // Response vector evaluated at u
    // Vector Pr; // Reference load vector
};
