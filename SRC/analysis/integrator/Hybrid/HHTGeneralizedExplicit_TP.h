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
// Description: This file contains the class definition for HHTGeneralizedExplicit_TP.
// HHTGeneralizedExplicit_TP is an algorithmic class for performing a transient analysis
// using the HHTGeneralizedExplicit_TP integration scheme based on the trapezoidal rule.
// Do not use this integrator for hybrid simulation. It updates the element displacements
// twice per time step because it needs the resisiting force at Ut and Upt
//
// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 10/05
// Revision: A
//
#pragma once
#include <TransientIntegrator.h>

class DOF_Group;
class FE_Element;
class Vector;

class HHTGeneralizedExplicit_TP : public TransientIntegrator
{
public:
    // constructors
    HHTGeneralizedExplicit_TP(double rhoB, double alphaF);
    HHTGeneralizedExplicit_TP(double alphaI, double alphaF,
        double beta, double gamma);
    
    // destructor
    ~HHTGeneralizedExplicit_TP();
    
    // method to set up the system of equations
    int formUnbalance(Vector& G) override;
    
    // methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object.
    int formEleTangent(FE_Element *);
    int formNodTangent(DOF_Group *);
    int formEleResidual(FE_Element *);
    int formNodUnbalance(DOF_Group *);
    
    int domainChanged();
    int newStep(double deltaT);
    int revertToLastStep();
    int update(const Vector &aiPlusOne);
    int commit();

    const Vector &getVel() override;

    void Print(OPS_Stream &s, int flag = 0);
    
private:
    double alphaI;
    double alphaF;
    double beta;
    double gamma;
    double deltaT;
    
    int updateCount;                            // method should only have one update per step
    double c1, c2, c3;                          // some constants we need to keep
    double alphaM, alphaD, alphaR, alphaP;      // weighting factors we need to keep
    Vector *Ut, *Utdot, *Utdotdot;              // response quantities at time t
    Vector *U, *Udot, *Udotdot;                 // response quantities at time t + deltaT
    Vector *Put;                                // unbalance at time t
};

