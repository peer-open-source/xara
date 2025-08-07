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
#ifndef Newmark_h
#define Newmark_h

// Written : fmk 
// Created : 11/98
// Modified: 02/05 ahs
// Revision: A
//
// Description: This file contains the class definition for Newmark.
// Newmark is an algorithmic class for performing a transient analysis
// using the Newmark integration scheme.
//
// What: "@(#) Newmark.h, revA"

#include <TransientIntegrator.h>
#include <Vector.h>

class DOF_Group;
class FE_Element;

class Newmark : public TransientIntegrator
{
public:
    // constructors
    Newmark(int classTag=INTEGRATOR_TAGS_Newmark);
    Newmark(double gamma, double beta, 
            int uflag=1,                   // choose which "u"nknown is solved for: d, v or a
            int iflag=3,                   // choose how to "i"nitialize the unknown: Dd=0, Dv=0 or Da=0
            bool aflag=false,
            int classTag=INTEGRATOR_TAGS_Newmark);
    ~Newmark();
    
    // methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object.
    int formEleTangent(FE_Element *theEle)  final;
    int formNodTangent(DOF_Group *theDof)   final;
    int formEleResidual(FE_Element* theEle) final;
    int formNodUnbalance(DOF_Group* theDof) final;
    
    int domainChanged() final;    
    int newStep(double deltaT) final;
    int revertToLastStep() final;
    int update(const Vector &deltaU) final;

    double getCFactor();

    const Vector &getVel();
    
    int sendSelf(int commitTag, Channel &) final;
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &) final;
    
    void Print(OPS_Stream &s, int flag) final;        
    
    // AddingSensitivity:BEGIN //////////////////////////////////
    int revertToStart() override;
    int saveSensitivity   (const Vector &v, int gradNum, int numGrads);
    int computeSensitivities() override;

protected:

private:
    int commitSensitivity (int gradNum, int numGrads);  
    int formSensitivityRHS(int gradNum);
    int formIndependentSensitivityRHS();

    enum Unknown {
      Displacement=1,
      Velocity=2,
      Acceleration=3
    };
    int unknown;                    // flag indicating whether displ(1), vel(2) or accel(3) increments
    int unknown_initialize = 1;     //

    double gamma;
    double beta;

    double c1, c2, c3;              // some constants we need to keep
    Vector *Ut, *Utdot, *Utdotdot;  // response quantities at time t
    Vector *U, *Udot, *Udotdot;     // response quantities at time t+deltaT
    bool determiningMass;           // flag to check if just want the mass contribution

    // Adding Sensitivity
    int sensitivityFlag;
    int gradNumber;
    Vector *massMatrixMultiplicator;
    Vector *dampingMatrixMultiplicator;
    int assemblyFlag;
    Vector independentRHS;
    Vector dUn, dVn, dAn;
    //////////////////////

};

#endif
