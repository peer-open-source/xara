//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the class definition for GeneralizedNewmark.
// GeneralizedNewmark is an algorithmic class for performing a transient analysis
// using the GeneralizedNewmark integration scheme.
//
// Written : cmp
// Created : 06/2024
// Adapted from Newmark.cpp
//
#ifndef GeneralizedNewmark_h
#define GeneralizedNewmark_h

#include <TransientIntegrator.h>
#include <Vector.h>

class DOF_Group;
class FE_Element;

class GeneralizedNewmark : public TransientIntegrator
{
public:
    GeneralizedNewmark(double gamma, double beta, 
                       double alphaF, double alpahM,
                       int uflag=1,                   // choose which "u"nknown is solved for: d, v or a
                       int iflag=3,                   // choose how to "i"nitialize the unknown: Dd=0, Dv=0 or Da=0
                       bool aflag=false);

    ~GeneralizedNewmark();
    
    // methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object.
    int formEleTangent(FE_Element *theEle)  final;
    int formNodTangent(DOF_Group *theDof)   final;
    int formEleResidual(FE_Element* theEle) final;
    int formNodUnbalance(DOF_Group* theDof) final;

    int domainChanged();    
    int newStep(double deltaT);
    virtual int update(const Vector &deltaU);
    int revertToLastStep();

    double getCFactor();

    const Vector &getVel();

    // MovableObject
    virtual int sendSelf(int commitTag, Channel &);
    virtual int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);
    
    void Print(OPS_Stream &s, int flag);        
    
    // Sensitivity
    int revertToStart();
    int saveSensitivity   (const Vector &v, int gradNum, int numGrads);
    int commitSensitivity (int gradNum, int numGrads);  
    int computeSensitivities();


private:
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
    double alphaF;
    double alphaM;

    int step;                       // track step number to initialize accelerations
    double dt;                      // store time step to determine step number
    double cu, cv, ca;              // some constants we need to keep
    Vector *Uo, *Vo, *Ao;           // solution at time t
    Vector *Ua, *Va, *Aa;           // solution at time t+alpha
    Vector *Un, *Vn, *An;           // solution at time t+deltaT
    bool determiningMass;           // flag to check if just want the mass contribution

    // Sensitivity
    int sensitivityFlag;
    int gradNumber;
    Vector *dAa;
    Vector *dVa;
    int assemblyFlag;
    Vector independentRHS;
    Vector dUn, dVn, dAn;
};

#endif
