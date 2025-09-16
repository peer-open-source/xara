//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, Claudio M. Perez
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the class definition for GS4.
// GS4 is an algorithmic class for performing a transient analysis
// using the GS4 integration scheme.
//
// Written : cmp
// Created : 10/2024
// Adapted from GeneralizedNewmark.cpp
//
#pragma once

#include <TransientIntegrator.h>
#include <Vector.h>
#include <array>

class DOF_Group;
class FE_Element;

class GS4 : public TransientIntegrator
{
public:
    GS4(double gamma, double beta, 
        double alphaF, double alpahM,
        int uflag=1,                   // choose which "u"nknown is solved for: d, v or a
        int iflag=3,                   // choose how to "i"nitialize the unknown: Dd=0, Dv=0 or Da=0
        bool aflag=false);

    ~GS4();

    //
    // Integrator
    //
    // methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object.
    int formEleTangent(FE_Element*)  final;
    int formNodTangent(DOF_Group*)   final;
    int formEleResidual(FE_Element*) final;
    int formNodUnbalance(DOF_Group*) final;

    //
    // IncrementalIntegrator
    //

    // Sensitivity
    int formSensitivityRHS(int gradNum);
    int updateGradient (const Vector &v, int gradNum, int numGrads);
    int commitGradient (int gradNum, int numGrads) {return -1;};
    int computeSensitivities();

    //
    // TransientIntegrator
    //
    int newStep(double deltaT) final;
    //
    int domainChanged() final;
    //
    // IncrementalIntegrator
    //
    int revertToLastStep() final;
    int revertToStart();
    int update(const Vector &deltaU) final;

    double getCFactor();
    const Vector &getVel();

    // MovableObject
    int sendSelf(int commitTag, Channel &) override;
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &) override;
    
    void Print(OPS_Stream &, int flag) final;

private:
    enum Unknown {
      Displacement=1,
      Velocity=2,
      Acceleration=3
    };

    struct GammaScheme {
      double g[3];
      double G[3][3];
    } G1, G2;

    int unknown;                    // flag indicating whether displ(1), vel(2) or accel(3) increments
    int unknown_initialize = 1;     //

    static int setConstants(int flag,
                            double dt, double gamma, double beta,
                            const std::array<double,3> &alpha,
                            GammaScheme& scheme);
    double gamma;
    double beta;
    double alphaF;
    double alphaM;
    std::array<double,3> alpha;

    int step;                       // track step number to initialize accelerations
    double dt;                      // store time step to determine step number
    double cu, cv, ca;              // some constants we need to keep
    Vector *Uo, *Vo, *Ao;           // solution at time t
    Vector *Ua, *Va, *Aa;           // solution at time t+alpha
    Vector *Un, *Vn, *An;           // solution at time t+deltaT
    bool determiningMass;           // flag to check if just want the mass contribution

    // Sensitivity
    int isSensitivityResidual;
    int gradNumber;
    Vector *dAa;
    Vector *dVa;
    int assemblyFlag;
    Vector dUn, dVn, dAn;
};
