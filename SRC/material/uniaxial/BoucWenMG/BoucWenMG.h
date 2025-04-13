#ifndef DEGBWMATERIALMOD_H
#define DEGBWMATERIALMOD_H

#include <array>
#include <vector>
#include <UniaxialMaterial.h>
class Channel;
class FEM_ObjectBroker;

class BoucWenMG: public UniaxialMaterial {
public:
    struct Params;
    BoucWenMG(int tag, const Params& params);

    // UniaxialMaterial
    int revertToStart();
    int revertToLastCommit();
    int setTrialStrain(double tstrain);
    int setTrialStrain(double tstrain, double) {
        return setTrialStrain(tstrain);
    }
    int commitState();
    double getStress();
    double getStrain();
    double getInitialTangent();
    double getTangent();
    UniaxialMaterial* getCopy();
    const char*getClassType() const {
      return "BoucWenMG";
    }

    // MovableObject
    int sendSelf(int, Channel&);
    int recvSelf(int, Channel&, FEM_ObjectBroker&);

    // TaggedObject
    void Print(OPS_Stream&, int);

    // Personal
    struct Params {
        double eta;
        double eta2;
        double k0;      // True elastic stiffness
        double sy0;     // True yield stress
        double sig;
        double lam;
        double mup;
        double sigp;
        double rsmax;
        double n;
        double alpha;
        double alpha1;
        double alpha2;
        double betam1;
    };

private:
    struct State {
        double strain = 0.0;
        double stress = 0.0;
        double stiffness = 0.0;
        double z = 0.0;
        double h = 0.0;
        double stressEl = 0.0;
        double stressSt = 0.0;
        double stressY = 1.0;
        double rs = 0.0;
        double rk = 1.0;
        double rkmin = 1.0;
        double emaxPos = 0.0;
        double emaxNeg = 0.0;
        double emax = 0.0;
        int regime = 0;
        std::array<double, 2> thetaMaxPos  = {0.0, 0.0};
        std::array<double, 2> thetaMaxNeg  = {0.0, 0.0};
        std::array<double, 2> momentMaxPos = {0.0, 0.0};
        std::array<double, 2> momentMaxNeg = {0.0, 0.0};
        double k2 = 1.0;
    };

    Params params;
    State cState;
    State tState;
    bool convFlag = true;

};


#endif // DEGBWMATERIALMOD_H

