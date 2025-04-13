//
//
//
#define _USE_MATH_DEFINES
#include <cmath>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <limits>
#include <Logging.h>
#include "BoucWenMG.h"

BoucWenMG::BoucWenMG(int tag, const Params& params)
  : UniaxialMaterial(tag, 0), params(params)
{
    cState.stiffness = params.k0;
    tState.stiffness = params.k0;
    revertToStart();
}

UniaxialMaterial*
BoucWenMG::getCopy()
{
    BoucWenMG* theCopy = new BoucWenMG(getTag(), params);
    return theCopy;
}

int
BoucWenMG::revertToStart() {
    cState = State();
    tState = State();
    cState.stiffness = params.k0;
    tState.stiffness = params.k0;
    cState.stressY   = params.sy0;
    tState.stressY   = params.sy0;
    return 0;
}

int
BoucWenMG::revertToLastCommit()
{
    tState = cState;
    return 0;
}

int
BoucWenMG::commitState() {
    cState = tState;
    return 0;
}

double
BoucWenMG::getStrain()
{
    return cState.strain;
}

double 
BoucWenMG::getInitialTangent()
{
    return params.k0;
}

double
BoucWenMG::getStress()
{
    return tState.stress;
}

double
BoucWenMG::getTangent()
{
    return tState.stiffness;
}


int 
BoucWenMG::setTrialStrain(double tstrain)
{
    // (1) Model Parameters
    double kTrue = params.k0;
    double syTrue = params.sy0;
    double eyTrue = syTrue / kTrue;

    double k0  = 1.0; // Normalized stiffness
    double sy0 = 1.0; // Normalized yield stress

    // Shape parameters
    double eta1 = params.eta;
    double eta2 = params.eta2;
    double n = params.n;
    double alpha = params.alpha;

    // Pinching parameters
    double sig = params.sig;
    double lam = params.lam;
    double mup = params.mup;
    double sigp = params.sigp;
    double rsmax = params.rsmax;

    // Degradation parameters
    double alpha1 = params.alpha1;
    double alpha2 = params.alpha2;
    double betam1 = params.betam1;

    // (2) Load current state variables
    double stressEl = cState.stressEl;
    double stressSt = cState.stressSt;
    double stressY = cState.stressY;
    double stressSig = sig * stressY;
    double stressBar = lam * stressY;

    double kEl = alpha * k0;
    double rk = cState.rk;
    double rkmin = cState.rkmin;
    double h = cState.h;
    double emaxPos = cState.emaxPos;
    double emaxNeg = cState.emaxNeg;
    double emax = cState.emax;
    double rs = cState.rs;

    int regime = cState.regime;
    auto& thetaMaxPos = cState.thetaMaxPos;
    auto& thetaMaxNeg = cState.thetaMaxNeg;
    auto& momentMaxPos = cState.momentMaxPos;
    auto& momentMaxNeg = cState.momentMaxNeg;
    double k2 = cState.k2;

    // Normalize responses
    double cStress = cState.stress / syTrue;
    double strainI = cState.strain / eyTrue;
    double strainIp1 = tstrain / eyTrue;
    double dStrain = strainIp1 - strainI;

    // Define regime
    if (cStress * dStrain < 0 && regime == 0) {
        regime = 1;
        if (cStress > 0) {
            thetaMaxPos = {strainI, strainI};
            momentMaxPos = {cStress, cStress};
        } else {
            thetaMaxNeg = {strainI, strainI};
            momentMaxNeg = {cStress, cStress};
        }
    }

    if (regime == 1 && cStress > 0 && dStrain > 0 && strainIp1 > 0) {
        regime = 2;
        thetaMaxPos[0] = strainI;
        momentMaxPos[0] = cStress;
        k2 = (thetaMaxPos[1] - thetaMaxPos[0]) != 0 ? 
             (momentMaxPos[1] - momentMaxPos[0]) / (thetaMaxPos[1] - thetaMaxPos[0]) : k0;
    
    } else if (regime == 1 && cStress < 0 && dStrain < 0 && strainIp1 < 0) {
        regime = 2;
        thetaMaxNeg[0] = strainI;
        momentMaxNeg[0] = cStress;
        k2 = (thetaMaxNeg[1] - thetaMaxNeg[0]) != 0 ? 
             (momentMaxNeg[1] - momentMaxNeg[0]) / (thetaMaxNeg[1] - thetaMaxNeg[0]) : k0;
    }

    double k;
    double Tz = stressSt;
    if (regime == 0 || regime == 1) {
        double TzOld = Tz;
        double TzNew = Tz + 0.01;
        int count = 0;
        int maxIter = 500;
        double kh,ks;

        while (std::fabs(TzOld - TzNew) > 1e-6 && count < maxIter) {
            double A = 1.0 - (eta1 * std::copysign(1.0, Tz * dStrain) + eta2) * std::pow(std::fabs(Tz / stressY), n);
            kh = (rk - alpha) * k0 * A;
            double s = std::max(rs * (emaxPos - emaxNeg), 0.0001);
            double B = std::exp(-0.5 * std::pow((Tz - stressBar * std::copysign(1.0, dStrain)) / stressSig, 2));

            if (B < 1e-20)
                B = 1e-20;

            ks = std::min(std::pow((1.0 / std::sqrt(2.0 * M_PI)) * (s / stressSig) * B, -1), 1000.0);
            double kr = (kh * ks) / (kh + ks);
            double f = kr * dStrain - Tz + stressSt;

            double A_z = - (eta1 * std::copysign(1.0, Tz * dStrain) + eta2) * n * 
                        std::pow(std::fabs(Tz / stressSig), n - 1) * std::copysign(1.0, Tz / stressSig);
            double B_z = - (Tz - stressBar * std::copysign(1.0, dStrain)) * B / (stressSig * stressSig);
            double ks_z = - (std::sqrt(2.0 * M_PI) * stressSig / s) * (B * B) * B_z;
            double kh_z = (rk - alpha) * k0 * A_z;
            double kr_z = (kh_z * ks * ks + ks_z * kh * kh) / (kh + ks) / (kh + ks);
            double f_z = -1.0 + kr_z * dStrain;

            TzNew = Tz - f / f_z;
            TzOld = Tz;
            Tz = TzNew;
            count++;

            if (std::isnan(Tz))
                Tz = 0.001;
            
            if (count == maxIter) {
                opserr << "Failed to converge in BoucWenMG::setTrialState()";
                convFlag = false;
            }
        }

        if (std::fabs(kh + ks) > std::numeric_limits<double>::epsilon()) {
            k = kEl + (kh * ks) / (kh + ks);  // Total stiffness of parallel/series spring system
        } else {
            k = k0;
            opserr << "Invalid stiffness calculation due to division by zero!";
        }
    }
    else {
        // regime != 0 or 1
        k = k2;
    }

    // Compute normalized stress and elastic/plastic stress
    double tStress = k * dStrain + cStress; // This is still a normalized stress value
    stressEl = stressEl + kEl * dStrain;   // Elastic normalized stress
    stressSt = tStress - stressEl;         // Elastic-plastic stress

    // How to go back to regime 0?
    if (regime == 2) {
        if (tStress > 0 && dStrain > 0 && tStress > momentMaxPos[1]) { 
            // If on Q1 and moment exceed pivot, then go back to regime 0
            regime = 0;
            tStress = 0.1 * tStress + 0.9 * momentMaxPos[1];
        } else if (tStress > 0 && dStrain < 0 && tStress < momentMaxPos[0]) {
            regime = 1;
            momentMaxPos[0] = cStress;
        } else if (tStress < 0 && dStrain < 0 && tStress < momentMaxNeg[1]) {
            regime = 0;
            tStress = 0.1 * tStress + 0.9 * momentMaxNeg[1];
        } else if (tStress < 0 && dStrain > 0 && tStress > momentMaxNeg[0]) {
            regime = 1;
            momentMaxNeg[0] = cStress;
        }
    }

    if (tStress * strainIp1 < 0) {
        regime = 0;
    }

    // Stiffness Degradation
    double rkTrial = (std::abs(tStress) + alpha1 * stressY) / 
                    (k0 * std::abs(strainI) + alpha1 * stressY);

    rkmin = std::min(rkTrial, rkmin);
    rk = rkTrial + (1 - alpha2) * (rkmin - rkTrial);

    // Strength Degradation
    double dh = tStress * dStrain;
    h = h + dh;
    stressY = sy0 / (1 + betam1 * h);
    emaxPos = std::max(emaxPos, strainI);
    emaxNeg = std::min(emaxNeg, strainI);

    // Msig = sig * My
    // Mbar = lam * My

    if (std::abs(strainI) - emax > 0) {
        double drs = (1 / (std::sqrt(2 * M_PI) * sigp)) * 
                    std::exp(-0.5 * std::pow((emax - mup) / sigp, 2)) * 
                    (std::abs(strainI) - emax);
        rs = rs + drs * rsmax;
        emax = std::abs(strainI);
    }

    // De-normalize the state parameters
    k = k * kTrue;
    tStress = tStress * syTrue;
    double tStrain = strainIp1 * eyTrue;

    // Update trial state
    tState.strain = tStrain;
    tState.stress = tStress;
    tState.stiffness = k;
    tState.z = Tz;
    tState.h = h;
    tState.stressEl = stressEl;
    tState.stressSt = stressSt;
    tState.stressY = stressY;
    tState.rs = rs;
    tState.rk = rk;
    tState.rkmin = rkmin;
    tState.emaxPos = emaxPos;
    tState.emaxNeg = emaxNeg;
    tState.emax = emax;
    tState.regime = regime;
    tState.thetaMaxPos = thetaMaxPos;
    tState.thetaMaxNeg = thetaMaxNeg;
    tState.momentMaxPos = momentMaxPos;
    tState.momentMaxNeg = momentMaxNeg;
    tState.k2 = k2;

    return convFlag? 0 : -1;
}


int
BoucWenMG::sendSelf(int, Channel&)
{
    // TODO
    return -1;
}

int
BoucWenMG::recvSelf(int, Channel&, FEM_ObjectBroker&)
{
    // TODO
    return -1;
}

void 
BoucWenMG::Print(OPS_Stream& s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
      s << OPS_PRINT_JSON_MATE_INDENT << "{";
      s << "\"name\": \"" << this->getTag() << "\", ";
      s << "\"type\": \"" << this->getClassType() << "\", ";
      s << "\"eta\": "    << params.eta << ", ";
      s << "\"eta2\": "   << params.eta2 << ", ";
      s << "\"k0\": "     << params.k0 << ", ";     // True elastic stiffness
      s << "\"sy0\": "    << params.sy0 << ", ";    // True yield stress
      s << "\"sig\": "    << params.sig << ", ";
      s << "\"lam\": "    << params.lam << ", ";
      s << "\"mup\": "    << params.mup << ", ";
      s << "\"sigp\": "   << params.sigp << ", ";
      s << "\"rsmax\": "  << params.rsmax << ", ";
      s << "\"n\": "      << params.n << ", ";
      s << "\"alpha\": "  << params.alpha << ", ";
      s << "\"alpha1\": " << params.alpha1 << ", ";
      s << "\"alpha2\": " << params.alpha2 << ", ";
      s << "\"betam1\": " << params.betam1;
      s << "}";
    }
    return;
}
