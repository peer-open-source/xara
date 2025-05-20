//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
//
//
// Written: cmp
// Created: May 2025
//
#ifndef BWBF_h
#define BWBF_h


#include <UniaxialMaterial.h>


static inline double 
signum(double value)
{
  if (value > 0.0)
    return  1.0;
  else
    return -1.0;
}

class BWBF : public UniaxialMaterial
{
  public:

    BWBF(int tag, 
        double E,
        double Fy,
        double alpha,
        double n,
        double beta,
        //
        double delta_a,
        double delta_v,
        double delta_n,
        //
        double pinch_slope,
        double pinch_slip,
        double pinch_start,
        double pinch_rate,
        double pinch_size,
        double pinch_lamda,
        //
        double tolerance,
        int maxNumIter);
    ~BWBF();

    const char *getClassType() const {return "BWBF";}

    UniaxialMaterial *getCopy();

    int setTrialStrain(double strain, double strainRate = 0.0);
    int commitState();
    int revertToLastCommit();    
    int revertToStart();

    double getStrain();          
    double getStress();
    double getTangent();
    double getInitialTangent();

    int sendSelf(int commitTag, Channel &);  
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);    
    void Print(OPS_Stream &s, int flag);


    int setParameter (const char **argv, int argc, Parameter &);
    int updateParameter          (int parameterID, Information &);
    int activateParameter        (int parameterID);
    // double getStressSensitivity     (int gradIndex, bool conditional);
    // double getStrainSensitivity     (int gradIndex);
    // double getTangentSensitivity    (int gradIndex);
    // double getDampTangentSensitivity(int gradIndex);
    // double getRhoSensitivity        (int gradIndex);
    // int    commitSensitivity        (double strainGradient, int gradIndex, int numGrads);
    // double  getInitialTangentSensitivity(int gradIndex);

  private:
    
    double wen(double z, double Psi, double A, double nu) {
      return (A - pow(fabs(z),n)*Psi*nu)*Ko/Fy;
    }

    double dwen(double z, double Psi, double A, double nu, 
                double dz, double dA, double dnu) {
      double pow1 = (z==0.0)? 0.0 : pow(fabs(z), (n-1));
      return (dA
             -n*pow1*signum(z)*Psi*nu*dz
             - std::pow(std::fabs(z), n)*Psi*dnu)*Ko/Fy;
    }

    // Material parameters
    double Ko;
    double Fy;
    double alpha;
    double n;
    double gamma;
    double beta;
    double Ao;

    double delta_a, delta_v, delta_n;

    
    // History variables (trial and committed)
    struct {
      double strain;
      double energy;
      double tangent;
      double z;
    } pres, past;
    // double Tz, Cz;
    // double Te, Ce;

    // Other variables
    double Tstress, Ttangent;
    
    double tolerance;
    int maxNumIter;

    int parameterID;


    struct Pinch {

      Pinch(double q, double zetas, double p, double psi, double delta_psi, double lamda) 
        : q(q), zetas(zetas), p(p), psi(psi), delta_psi(delta_psi), lamda(lamda) {
        pres.zeta_1 = zetas;
        pres.zeta_2 = 0.0;
      }

      double update(double sgn, double z, double e, double u) {
        pres.e = e;
        pres.u = u;
        pres.z = z;
        pres.sgn = sgn;

        pres.zeta_1 = -zetas*std::expm1(-p*e);
        pres.zeta_2 = (psi + delta_psi*e)*(lamda + pres.zeta_1);
        if (std::fabs(pres.zeta_1) < 1e-10 || std::fabs(pres.zeta_2) < 1.0e-10) {
          return 1.0;
        }
        return 1.0 - pres.zeta_1*std::exp(-std::pow(z*sgn - q*u, 2)/(pres.zeta_2*pres.zeta_2));
      }

      double tangent(double de, double ue, double dz) {
        double Tzeta1 = pres.zeta_1; // *(1.0 - std::exp(-p*e));
        double Tzeta2 = pres.zeta_2; // (psi + delta_psi*e)*(lamda + Tzeta1);
      
        if (std::fabs(pres.zeta_1) < 1e-10 || std::fabs(pres.zeta_2) < 1.0e-10) {
          return 0.0;
        }

        double e = pres.e;
        double u = pres.u;
        double z = pres.z;
        double sgn = pres.sgn;

        double m = (z*sgn - q*u)/pres.zeta_2;
        double a3 = std::exp(-std::pow(m,2));

        double dzeta1 = zetas*p*std::exp(-p*e);
        double dzeta2 = psi*dzeta1 + lamda*delta_psi + delta_psi*e*dzeta1 + delta_psi*Tzeta1;


        double me = -z*sgn/(Tzeta2*Tzeta2)*dzeta2 - q*(ue*Tzeta2 - u*dzeta2)/(Tzeta2*Tzeta2);
        double he = -dzeta1*a3 + Tzeta1*2*m*me*a3;

        double dh = he*de;

        if (dz) {
          double mz = sgn/Tzeta2;
      
          double hze = Tzeta1*2*m*mz*a3;
          dh += hze*dz;
        } 
        return dh;
      }
      double q;
      double zetas;
      double p;
      double psi;
      double delta_psi;
      double lamda;

    private:
      struct {
        double u,e,z, sgn;
        double zeta_1, zeta_2;
      } pres;
    } pinch;
};


#endif



