//
// Written: Quan Gu and Zhijian Qiu 
// Created: 2013.7 
// Reference:JP Conte, MK. Jagannath, \Seismic relibility analysis of concrete 
// gravity dams, A Report on Research, Rice University, 1995.
// 3D J2 plasticity model with linear isotropic and kinematic hardening
//  
// -------------------


#ifndef SimplifiedJ2_h
#define SimplifiedJ2_h

#include <NDMaterial.h>
#include <Matrix.h>
#include <MatrixND.h>
#include <Vector.h>


class SimplifiedJ2 : public NDMaterial
{
 public:
  
  SimplifiedJ2 (int tag, 
                int nd,
                double G,
                double K,
                double sigmaY0,
                double H_kin,
                double H_iso);
  SimplifiedJ2();
  
  SimplifiedJ2(const SimplifiedJ2 &);
  virtual ~SimplifiedJ2 ();
  
  const char *getClassType(void) const {return "SimplifiedJ2";}
  const char *getType(void) const { return "ThreeDimensional";}

  int setTrialStrain (const Vector &strain);
  int setTrialStrain(const Vector &v, const Vector &r);
  int setTrialStrainIncr(const Vector &v);
  int setTrialStrainIncr(const Vector &v, const Vector &r);
  
  // Calculates current tangent stiffness.
  const Matrix &getTangent ();
  const Matrix &getInitialTangent ();
  
  // Calculates the corresponding stress increment (rate), for a given strain increment. 
  const Vector &getStress ();
  const Vector &getStrain ();
  const Vector &getCommittedStress ();
  const Vector &getCommittedStrain ();
  
  /*int setTrialStrain (const Tensor &v) {return 0;}
    int setTrialStrain (const Tensor &v, const Tensor &r) {return 0;}
    int setTrialStrainIncr (const Tensor &v) {return 0;}
    int setTrialStrainIncr (const Tensor &v, const Tensor &r) {return 0;}
  */
  
  int commitState ();
  int revertToLastCommit ();
  int revertToStart();
  
  NDMaterial *getCopy ();
  NDMaterial *getCopy (const char *code);
  
  
  int sendSelf(int commitTag, Channel &);  
  int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);    
  
  Response *setResponse (const char **argv, int argc, OPS_Stream &s);
  int getResponse (int responseID, Information &matInformation);
  void Print(OPS_Stream &s, int flag =0);
  
  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int responseID, Information &eleInformation);
  
 private:
  
  int plastIntegrator();
  
  int ndm;
  double G;
  double K;
  double sigmaY0;
  double H_kin;
  double H_iso;
  
  
  // --- trial variables
  Vector stress;
  Vector strain;
  
  
  Vector plastStrainDev;
  Vector backStress;
  double sigmaY;
  double cumPlastStrainDev;
  
  // ---- committed variables
  
  Vector Cstress;
  Vector Cstrain;
  Vector CplastStrainDev;
  Vector CbackStress;
  double CsigmaY;
  double CcumPlastStrainDev;
  
  // --
  double lambda;  // positive plastStrainDevInc
  
  
  OpenSees::MatrixND<6,6> tangent;
  Matrix theTangent; 
  
  
  // ---  define classwide variables
  
  static Vector tmpVector;
  static Matrix tmpMatrix;
};

#endif
