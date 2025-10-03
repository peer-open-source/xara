//
// Created by Alex Hartloper on 09.07.18.
//

#ifndef CPP_UVC_MA_H
#define CPP_UVC_MA_H

#include <vector>
#include "NDMaterial.h"
#include "Matrix.h"
#include "Vector.h"
#include "MatrixND.h"
#include "VectorND.h"
using namespace OpenSees;


/* ------------------------------------------------------------------------ */

class UVCmultiaxial : public NDMaterial {

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

public:
   UVCmultiaxial();
  ~UVCmultiaxial();
  UVCmultiaxial(int tag,
                double E, double poissonRatio, 
                double sy0, 
                double qInf, double b,
                double dInf, double a, 
                const std::vector<double>& cK, 
                const std::vector<double>& gammaK);


  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

public:
  const char*
  getClassType() const override
  {
    return "UVCmultiaxial";
  }

  const char*
  getType() const override
  {
    return "ThreeDimensional";
  }

  //! Returns the number of vector components
  int
  getOrder() const override
  {
    return 6;
  }

  //! Calculates the trial strain and stress, provided the total strain
  int setTrialStrain(const Vector& v);
  int setTrialStrain(const Vector& v, const Vector& r);

  //! Calculates the trial strain and stress, provided the strain increment
  int setTrialStrainIncr(const Vector& v);
  int setTrialStrainIncr(const Vector& v, const Vector& r);


  //! Returns the trial strain
  const Vector& getStrain();

  //! Returns the trial stress
  const Vector& getStress();

  //! Returns the trial elastoplastic tangent modulus
  const Matrix& getTangent();

  //! Returns the tangent modulus in the undeformed configuration
  const Matrix& getInitialTangent();

  //! Returns the mass density of the material - zero mass assumed
  double
  getRho()
  {
    return 0.;
  }

  //! Sets the converged state to be the current trial state
  int commitState();

  //! Sets the trial state to be the converged state
  int revertToLastCommit();

  //! Sets the converged state to the undeformed configuration
  int revertToStart();

  //! Returns a copy of the material in the current state
  NDMaterial* getCopy();

  //! Returns a copy of the material without copying the state variables
  NDMaterial* getCopy(const char* code);

  //! todo: fill out
  int sendSelf(int commitTag, Channel&);

  //! todo: fill out
  int recvSelf(int commitTag, Channel&, FEM_ObjectBroker&);

  //! Adds the print information to the stream
  void Print(OPS_Stream& s, int flag);

private:
  //! Determines the trial stress for the given strain increment
  int returnMapping();

  //! Sets the elastoplastic tangent modulus based on the trial state
  void calculateStiffness(double plasticMultiplier, double stressRelativeNorm, Vector alphaDiff);

  //! Returns the equivalent dot product of a 2nd order symmetric tensor
  double dotprod6(Vector v1, Vector v2);

  //! Calculates the elastic stiffness matrix
  void calculateElasticStiffness();

  //! Returns the current yield stress
  double calculateYieldStress();

  //! Returns the isotropic hardening modulus
  double calculateIsotropicModulus();

  //! Returns the current eK value
  double calculateEk(unsigned int i);

  /* ------------------------------------------------------------------------ */
  /* Members                                                                  */
  /* ------------------------------------------------------------------------ */

private:
  // Parameters
  static constexpr unsigned int N_BASIC_PARAMS     = 5;
  static constexpr unsigned int N_PARAM_PER_BACK   = 2;
  static constexpr double RETURN_MAP_TOL           = 1.0e-10;
  static constexpr unsigned int MAXIMUM_ITERATIONS = 1000;
  static constexpr unsigned int N_DIRECT           = 3;
  static constexpr unsigned int N_DIMS             = 6;

  // Material properties, set by the constructor
  double elasticModulus;
  double shearModulus;
  double bulkModulus;
  double poissonRatio;
  double initialYield;
  double qInf;
  double bIso;
  double dInf;
  double aIso;
  Matrix stiffnessInitial;
  Matrix elasticMatrix;
  std::vector<double> cK;
  std::vector<double> gammaK;
  unsigned int nBackstresses;

  // Internal variables
  Vector strainConverged;
  Vector strainTrial;
  Vector strainPlasticConverged;
  Vector strainPlasticTrial;
  double strainPEqConverged; // Equivalent plastic strain
  double strainPEqTrial;
  Vector stressConverged;
  Vector stressTrial;
  std::vector<VectorND<6>> alphaKConverged;
  std::vector<VectorND<6>> alphaKTrial;
  Matrix stiffnessConverged;
  Matrix stiffnessTrial;
  VectorND<6> flowNormal;
  bool plasticLoading;
};

/* ------------------------------------------------------------------------ */

#endif //CPP_UVC_MA_H
