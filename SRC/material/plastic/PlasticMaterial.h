/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// Written: Ed "C++" Love
//
// J2 isotropic hardening material class
//
//  Elastic Model
//  sigma = K*trace(epsilion_elastic) + (2*G)*dev(epsilon_elastic)
//
//  Yield Function
//  phi(sigma,q) = || dev(sigma) ||  - sqrt(2/3)*q(xi)
//
//  Saturation Isotropic Hardening with linear term
//  q(xi) = simga_0 + (sigma_infty - sigma_0)*exp(-delta*xi) + H*xi
//
//  Flow Rules
//  \dot{epsilon_p} =  gamma * d_phi/d_sigma
//  \dot{xi}        = -gamma * d_phi/d_q
//
//  Linear Viscosity
//  gamma = phi / eta  ( if phi > 0 )
//
//  Backward Euler Integration Routine
//  Yield condition enforced at time n+1
//
//  set eta := 0 for rate independent case
//
#ifndef PlasticMaterial_h
#define PlasticMaterial_h


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <Vector.h>
#include <Matrix.h>
#include <Mate.h>
#include <Constants.h>
#include <ElasticIsotropic.h>


namespace OpenSees {
// template <int n>
// using Material = Mate<n>;

// matrix index to tensor index mapping
class TensorIndexing {
  static void map(int matrix_index, int& i, int& j);
};

template <int n, PlaneType type, typename index>
class PlasticMaterial : public Mate<n> {
public:
  // null constructor
  PlasticMaterial();

  // full constructor
  PlasticMaterial(int tag,
               int classTag,
               double K,
               double G,
               double yield0,
               double yield_infty,
               double d,
               double H,
               double viscosity = 0,
               double rho       = 0.0);

  // elastic constructor
  PlasticMaterial(int tag, int classTag, double K, double G);

  // destructor
  virtual ~PlasticMaterial();

  virtual const char*
  getClassType() const
  {
    return "PlasticMaterial";
  }

  //swap history variables
  virtual int commitState();

  //revert to last saved state
  virtual int revertToLastCommit();

  //revert to start
  virtual int revertToStart();

  //get the strain and integrate plasticity equations
  virtual int setTrialStrain(const MatrixSD<n,true>& e) final;
  virtual int setTrialStrain(const MatrixND<n,n>& F) final;

  //unused trial strain functions
//int setTrialStrain(const MatrixSD<n,true>& v, const Vector& r);
//int setTrialStrainIncr(const MatrixSD<n,true>& v);
//int setTrialStrainIncr(const MatrixSD<n,true>& v, const Vector& r);


  MatrixSD<n,true> getStrain();
  const MatrixSD<n>& getStress();
  const Matrix& getTangent();
  const Matrix& getInitialTangent();

#if 0
  //sending and receiving
  virtual int sendSelf(int tag, Channel&);
  virtual int recvSelf(int tag, Channel&, FEM_ObjectBroker& theBroker);
#endif

  //print out material data
  void Print(OPS_Stream& s, int flag = 0);

  virtual Mate<n>* getCopy();
  virtual PlaneType getType() const;
  virtual int getOrder() const;

  double getRho()
  {
    return rho;
  }

#if 0
  virtual int setParameter(const char** argv, int argc, Parameter& param);
#endif
  virtual int updateParameter(int parameterID, Information& info);
  virtual int activateParameter(int paramID);

protected:
  // zero internal variables
  void zero();

  //plasticity integration routine
  void plastic_integrator();

  void doInitialTangent();

  // hardening function
  double q(double xi);

  //hardening function derivative
  double qprime(double xi);


  //material parameters
  double bulk;        //bulk modulus
  double shear;       //shear modulus
  double sigma_0;     //initial yield stress
  double sigma_infty; //final saturation yield stress
  double delta;       //exponential hardening parameter
  double Hard;        //linear hardening parameter
  double eta;         //viscosity

  // internal variables
  MatrixND<3,3> epsilon_p_n;      // plastic strain time n
  MatrixND<3,3> epsilon_p_nplus1; // plastic strain time n+1
  double xi_n;             // xi time n
  double xi_nplus1;        // xi time n+1

  // material response
  MatrixSD<3> stress;                       // stress tensor
  double tangent[3][3][3][3];               // material tangent
  double initialTangent[3][3][3][3];        // initial tangent

  MatrixSD<6> tangent_matrix;

  // material input
  MatrixSD<n> strain; //strain tensor

  static constexpr double root23
      = OpenSees::Constants::sqrt2/OpenSees::Constants::sqrt3;


  double rho;

  int parameterID;

}; // class PlasticMaterial
}

#include "PlasticMaterial.tpp"
// namespace OpenSees

#endif

