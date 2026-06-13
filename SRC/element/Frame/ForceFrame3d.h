//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
#ifndef ForceFrame3d_h
#define ForceFrame3d_h
//
#include <array>
#include <vector>
#include <element/Frame/BasicFrame3d.h>
#include <Vector.h>
#include <VectorND.h>
#include <MatrixND.h>
#include <FrameSection.h>
#include <BasicFrameTransf.h>
#define BASIC_TRANSFORM 1

class Matrix;
class Channel;
class Response;
class ElementalLoad;
class BeamIntegration;
namespace OpenSees {
class FrameTransformBuilder;
}

using namespace OpenSees;

template <int NIP, int nsr, int nwm, int shear_flag>
class ForceFrame3d: public BasicFrame3d, 
                    public FiniteElement<2, 3, 6+nwm>
{
 public:
  ForceFrame3d(int tag,
               std::array<int,2>& nodes,
               std::vector<FrameSection*>& sections,
               BeamIntegration &,
               FrameTransformBuilder &, 
               double density, int mass_flag, bool use_density,
               int max_iter, double tolerance
  );

  ~ForceFrame3d();

private:
    constexpr static int 
        NDF = 6+nwm,
        ndm = 3,        // dimension of the problem (3D)
        NEN = 2,        // number of element nodes
        NBV = 6+nwm*2;  // number of element DOFs in the basic system

public:

  const char *
  getClassType() const final {
    return "ForceFrame";
  }

  int setNodes();
  int commitState() override;
  int revertToLastCommit() override;        
  int revertToStart() override;
  int update() override;

  virtual const Matrix &getMass() final;
  virtual const Matrix &getTangentStiff() final;
  const Matrix &getInitialStiff() final;
  const Vector &getResistingForce() final;

  void zeroLoad() override {
    this->BasicFrame3d::zeroLoad();
    this->FiniteElement<2, 3, NDF>::zeroLoad();
  }

  int addLoad(ElementalLoad *theLoad, double loadFactor) final {
    return this->BasicFrame3d::addLoad(theLoad, loadFactor);
  }

  // bool hasMass() const final {
  //   return (mass_flag != 0) || use_density;
  // }

  /*
  const Vector &getResistingForceIncInertia();
  int addInertiaLoadToUnbalance(const Vector &accel); 
  */
  
  Response *setResponse(const char **argv, int argc, OPS_Stream &) final;
  int getResponse(int responseID, Information &) final;
  
  // Element: Parameters
  int setParameter(const char **argv, int argc, Parameter &) final;
  int updateParameter(int parameterID, Information &) final;
  // int activateParameter(int parameterID);

  // Element: Sensitivity
  const Vector &getResistingForceSensitivity(int gradNumber) final;
  int commitSensitivity(int gradNumber, int numGrads) final;
  int getResponseSensitivity(int responseID, int gradNumber, Information &) final;


  // MovableObject
  int sendSelf(int cTag, Channel &) override;
  int recvSelf(int cTag, Channel &, FEM_ObjectBroker &) override;
  
  // TaggedObject
  void Print(OPS_Stream &, int flag) override;    
  

 private:
  //
  // Constexpr
  //

  constexpr static int NNW = 6; // number of non-warping basic DOFs

  // static constexpr int shear_flag = (nsr-2*nwm == 6) ? 1 : 0;

  static constexpr FrameStressLayout scheme = {
    FrameStress::N,
    FrameStress::T,
    FrameStress::My,
    FrameStress::Mz,
    shear_flag ? FrameStress::Vy : FrameStress::Bimoment,
    shear_flag ? FrameStress::Vz : FrameStress::Bishear,
    FrameStress::Bimoment,
    FrameStress::Bishear
  };
  
  enum : int {
#if BASIC_TRANSFORM == 0
    inx = -12, //  0
    iny = -12, //  1
    inz = -12, //  2
    imx = -12, //  3
#endif
    imy =   3, //  4
    imz =   1, //  5
    iwx =   6, //
    //
    jnx =   0, //  6
#if BASIC_TRANSFORM == 0
    jny = -12, //  7
    jnz = -12, //  8
#endif
    jmx =   5, //  9
    jmy =   4, // 10
    jmz =   2, // 11
    jwx =   7,
  };

  enum Respond: int {
    GlobalForce = 1,
    BasicPlasticDeformation = 4,
    LocalForce  = 2,
    BasicForce  = 7,
    BasicStiff  =19,
    ResultantGradient=76,
    ResultantStiffness=1001
  };

  struct SectionState {
      VectorND<nsr>       es;      // section deformations
      MatrixND<nsr,nsr>   Fs;      // section flexibility
  };
  //
  // Functions
  //
  int update01();
  // int update02();
  // int updateMixed02();
  // int solveMixed02(const VectorND<NBV> &v_trial,
  //                 VectorND<NBV> &q_pres,
  //                 std::array<SectionState, NIP> &trial,
  //                 double L) noexcept;
  int getInitialFlexibility(MatrixND<NBV,NBV> &Fe);
  int getInitialDeformations(Vector &v0);

  void addLoadAtSection(VectorND<nsr> &sp, double x);
  void addLoadTangent(MatrixND<2*NDF,2*NDF>& K, double c);

  int setSectionPointers(std::vector<FrameSection*>&);
  void initializeSectionHistoryVariables();
  int getIntegral(Field field, State state, double& total);

  // Sensitivity
  int parameterID;
  VectorND<6+nwm*2> getBasicForceGrad(int gradNumber);
  void getStressGrad(VectorND<nsr> &dspdh, int isec, int gradNumber);

  //
  // Data
  //

  //
  // Element State
  //
  // Parameters
  double density;                // mass density per unit length
  double twist_mass;
  double total_mass;
  int    mass_flag;
  bool   mass_initialized;
  bool   use_density;

  int    max_iter;               // maximum number of local iterations
  double tol;	                   // tolerance for relative energy norm for local iterations


  MatrixND<NBV,NBV> K_pres,      // stiffness matrix in the basic system 
                    K_save;      // committed stiffness matrix in the basic system
  VectorND<NBV> q_pres,          // element resisting forces in the basic system
                q_save,          // committed element end forces in the basic system
                v_past;
  
  int    state_flag;             // indicate if the element has been initialized


  //
  // Section State
  //
  struct GaussPoint {
    double point,
           weight;
    FrameSection* material;
    MatrixND<nsr,nsr> Fs;         // Section flexibility
    VectorND<nsr>     es;         // Section deformations
    VectorND<nsr>     sr;         // Section stress resultants
    VectorND<nsr> es_save;        // Committed section deformations
  };

  std::vector<GaussPoint> points;
  BeamIntegration*        stencil;

  BasicFrameTransf3d<NDF> *basic_system;
  Matrix *Ki;
};

#define THREAD_LOCAL  static
#define ALWAYS_STATIC static // used when we need to do things like return a Matrix reference
#include "ForceFrame3d.tpp"
#undef THREAD_LOCAL
#undef ALWAYS_STATIC
#endif
