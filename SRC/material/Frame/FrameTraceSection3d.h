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
//===----------------------------------------------------------------------===//
//
// Description: FrameTraceSection3d provides the abstraction of a 
// 3D beam section discretized by fibers. The section stiffness and
// stress resultants are obtained by summing fiber contributions.
//
// Written: cmp
// Created: Jan. 2026
//
#ifndef FrameTraceSection3d_h
#define FrameTraceSection3d_h
#include <array>
#include <memory>
#include <FrameSection.h>
#include <Vector.h>
#include <vector>
#include <Matrix.h>
#include <MatrixND.h>
#include <VectorND.h>
#include <Matrix3D.h>

class NDMaterial;
class MaterialBuilder;
class Response;

namespace OpenSees {
class FrameTraceSection3d : public FrameSection
{
  public:
    FrameTraceSection3d(int tag, int reserve, bool wagner);
  private:
    FrameTraceSection3d(const FrameTraceSection3d &);
  public:
    ~FrameTraceSection3d();

    int addFiber(MaterialBuilder&, double area, double y, double z=0.0);

    //
    const char *getClassType() const {
      return "FrameTraceSection3d";
    }

    int   setTrialSectionDeformation(const Vector &e); 
    const Vector &getSectionDeformation();

    int   getIntegral(Field field, State state, double& value) const final;
    const Vector &getStressResultant() override;
    const Matrix &getSectionTangent() override;
    const Matrix &getInitialTangent() override;
    MatrixND<12,12> getFullTangent(State state) noexcept override;

    int   commitState() override;
    int   revertToLastCommit() override;    
    int   revertToStart() override;
 
    FrameSection *getFrameCopy() override;
    const ID &getType() override;
    int getOrder () const override;
    
    // MovableObject
    int sendSelf(int tag, Channel &) override;
    int recvSelf(int tag, Channel &, FEM_ObjectBroker &) override;
    // TaggedObject
    void Print(OPS_Stream &s, int flag) override;
	    
    Response *setResponse(const char **argv, int argc, OPS_Stream &) override;
    int getResponse(int responseID, Information &) override;


    // Sensitivity
    int setParameter(const char **argv, int argc, Parameter &) override;
    int updateParameter(int id, Information &) override;
    int activateParameter(int id) override;
    const Vector& getStressResultantSensitivity(int gradIndex, bool conditional);
    const Vector& getSectionDeformationSensitivity(int gradIndex);
    const Matrix& getInitialTangentSensitivity(int gradIndex);
    int commitSensitivity(const Vector& strainGrad, int gradIndex, int numGrads);


  private:
    constexpr static int nsr = 12;
    constexpr static int nwm =  3; // Number of warping modes

    struct Tangent {
      OpenSees::MatrixND<6,6> se;
      OpenSees::MatrixND<6,3> sw, sv;
      OpenSees::MatrixND<3,3> ww, vv, wv;
      void zero() {
        se.zero(); sw.zero(); sv.zero();
                   ww.zero(); wv.zero();
                              vv.zero();
      }
    };
    Tangent K_pres;
    std::shared_ptr<Tangent> K_init;

    int form_shifts(MatrixND<3,6>& Lw, MatrixND<6,6>& Lr) const;
  
    int stateDetermination(Tangent& K, VectorND<nsr>* s_trial, const VectorND<nsr> * const e_trial, int tangentFlag);

    struct Param {
      enum : int {
        FiberWarpX=0,
        FiberWarpXY,
        FiberWarpXZ,
        FiberWarpY,
        FiberWarpYY,
        FiberWarpYZ,
        FiberWarpZ,
        FiberWarpZY,
        FiberWarpZZ,
        FiberY,
        FiberZ,
        FiberArea,
        ShearAlignYY,
        ShearAlignZZ,
        ShearAlignYZ,
        ShearAlignZY,
        ShiftTwistY,
        ShiftTwistZ,
        ShiftAxialY,
        ShiftAxialZ,
        FiberFieldBase=10000
      };
    };

    enum: int {
      CurrentTangent, InitialTangent
    };

    enum : int {
      inx = 0, 
      iny,
      inz,

      imx,
      imy,
      imz,

      iwx,     
      iwy,
      iwz,

      ivx,     
      ivy,
      ivz
    };
    
    struct FiberData {
      double area;
      using WarpArray = std::array<std::array<double,3>,nwm>;
      WarpArray warp{{{0}}};
      OpenSees::VectorND<2> r;
    };
    const std::shared_ptr<std::vector<FiberData>> fibers;
    std::vector<NDMaterial*> materials;

    VectorND<nsr> s, e;
    Vector s_wrap, e_wrap;
    // Trace matrix blocks
    Matrix3D    shear_align;
    VectorND<3> shift_twist, shift_axial;
    // Centroid
    Vector3D centroid;
    // Average material Poisson's ratio
    double nubar;
    enum class FiberState {
      Dirty, Clean
    } fiber_state;

    const bool wagner;
    int parameterID;

    static ID code;

    Vector dedh;

    static constexpr int MaxThreads = 12;
    void *pool;        // thread pool
    bool thread_safe = true;
};

#endif
} // namespace OpenSees