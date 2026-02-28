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
// Description: FrameSolidSection3d provides the abstraction of a 
// 3D beam section discretized by fibers. The section stiffness and
// stress resultants are obtained by summing fiber contributions.
//
// Written: cmp
// Created: Spring 2025
//
#ifndef FrameSolidSection3d_h
#define FrameSolidSection3d_h
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
class Response;

namespace OpenSees {
class FrameSolidSection3d : public FrameSection
{
  public:
    FrameSolidSection3d(); 
    FrameSolidSection3d(int tag, int numFibers);
  private:
    FrameSolidSection3d(const FrameSolidSection3d &);
  public:
    ~FrameSolidSection3d();

    int addFiber(NDMaterial&, double area, double y, double z=0.0);

    // Element
    const char *getClassType() const {
      return "FrameSolidSection3d";
    }

    int   setTrialSectionDeformation(const Vector &e); 
    const Vector &getSectionDeformation();

    int   getIntegral(Field field, State state, double& value) const final;
    const Vector &getStressResultant() override;
    const Matrix &getSectionTangent() override;
    const Matrix &getInitialTangent() override;

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
        alpha,
        ShearAlignYY,
        ShearAlignZZ,
        ShearAlignYZ,
        ShearAlignZY,
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

    struct Tangent {
      OpenSees::MatrixND<3,3> nn, nm, nw, nv,
                              mn, mm, mw, mv, 
                              wn, wm, ww, wv,
                              vn, vm, vw, vv;
      void zero() {
        nn.zero(); nm.zero(); nw.zero(); nv.zero();
        mn.zero(); mm.zero(); mw.zero(); mv.zero();
        wn.zero(); wm.zero(); ww.zero(); wv.zero();
        vn.zero(); vm.zero(); vw.zero(); vv.zero();
      }
    };

    Tangent K_pres;
    std::shared_ptr<Tangent> K_init;
  
    int stateDetermination(Tangent& K, VectorND<nsr>* s_trial, const VectorND<nsr> * const e_trial, int tangentFlag);

    
    struct FiberData {
      // NOTE: this may be faster if qualified with const, but this would
      // prevent parameterizing the fiber data.
      double y;
      double z;
      double area;
      std::array<std::array<double,3>,nwm> warp{{{0}}};
      // std::array<double, 3> wmix{{{0}}};
      OpenSees::VectorND<3> r;
    };
    const std::shared_ptr<std::vector<FiberData>> fibers;
    std::vector<NDMaterial*> materials;

    VectorND<nsr> s, e;
    Vector s_wrap, e_wrap;
    Matrix3D shear_align;
    Vector3D centroid;
    double nubar;
    enum class FiberState {
      Dirty, Clean
    } fiber_state;

    bool wagner;
    int parameterID;

    static ID code;

    Vector dedh;
};

#endif
} // namespace OpenSees