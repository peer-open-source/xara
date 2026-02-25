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
// Description: MixedFrameSection provides the abstraction of a 
// 3D beam section discretized by fibers. The section stiffness and
// stress resultants are obtained by summing fiber contributions.
//
// Written: cmp
// Created: Jan. 2026
//
#ifndef MixedFrameSection_h
#define MixedFrameSection_h
#include <array>
#include <memory>
#include <FrameSection.h>
#include <Vector.h>
#include <vector>
#include <Matrix.h>
#include <MatrixND.h>
#include <VectorND.h>
#include <Matrix3D.h>
#include "FrameSectionConstants.h"

class NDMaterial;
class MaterialBuilder;
class Response;

namespace OpenSees {
class MixedFrameSection : public FrameSection
{
  public:
    using MixedType = Frame::Prism::MixedType;
    MixedFrameSection(int tag, int reserve, MixedType type);
  private:
    MixedFrameSection(const MixedFrameSection &);
  public:
    ~MixedFrameSection();

    int addFiber(MaterialBuilder&, double area, double y, double z=0.0);

    //
    const char *getClassType() const {
      return "MixedFrameSection";
    }

    int   setTrialSectionDeformation(const Vector &e);
    const Vector &getSectionDeformation();

    int   getIntegral(Field field, State state, double& value) const final;
    const Vector &getStressResultant() override;
    const Matrix &getSectionTangent() override;
    const Matrix &getInitialTangent() override;
    MatrixND<12,12> getFullTangent(State state) override;

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
    constexpr static int nwm =  3; // Number of warping shapes
    constexpr static int nem =  0; // number of enhanced modes
    constexpr static int nep =  3;

    VectorND<nsr> s, e;
    Vector s_wrap, e_wrap;
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

    int formMixedUniformL(Matrix3D& Lr, Matrix3D& Lw) const;

    int solveMixed(const VectorND<nsr>& e, MatrixND<6,6>& Kee, Tangent& Ks);
  
    int stateDetermination(Tangent& K, VectorND<nsr>* s_trial, const VectorND<nsr> * const e_trial, int tangentFlag);

    int checkFiberState();

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


    // Mixed formulation data
    MixedType mixed_type;
    enum MixedShapes: int {
      TwistX = 1<< 0,
      ShearY = 1<< 1,
      ShearZ = 1<< 2,
      TwistE = 1<< 3,
    };
    int mixed_shapes = 0;
    Matrix3D    shear_align;
    VectorND<3> shift_twist, shift_axial;

    // Centroid
    Vector3D centroid;
    // Average material Poisson's ratio
    double nubar;
    enum class FiberState {Dirty, Clean} fiber_state;
    State tangent_state = State::None;

    const bool wagner;
    int parameterID;

    static ID code;

    Vector dedh;


    static constexpr int MaxThreads = 12;
    int num_threads = 8;
    void *pool;        // thread pool

    inline int 
    RigidShape(const FiberData& fiber, double aw, MatrixND<3,6>& Ae) const noexcept {
      const VectorND<2>& r = fiber.r;
      Ae(0,0) = 1.0;
      Ae(1,1) = 1.0;
      Ae(2,2) = 1.0;
      Ae(0,4) =  r[1];
      Ae(0,5) = -r[0];
      Ae(1,3) = -r[1];
      Ae(2,3) =  r[0];
      return 0;
    }

    inline int
    WarpShape(const FiberData& fiber, Matrix3D& iow, Matrix3D& iodw) const noexcept {

      const FiberData::WarpArray& w = fiber.warp;
      switch (mixed_type) {
        case MixedType::UT:
        case MixedType::Energetic:
        case MixedType::Constant:
          return 0;
        case MixedType::None:
          iow(0,1)  = w[1][0];
          iow(0,2)  = w[2][0];
          iodw(1,1) = w[1][1]; // v_y,y
          iodw(2,1) = w[1][2];
          iodw(1,2) = w[2][1];
          iodw(2,2) = w[2][2];
          // fallthrough
        case MixedType::Equilibrium:
          iow(0,0)  = w[0][0];
          iodw(1,0) = w[0][1];
          iodw(2,0) = w[0][2];
          return 1;
      }
    }

    inline void 
    MixedShape(FiberData& fiber, const Matrix3D& Gr, const Matrix3D& Gw, Matrix3D& An) const {
      
      constexpr static Matrix3D oneS {{
        0.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0
      }};

      const FiberData::WarpArray& w = fiber.warp;
      if (mixed_type == MixedType::None)
        return;
      else if (mixed_type == MixedType::UT) {
        An(1,2) = w[0][1];
        An(2,2) = w[0][2];
      }
      else if (mixed_type == MixedType::Equilibrium) {
        double b = Gw(2,2);
        An(1,2) = w[0][1] + b*w[1][1];
        An(2,2) = w[0][2] + b*w[1][2];
      }
      else {
        const Vector3D r = {0.0, fiber.r[0], fiber.r[1]};
        const Matrix3D Anr {{
          0.0,   1.0,  0.0, // Vy
          0.0,   0.0,  1.0, // Vz
          0.0, -r[2], r[1], // T
        }};

        An = Anr*Gr;
        {
          Matrix3D Anwo {{
            0.0, w[0][1], w[0][2],
            0.0, w[1][1], w[1][2],
            0.0, w[2][1], w[2][2]
          }};
          // -(dev ror)Xi
          Anwo.addTensorProduct(r, r, -nubar);
          Anwo.addMatrix(oneS, 0.5*r.dot(r)*nubar);
          // + roc 
          Anwo.addTensorProduct(r, centroid, nubar);
          An.addMatrixProduct(Anwo, Gw, 1.0);
        }
      }
    }
};

#endif
} // namespace OpenSees