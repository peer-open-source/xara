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
// Description: This file contains the class definition for 
// FrameFiberSection3d.h. FrameFiberSection3d provides the abstraction of a 
// 3d beam section discretized by fibers. The section stiffness and
// stress resultants are obtained by summing fiber contributions.
//
// Written: cmp
// Created: 2024
//
#ifndef FrameFiberSection3d_h
#define FrameFiberSection3d_h

#define SEC_TAG_FrameFiberSection3d 0

#include <FrameSection.h>
#include <Vector.h>
#include <Matrix.h>
#include <VectorND.h>
#include <MatrixND.h>
#include <vector>
#include <memory>
#include <Frame/Shape.h>

class Response;
class UniaxialMaterial;

namespace OpenSees {

class FrameFiberSection3d : public FrameSection
{
  public:
    explicit
    FrameFiberSection3d(int tag,
                        int numFibers,
                        const Frame::Shape& shape,
                        double mass, 
                        bool use_mass
                      );
  private:
    FrameFiberSection3d(const FrameFiberSection3d &); 
  public:
    ~FrameFiberSection3d();

    const char *getClassType() const {
      return "UniaxialFiberSection3d";
    }

    int   setTrialSectionDeformation(const Vector &e);
    const Vector &getSectionDeformation();

    int   getIntegral(Field field, State state, double& value) const override final;
    const Vector &getStressResultant();
    const Matrix &getSectionTangent();
    const Matrix &getInitialTangent();

    VectorND<12> getFullStress() noexcept final {
      return sr;
    }

    MatrixND<12,12>
    getFullTangent(State state) noexcept final {
      if (state == State::Pres)
        return tangent.matrix;
      else {
        Tangent T{};
        VectorND<nsr> s;
        this->initializeShear(T, *shape);
        this->updateAxial(state, T, s);
        return T.matrix;
      }
    }

    int   commitState();
    int   revertToLastCommit();    
    int   revertToStart();
 
    FrameSection *getFrameCopy();
    const ID &getType();
    int getOrder () const; // {return 4;};
 
    int sendSelf(int cTag, Channel &) override;
    int recvSelf(int cTag, Channel &, FEM_ObjectBroker &) override;
    void Print(OPS_Stream &s, int flag) override;
	    
    Response *setResponse(const char **argv, int argc, OPS_Stream &) override;
    int getResponse(int responseID, Information &) override;

    int addFiber(UniaxialMaterial &theMat, double area, double y, double z=0.0);

    int setParameter(const char **argv, int argc, Parameter &);
    int updateParameter(int id, Information &) override;
    const Vector & getStressResultantSensitivity(int gradIndex, bool conditional);
    const Matrix & getSectionTangentSensitivity(int gradIndex);
    int   commitSensitivity(const Vector& sectionDeformationGradient, int gradIndex, int numGrads);

    const Vector & getSectionDeformationSensitivity(int gradIndex);

    double getEnergy() const;


  protected:
    constexpr static int nsr = 12;
    constexpr static int nwm = 3;

  private:
    struct FiberData {
      double y;
      double z;
      double area;
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

    struct Tangent {
      MatrixND<nsr,nsr> matrix;
      void zeroAxial() noexcept {
        matrix(inx, inx) = 0.0;
        matrix(imz, imz) = 0.0;
        matrix(imy, imy) = 0.0;
        matrix(imz, imy) = 0.0;
        matrix(imy, imz) = 0.0;
        matrix(inx, imz) = 0.0;
        matrix(imz, inx) = 0.0;
        matrix(inx, imy) = 0.0;
        matrix(imy, inx) = 0.0;
      }
      void zeroShear() noexcept {
        matrix(iny, iny) = 0.0;
        matrix(inz, inz) = 0.0;
        matrix(imx, imx) = 0.0;
      }
    } tangent;
    
    std::shared_ptr<Frame::Shape> shape;
    std::shared_ptr<std::vector<FiberData>> fibers;
    std::vector<UniaxialMaterial*> materials;

    Matrix K_wrap;


    static ID code;

    VectorND<nsr> es, sr;
    Vector  e;         // trial section deformations 
    Vector  s;         // section resisting forces  (axial force, bending moment)

    int parameterID =0;

    int updateAxial(const State, Tangent&, VectorND<nsr>& s) const;
    static int initializeShear(Tangent&, Frame::Shape&);

    inline int 
    FiberGrad(int i, int j, double& dA, double& dy, double& dz) const noexcept
    {
      dA = 0.0; dy = 0.0; dz = 0.0;
      if (parameterID >= Param::FiberFieldBase) {
        int fiberID = (parameterID - Param::FiberFieldBase) / 100;
        int field   = (parameterID - Param::FiberFieldBase) % 100;
        if (i == fiberID) {
          switch (field) {
          case Param::FiberArea:   dA = 1.0; break;
          case Param::FiberY:      dy = 1.0; break;
          case Param::FiberZ:      dz = 1.0; break;
          default:
            break;
          }
        }
      }
      return 0;
    }
};

} // namespace OpenSees
#endif
