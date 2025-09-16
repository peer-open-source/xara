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

class NDMaterial;
class Response;

namespace OpenSees {
class FrameSolidSection3d : public FrameSection
{
  public:
    FrameSolidSection3d(); 
    FrameSolidSection3d(int tag, int numFibers);
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
 
    FrameSection *getFrameCopy();
    const ID &getType();
    int getOrder () const;
    
    int sendSelf(int tag, Channel &) override;
    int recvSelf(int tag, Channel &, FEM_ObjectBroker &) override;
    void Print(OPS_Stream &s, int flag = 0) override;
	    
    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &info);


    // Sensitivity
    int setParameter(const char **argv, int argc, Parameter &);
    int updateParameter(int id, Information &);
    int activateParameter(int id);
    const Vector& getStressResultantSensitivity(int gradIndex, bool conditional);
    const Vector& getSectionDeformationSensitivity(int gradIndex);
    const Matrix& getInitialTangentSensitivity(int gradIndex);
    int commitSensitivity(const Vector& strainGrad, int gradIndex, int numGrads);


  private:
    constexpr static int nsr = 12;
    constexpr static int nwm =  3; // Number of warping modes

    FrameSection::Tangent K_pres;
    std::shared_ptr<Tangent> K_init;
  
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
        alpha,
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
      // NOTE: this may be faster if qualified with const, but this would
      // prevent parameterizing the fiber data.
      double y;
      double z;
      double area;
      std::array<std::array<double,3>,nwm> warp{{{0}}};
      // std::array<double, 3> wmix{{{0}}};
      OpenSees::VectorND<3> r;
    };
    std::shared_ptr<std::vector<FiberData>> fibers;
    std::vector<NDMaterial*> materials;

    VectorND<nsr> s, e;
    Vector s_wrap, e_wrap;

    Vector3D centroid;
    double yBar;                      // Section centroid
    double zBar;                      // Section centroid

    bool wagner;

    static ID code;

    int parameterID;
    Vector dedh;
};

#endif
} // namespace OpenSees