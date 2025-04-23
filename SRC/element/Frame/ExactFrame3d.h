//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Claudio M. Perez
//
#pragma once
#include <set>
#include <array>
#include <Flag.h>
#include <Element.h>
#include <VectorND.h>
#include <MatrixND.h>

#include <Frame/FrameMass.hpp>
#include <FrameSection.h>
#include <FrameTransform.h>
#include <Frame/FiniteElement.h>
#include <Rotations.hpp>

namespace OpenSees {

template<std::size_t nen, int nwm=0>
class ExactFrame3d: 
  public FiniteElement<nen, 3, 6+nwm>
{
public:
  enum Logarithm {
    None,
    LogLeft,
    LogRight
  };

  ExactFrame3d(int tag, std::array<int,nen>& nodes,
               FrameSection *section[nen-1], 
               FrameTransform3d& transf
//             int mass_type
  );

  ~ExactFrame3d();

  const char*
  getClassType() const final
  {
    return "ExactFrame3d";
  }

  virtual int setNodes();
  virtual const Vector &getResistingForce();
  virtual const Matrix &getTangentStiff();
  virtual const Matrix &getMass();
  virtual const Matrix &getInitialStiff();

  virtual int update() final;
  virtual int addLoad(ElementalLoad* , double scale) override final;

  virtual int revertToStart();
  virtual int revertToLastCommit();
  virtual int commitState() final {

    if (this->Element::commitState() != 0)
      opserr << "ExactFrame3d::commitState () - failed in base class";

    past = pres;

    for (GaussPoint& point : pres) {
      if (point.material->commitState() != 0)
        return -1;
    }
    return OpenSees::Flag::Success;
  }

  // Element: Parameters
  virtual int setParameter(const char **argv, int argc, Parameter &);
  virtual int updateParameter(int parameterID, Information &);
  virtual int activateParameter(int parameterID);

  Response *setResponse(const char **argv, int argc, OPS_Stream &s) final;
  virtual int getResponse(int responseID, Information &) final;

  // MovableObject
  int sendSelf(int cTag, Channel&);
  int recvSelf(int cTag, Channel&, FEM_ObjectBroker&);

  // TaggedObject
  void Print(OPS_Stream& s, int flag = 0);

  private:
    //
    // Constexpr
    //
    constexpr static int 
          nsr = 6+2*nwm,        // Number of section resultants
          ndm = 3,              // Dimension of the problem (3D)
          nip = nen-1,          // Number of integration points
          ndf = 6+nwm;          // Degrees of freedom per node

    enum Respond: int {
      GlobalForce = 1,
      BasicPlasticDeformation = 4,
      LocalForce  = 2,
      BasicForce  = 7,
      BasicStiff  =19,
    };
    
    // Layout of stress resultants
    static constexpr FrameStressLayout scheme = {
      FrameStress::N,
      FrameStress::Vy,
      FrameStress::Vz,
      FrameStress::T,
      FrameStress::My,
      FrameStress::Mz,
      FrameStress::Bimoment,
  //  FrameStress::By,
  //  FrameStress::Bz,
      FrameStress::Bishear,
  //  FrameStress::Qy,
  //  FrameStress::Qz
    };

    //
    //
    //
    struct GaussPoint {
      double point,
             weight;
      FrameSection* material;

      double   shape[2][nen];
      Matrix3D rotation;
      Vector3D curvature;
    };
    // Node locations in local (scalar) coordinate
    double xn[nen];
    double jxs;

    Matrix3D R0;

    std::set<FrameLoad*> frame_loads;

    std::array<GaussPoint,nip> pres;
    std::array<GaussPoint,nip> past;
    FrameTransform3d*       transform;
    Logarithm               logarithm;
    BeamIntegration*        stencil;

    //
    //
    VectorND<nen*ndf> p;
    MatrixND<nen*ndf,nen*ndf> K;
    // Parameters
    int parameterID;
};

} // namespace OpenSees
#include "ExactFrame3d.tpp"
