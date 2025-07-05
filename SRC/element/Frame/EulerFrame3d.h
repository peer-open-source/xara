//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//----------------------------------------------------------------------------//
//
//
#ifndef EulerFrame3d_h
#define EulerFrame3d_h

#include <array>
#include <Matrix.h>
#include <Vector.h>
#include <Frame/BasicFrame3d.h>
#include <FrameSection.h>
#include <BasicFrameTransf.h>

class Node;
class BeamIntegration;
class Response;

namespace OpenSees {
  class FrameTransformBuilder;

class EulerFrame3d : public BasicFrame3d,
                     public FiniteElement<2, 3, 6>
{
  public:
    EulerFrame3d(int tag, 
                 std::array<int,2>& nodes,
                 int numSections,
                 FrameSection **s,
                 BeamIntegration &bi, 
                 FrameTransformBuilder& tb,
                 double rho, int mass_flag
    );
    ~EulerFrame3d();

    const char *getClassType() const {return "EulerFrame3d";}


    // public methods to set the state of the element   
    int commitState();
    int revertToLastCommit();
    int revertToStart();

    // public methods to obtain stiffness, mass, damping and residual information    
    int update() final;
    const Vector &getResistingForce() final;
    const Matrix &getMass() final;
    const Matrix &getTangentStiff() final;
    const Matrix &getInitialStiff() final;

    void zeroLoad() final {
      this->BasicFrame3d::zeroLoad();
      this->FiniteElement<2, 3, 6>::zeroLoad();
    }
    
    int addLoad(ElementalLoad *theLoad, double loadFactor) final {
      return this->BasicFrame3d::addLoad(theLoad, loadFactor);
    }
//  int addInertiaLoadToUnbalance(const Vector &accel);
//  const Vector &getResistingForceIncInertia();

    // public methods for element output
    int sendSelf(int commitTag, Channel &);
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker&);

    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &);


    virtual int setParameter(const char **argv, int argc, Parameter &param) final;
//  virtual int updateParameter(int parameterID, Information &);
//  virtual int activateParameter(int parameterID) final;

    // Sensitivity
    virtual const Vector & getResistingForceSensitivity(int gradNumber);
    virtual const Matrix & getInitialStiffSensitivity(int gradNumber);
//  virtual const Matrix & getMassSensitivity(int gradNumber);
    virtual int            commitSensitivity(int gradNumber, int numGrads);


    virtual int getIntegral(Field field, State state, double& total);


protected:
    // For FiniteElement
    virtual int setNodes() final;

private:
  //
  // Constexpr
  //
  constexpr static int 
        nsr = 4,              // number of section resultants
        nen = 2,              // number of element nodes
        ndm = 3,              // dimension of the problem (3D)
        nq  = 6,              // number of element dof's in the basic system; N,My,Mz
        NDF = 6,
        maxNumSections = 20;
  constexpr static int max_nip = 20;

  OpenSees::VectorND<6>&   getBasicForce();
  OpenSees::MatrixND<6,6>& stateDetermination(State state, int rate);

  VectorND<nq> getBasicForceGrad(int index);

//  const Matrix &getInitialBasicStiff();

    int numSections;

    //
    // Section State
    //
    struct GaussPoint {
      double point,
             weight;
      FrameSection* material;
    };
    double xi[max_nip];
    double wt[max_nip];

    std::vector<GaussPoint> points;
    BeamIntegration*        stencil;

    OpenSees::MatrixND<6,6> kb;
    OpenSees::VectorND<6>   q;

    double density;             // Mass density per unit length
    double total_mass,
           twist_mass;
    int    mass_flag;
    bool   use_density;
    bool   mass_initialized;

    static constexpr FrameStressLayout scheme = {
      FrameStress::N,
      FrameStress::Mz,
      FrameStress::My,
      FrameStress::T,
    };
    BasicFrameTransf3d<NDF> *basic_system;
};
} // namespace OpenSees
#endif

