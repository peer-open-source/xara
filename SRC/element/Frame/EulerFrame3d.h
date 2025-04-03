//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
//
#ifndef EulerFrame3d_h
#define EulerFrame3d_h

#include <array>
#include <Matrix.h>
#include <Vector.h>
#include <Frame/BasicFrame3d.h>
#include <FrameSection.h>

class Node;
class FrameTransform3d;
class BeamIntegration;
class Response;

class EulerFrame3d : public BasicFrame3d
{
  public:
    EulerFrame3d(int tag, std::array<int,2>& nodes,
                 int numSections, FrameSection **s,
                 BeamIntegration &bi, 
                 FrameTransform3d &coordTransf,
                 double rho, int mass_flag
    );
    EulerFrame3d();
    ~EulerFrame3d();

    const char *getClassType() const {return "EulerFrame3d";}
    static constexpr const char* class_name = "EulerFrame3d";


    // public methods to set the state of the element   
    int commitState();
    int revertToLastCommit();
    int revertToStart();

    // public methods to obtain stiffness, mass, damping and residual information    
    int update();
    virtual const Vector &getResistingForce() final;
    virtual const Matrix &getMass() final;
//  const Matrix &getTangentStiff();
//  const Matrix &getInitialStiff();

//  void zeroLoad();
//  int addLoad(ElementalLoad *theLoad, double loadFactor);
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
    // For BasicFrame3d
    virtual  OpenSees::VectorND<6>&   getBasicForce() final;
    virtual  OpenSees::MatrixND<6,6>& getBasicTangent(State state, int rate) final;
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
        maxNumSections = 20;
  constexpr static int max_nip = 20;

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
};

#endif

