//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
#ifndef CubicFrame3d_h
#define CubicFrame3d_h

#include <array>
#include <vector>
#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <VectorND.h>
#include <FrameSection.h>

#include <ID.h>

class Node;
class CrdTransf;
class BeamIntegration;
class Response;
class BasicFrame3d;

namespace OpenSees {

template <bool shear, int nwm>
class CubicFrame3d : public Element {
public:
  CubicFrame3d(int tag, 
               std::array<int, 2>&,
               std::vector<FrameSection*>&, 
               BeamIntegration&,
               CrdTransf&, 
               double rho);
  CubicFrame3d();
  ~CubicFrame3d();


  const char*
  getClassType() const final
  {
    return "CubicFrame3d";
  }

  int getNumExternalNodes() const;
  const ID& getExternalNodes();
  Node** getNodePtrs();

  int getNumDOF();
  void setDomain(Domain*);

  // Element: State
  int commitState();
  int revertToLastCommit();
  int revertToStart();

  // public methods to obtain stiffness, mass, damping and residual information
  int update();
  const Matrix& getTangentStiff();
  const Matrix& getInitialStiff();
  const Matrix& getMass();

  void zeroLoad();
  int addLoad(ElementalLoad* theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector& accel);

  const Vector& getResistingForce();
  const Vector& getResistingForceIncInertia();


  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information& eleInfo);

  // Element: Parameters
  int setParameter(const char** argv, int argc, Parameter& param);
  int updateParameter(int parameterID, Information& info);
  int activateParameter(int parameterID);

  // Element: Sensitivity
  const Vector& getResistingForceSensitivity(int gradNumber);
  const Matrix& getInitialStiffSensitivity(int gradNumber);
  const Matrix& getMassSensitivity(int gradNumber);
  int commitSensitivity(int gradNumber, int numGrads);

  // MovableObject
  int sendSelf(int commitTag, Channel&);
  int recvSelf(int commitTag, Channel&, FEM_ObjectBroker&);

  // TaggedObject
  void Print(OPS_Stream& s, int flag = 0);

protected:

private:
  void getBasicStiff(Matrix& kb, int initial = 0);
private:
  constexpr static int 
        nsr = shear? 6 : 4,   // number of section resultants
        ndm = 3,              // dimension of the problem (3D)
        nq  = 6+nwm*2,        // number of element dof's in the basic system
        NDF = 6+nwm,          //
        NEN = 2,              // number of element nodes
        maxNumSections = 20;


  int numSections;
  FrameSection** theSections;       // the materials
  CrdTransf* theCoordTransf; // coordinate transformation object
  BeamIntegration* beamInt;

  double xi[maxNumSections];
  double wt[maxNumSections];
  double phizs[maxNumSections]; // Shear term 12EIz/(GA L^2)
  double phiys[maxNumSections]; // Shear term 12EIy/(GA L^2)

  ID connectedExternalNodes;    // Tags of nodes

  Node* theNodes[NEN];


  Vector Q;       // Applied nodal loads
  VectorND<nq> q;

  double density;             // Mass density per unit length
  double total_mass,
         twist_mass;
  int    mass_flag;
  bool   use_density;
  bool   mass_initialized;

  BasicFrame3d* loads;

  int parameterID;

  // Matrix K; // Element stiffness, damping, and mass Matrix
  MatrixND<NEN*NDF,NEN*NDF> K;
  VectorND<NDF*NEN> P; // Element resisting force vector
  Matrix K_wrap;
  Vector P_wrap;

  static constexpr FrameStressLayout scheme = {
      FrameStress::N, 
      FrameStress::T, 
      FrameStress::My, 
      FrameStress::Mz,
      FrameStress::Vy, 
      FrameStress::Vz,
  };
};
}
#include <CubicFrame3d.tpp>
#endif
