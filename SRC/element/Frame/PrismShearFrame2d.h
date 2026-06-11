//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// PrismShearFrame2d is a shear-deformable frame member that takes
// rotational inertia effects into account.
//
// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 03/13
// Revision: A
//
#ifndef PrismShearFrame2d_h
#define PrismShearFrame2d_h
//
#include <Element.h>
#include <Matrix.h>
#include <MatrixND.h>
#include <VectorND.h>

class CrdTransf;

class PrismShearFrame2d : public Element {
public:
  PrismShearFrame2d(int tag, int Nd1, int Nd2, double E, double G, double A, double Iz, double Avy,
                    CrdTransf& theTransf, double rho = 0.0, int cMass = 0);
  PrismShearFrame2d();

  ~PrismShearFrame2d();

  const char*
  getClassType() const
  {
    return "PrismShearFrame2d";
  }

  // public methods to obtain information about dof & connectivity
  int getNumExternalNodes() const;
  const ID& getExternalNodes();
  Node** getNodePtrs();
  int getNumDOF();
  void setDomain(Domain* theDomain);

  // public methods to set the state of the element
  int commitState();
  int revertToLastCommit();
  int revertToStart();
  int update();

  // public methods to obtain stiffness, mass, damping and residual information
  const Matrix& getTangentStiff();
  const Matrix& getInitialStiff();
  const Matrix& getMass();

  void zeroLoad();
  int addLoad(ElementalLoad* theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector& accel);

  const Vector& getResistingForce();
  const Vector& getResistingForceIncInertia();

  // public methods for element output
  int sendSelf(int commitTag, Channel& theChannel);
  int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

  void Print(OPS_Stream& s, int flag = 0);

  // public methods for element recorder
  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information& info);

  int setParameter(const char** argv, int argc, Parameter& param);
  int updateParameter(int parameterID, Information& info);

private:
  // private methods
  void setUp();

  // private attributes - a copy for each object of the class
  ID connectedExternalNodes; // contains the tags of the end nodes
  Node* theNodes[2];         // array of nodes
  CrdTransf* theCoordTransf; // coordinate transformation

  // parameters
  double E;      // elastic modulus
  double G;      // shear modulus
  double A;      // cross section (axial) area
  double Iz;     // moment of inertia about local z axis
  double Avy;    // shear area along local y axis
  double phi;    // ratio of bending to shear stiffness
  double rho;    // mass per unit length
  int cMass;     // consistent mass flag
  int geom_flag; // nonlinear geometry flag


  // State variables
  double L;                  // element length
  Vector ul;                 // displacements in local system
  OpenSees::VectorND<6> ql;  // forces in local system
  OpenSees::VectorND<6> ql0; // fixed end forces due to loads in local system
  Matrix kl;                 // stiffness matrix in local system
  Matrix klgeo;              // geometric stiffness matrix in local system

  // Constant variables
  Matrix Tgl; // transformation matrix from global to local system
  Matrix Ki;  // initial stiffness matrix in global system
  Matrix M;   // mass matrix in global system

  static Matrix theMatrix; // a class wide Matrix
  static Vector theVector; // a class wide Vector
  Vector theLoad;
};

#endif
