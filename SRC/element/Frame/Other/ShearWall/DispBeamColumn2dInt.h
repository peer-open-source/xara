//
//
// Description: This file contains the class implementation of DispBeamColumn2dInt.
// Based on DispBeamColumn2d.cpp.
//
// 1. Massone, L. M., 2006; "RC Wall Shear - Flexure Interaction:
//   Analytical and Experimental Responses", Ph.D. Dissertation, University
//   of California, Los Angeles, June 2006, 398 pp.
// 2. Massone, L. M.; Orakcal, K.; and Wallace, J. W. , 2006; "Shear -
//   Flexure Interaction for Structural Walls"; SP-236, ACI Special
//   Publication - Deformation Capacity and Shear Strength of Reinforced
//   Concrete Members Under Cyclic Loading, editors: Adolfo Matamoros &amp;
//   Kenneth Elwood, p. 127-150.
//
// Created by: LMS
// Created: 07/04

#ifndef DispBeamColumn2dInt_h
#define DispBeamColumn2dInt_h

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <LegendreBeamIntegration.h>

class Node;
class SectionForceDeformation;
class FiberSection2dInt;
class CrdTransf;
class LinearCrdTransf2dInt;
class Response;

class DispBeamColumn2dInt : public Element

{

public:
  DispBeamColumn2dInt(int tag, int nd1, int nd2, int numSections,
                      SectionForceDeformation **s, CrdTransf &coordTransf, double C,
                      double rho = 0.0);

  DispBeamColumn2dInt();

  virtual ~DispBeamColumn2dInt();

  int getNumExternalNodes() const;
  const ID &getExternalNodes();
  Node **getNodePtrs();
  int getNumDOF();
  void setDomain(Domain *theDomain);

  // public methods to set the state of the element

  int commitState();
  int revertToLastCommit();
  int revertToStart();

  // public methods to obtain stiffness, mass, damping and residual information
  //
  int update();
  const Matrix &getTangentStiff();
  const Matrix &getInitialStiff();
  const Matrix &getMass();
  void zeroLoad();
  int addLoad(ElementalLoad *theLoad, double loadFactor);

  int addInertiaLoadToUnbalance(const Vector &accel);

  const Vector &getResistingForce();

  const Vector &getResistingForceIncInertia();

  // public methods for element output

  int sendSelf(int commitTag, Channel &);
  int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);


  void Print(OPS_Stream &s, int flag);

  Response *setResponse(const char **argv, int argc, OPS_Stream &s);

  int getResponse(int responseID, Information &eleInfo);

  // AddingSensitivity:BEGIN //////////////////////////////////////////
  const Vector &getResistingForceSensitivity(int gradNumber);
  const Matrix &getKiSensitivity(int gradNumber);
  const Matrix &getMassSensitivity(int gradNumber);
  int commitSensitivity(int gradNumber, int numGrads);
  // AddingSensitivity:END ///////////////////////////////////////////

private:
  const Matrix &getInitialBasicStiff();

  int numSections;

  FiberSection2dInt **theSections; // pointer to the ND material objects
  LinearCrdTransf2dInt *crdTransf; // pointer to coordinate transformation object

  double C1;
  ID connectedExternalNodes; // Tags of quad nodes

  Node *theNodes[2];

  static Matrix K;           // Element stiffness, damping, and mass Matrix
  static Vector P;           // Element resisting force vector

  Vector Q;                  // Applied nodal loads
  Vector q;                  // Basic force
  double q0[6];
  double p0[6];

  double rho; // Mass density per unit length

  static double workArea[];

  static LegendreBeamIntegration quadRule;

  // AddingSensitivity:BEGIN //////////////////////////////////////////
  int parameterID;
  // AddingSensitivity:END ///////////////////////////////////////////
};

#endif
