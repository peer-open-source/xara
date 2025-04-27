/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
//
// Description: This file contains the class definition for Truss. A Truss object
// provides the abstraction of the small deformation bar element. Each truss
// object is associated with a material object. This Truss element will work
// in 1d, 2d or 3d problems.
//
// Written: fmk
// Created: 07/98
// Revision: A
//
#ifndef Truss_h
#define Truss_h
//
#include <Element.h>
#include <Matrix.h>

class Node;
class Channel;
class UniaxialMaterial;

class Truss : public Element {
public:
  Truss(int tag, int dimension, int Nd1, int Nd2, UniaxialMaterial& theMaterial, double A,
        double rho = 0.0, int doRayleighDamping = 0, int cMass = 0, bool initDisp = true);

  Truss();
  ~Truss();

  const char*
  getClassType() const
  {
    return "Truss";
  }
  static constexpr const char* class_name = "Truss";

  // methods to obtain information about dof & connectivity
  int getNumExternalNodes() const;
  const ID& getExternalNodes();
  Node** getNodePtrs();

  int getNumDOF();
  void setDomain(Domain* theDomain);

  // methods to set the state of the element
  int commitState();
  int revertToLastCommit();
  int revertToStart();
  int update();

  // methods to obtain stiffness, mass, damping and residual information
  const Matrix& getKi();
  const Matrix& getTangentStiff();
  const Matrix& getInitialStiff();
  const Matrix& getDamp();
  const Matrix& getMass();

  void zeroLoad();
  int addLoad(ElementalLoad* theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector& accel);

  const Vector& getResistingForce();
  const Vector& getResistingForceIncInertia();


  // Element Sensitivity
  int setParameter(const char** argv, int argc, Parameter& param);
  int updateParameter(int parameterID, Information& info);
  int activateParameter(int parameterID);

  int addInertiaLoadSensitivityToUnbalance(const Vector& accel, bool tag);
  const Vector& getResistingForceSensitivity(int gradNumber);
  const Matrix& getKiSensitivity(int gradNumber);
  const Matrix& getMassSensitivity(int gradNumber);
  int commitSensitivity(int gradNumber, int numGrads);


  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information& eleInformation);

  // MovableObject
  int sendSelf(int commitTag, Channel& theChannel);
  int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

  // TaggedObject
  void Print(OPS_Stream& s, int flag = 0);

protected:
private:
  double computeCurrentStrain() const;
  double computeCurrentStrainRate() const;

  // private attributes - a copy for each object of the class
  UniaxialMaterial* theMaterial; // pointer to a material
  ID connectedExternalNodes;     // contains the tags of the end nodes
  int dimension;                 // truss in 2 or 3d domain
  int numDOF;                    // number of dof for truss

  Vector* theLoad;   // pointer to the load vector P
  Matrix* theMatrix; // pointer to objects matrix (a class wide Matrix)
  Vector* theVector; // pointer to objects vector (a class wide Vector)

  double L;              // length of truss based on undeformed configuration
  double A;              // area of truss
  double rho;            // rho: mass density per unit length
  int doRayleighDamping; // flag to include Rayleigh damping
  int cMass;             // consistent mass flag
  bool useInitialDisp;

  double cosX[3]; // direction cosines

  Node* theNodes[2];
  double* initialDisp;


  // Sensitivity
  int parameterID;
  Vector* theLoadSens;


  //
  // static data - single copy for all objects of the class
  //
  static Matrix trussM2;  // class wide matrix for 2*2
  static Matrix trussM4;  // class wide matrix for 4*4
  static Matrix trussM6;  // class wide matrix for 6*6
  static Matrix trussM12; // class wide matrix for 12*12
  static Vector trussV2;  // class wide Vector for size 2
  static Vector trussV4;  // class wide Vector for size 4
  static Vector trussV6;  // class wide Vector for size 6
  static Vector trussV12; // class wide Vector for size 12
};

#endif
