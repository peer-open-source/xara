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

// Written: fmk
// Created: 07/98
// Revision: A
//
// Description: This file contains the class definition for Truss. A Truss
// object provides the abstraction of the small deformation bar element.
// Each truss object is associated with a section object. 
// This Truss element will work in 1d, 2d or 3d problems.
//
#ifndef TrussSection_h
#define TrussSection_h

#include <Element.h>
#include <Matrix.h>
#include <FrameSection.h>

class Node;
class Channel;
class FrameSection;


class TrussSection : public Element {
public:
  TrussSection(int tag, 
               int dimension, 
               int Nd1, int Nd2, 
               FrameSection& theSection,
               double rho = 0.0, 
               int doRayleighDamping = 0, 
               int cMass = 0);

  TrussSection();
  ~TrussSection();

  const char*
  getClassType() const override
  {
    return "TrussSection";
  }

  // public methods to obtain information about dof & connectivity
  int getNumExternalNodes() const override;
  const ID& getExternalNodes() override;
  Node** getNodePtrs() override;

  int getNumDOF() override;
  void setDomain(Domain*);

  // Alter the state of the element
  int commitState() override;
  int revertToLastCommit() override;
  int revertToStart() override;
  int update() override;

  // Obtain stiffness, mass, damping and residual information
  const Matrix& getTangentStiff() override;
  const Matrix& getInitialStiff() override;
  const Matrix& getDamp() override;
  const Matrix& getMass() override;
  const Vector& getResistingForce() override;
  const Vector& getResistingForceIncInertia() override;

  void zeroLoad() override;
  int addLoad(ElementalLoad* theLoad, double loadFactor) override;
  int addInertiaLoadToUnbalance(const Vector& accel) override;

  // MovableObject interface
  int sendSelf(int commitTag, Channel&) override;
  int recvSelf(int commitTag, Channel&, FEM_ObjectBroker&) override;

  // TaggedObject interface
  void Print(OPS_Stream& s, int flag) override;

  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information& );

  // Sensitivity
  int addInertiaLoadSensitivityToUnbalance(const Vector& accel, bool tag);
  int setParameter(const char** argv, int argc, Parameter&);
  int updateParameter(int parameterID, Information& info);
  int activateParameter(int parameterID);
  const Vector& getResistingForceSensitivity(int gradNumber);
  const Matrix& getKiSensitivity(int gradNumber);
  const Matrix& getMassSensitivity(int gradNumber);
  int commitSensitivity(int gradNumber, int numGrads);


private:
  double computeCurrentStrain() const;

  // Layout of stress resultants
  static constexpr FrameStressLayout section_layout = {
    FrameStress::N,
  };

  // private attributes - a copy for each object of the class
  ID connectedExternalNodes; // contains the tags of the end nodes
  int dimension;             // truss in 2 or 3d domain
  int numDOF;                // number of dof for truss

  Vector* theLoad;   // pointer to the load vector P
  Matrix* theMatrix; // pointer to objects matrix (a static Matrix)
  Vector* theVector; // pointer to objects vector (a static Vector)

  double cosX[3]; // direction cosines

  double L;              // length of truss based on undeformed configuration
  double rho;            // mass density per unit length
  int doRayleighDamping; // flag to include Rayleigh damping
  int cMass;             // consistent mass flag

  Node*   theNodes[2];
  double* initialDisp;

  FrameSection* theSection;

  // Sensitivity
  int parameterID;
  Vector* theLoadSens;

  // static data - single copy for all objects of the class
  static Matrix trussM2;  // class wide matrix for 2*2
  static Matrix trussM3;  // class wide matrix for 3*3
  static Matrix trussM4;  // class wide matrix for 4*4
  static Matrix trussM6;  // class wide matrix for 6*6
  static Matrix trussM12; // class wide matrix for 12*12
  static Vector trussV2;  // class wide Vector for size 2
  static Vector trussV3;  // class wide Vector for size 3
  static Vector trussV4;  // class wide Vector for size 4
  static Vector trussV6;  // class wide Vector for size 6
  static Vector trussV12; // class wide Vector for size 12
};

#endif
