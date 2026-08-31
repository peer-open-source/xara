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

#ifndef CorotTrussSection_h
#define CorotTrussSection_h

// Written: MHS
// Created: May 2001
//
// Description: This file contains the class definition for CorotTrussSection,
// a small strain, large displacement corotational space truss element,
// as described by Crisfield in [1]
//
//
// [1] Crisfield, M.A. "Nonlinear Finite Element Analysis of Solids and Structures", Vol. 1, 1991, J.T. Wiley.
//
#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <Matrix3D.h>
#include <FrameSection.h>
#include <Node.h>

class Node;

class CorotTrussSection : public Element {
public:
  CorotTrussSection(int tag, int dim, int Nd1, int Nd2, FrameSection&,
                    double rho = 0.0, int doRayleighDamping = 0, int cMass = 0);

  ~CorotTrussSection();

  const char*
  getClassType() const
  {
    return "CorotTrussSection";
  }

  // public methods to obtain information about dof & connectivity
  int getNumExternalNodes() const override;
  const ID& getExternalNodes();
  Node** getNodePtrs();

  int getNumDOF() override;
  void setDomain(Domain* theDomain);

  // public methods to set the state of the element
  int commitState();
  int revertToLastCommit();
  int revertToStart();
  int update();

  // public methods to obtain stiffness, mass, damping and residual information
  const Matrix& getTangentStiff();
  const Matrix& getInitialStiff();
  const Matrix& getDamp();
  const Matrix& getMass();

  void zeroLoad();
  int addLoad(ElementalLoad* theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector& accel);

  const Vector& getResistingForce();
  const Vector& getResistingForceIncInertia();

  // Sensitivity 
  int setParameter(const char** argv, int argc, Parameter& param) final;
  int updateParameter(int parameterID, Information& info) final;
  int activateParameter(int param) final;

  // public methods for element output
  void Print(OPS_Stream& s, int flag) final;

  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information& eleInformation);


private:
  double computeCurrentStrain();

  // Layout of stress resultants
  static constexpr FrameStressLayout section_layout = {
    FrameStress::N,
  };

  // private attributes - a copy for each object of the class
  FrameSection* theSection; // material (section) response
  ID connectedExternalNodes;           // contains the tags of the end nodes
  int numDOF;                          // number of dof for CorotTrussSection
  int numDIM;                          // number of dimensions

  double Lo;             // initial length of truss
  double Ln;             // current length of truss
  Vector3D d21;          // current displacement offsets in basic system
  double rho;            // mass density per unit length
  int doRayleighDamping; // flag to include Rayleigh damping
  int cMass;             // consistent mass flag

  Node* theNodes[2];

  Matrix3D R; // Rotation matrix
  int parameterID;

  Vector* theLoad;   // pointer to the load vector P
  Matrix* theMatrix; // pointer to objects matrix (a class wide Matrix)
  Vector* theVector; // pointer to objects vector (a class wide Vector)

  static Matrix M2;
  static Matrix M4;
  static Matrix M6;
  static Matrix M12;

  static Vector V2;
  static Vector V4;
  static Vector V6;
  static Vector V12;
};

#endif
