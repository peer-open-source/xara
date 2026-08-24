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
// by Jinchi Lu and Zhaohui Yang (May 2004)
//
// 20NodeBrick element
//
#pragma once

#include <array>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <NDMaterial.h>

class Twenty_Node_Brick : public Element {

public:
  Twenty_Node_Brick();

  Twenty_Node_Brick(int tag, 
                    const std::array<int, 20>& node_tags,
                    NDMaterial& theMaterial, 
                    double b1 = 0.0, 
                    double b2 = 0.0,
                    double b3 = 0.0);

  virtual ~Twenty_Node_Brick();

  const char*
  getClassType() const
  {
    return "Twenty_Node_Brick";
  }

  void setDomain(Domain*);
  int getNumExternalNodes() const;

  // return connected external nodes
  const ID& getExternalNodes();
  Node** getNodePtrs();

  int getNumDOF();

  int commitState();
  int revertToLastCommit();
  int revertToStart();

  int update();

  const Matrix& getTangentStiff();
  const Matrix& getInitialStiff();
  const Matrix& getDamp();
  const Matrix& getMass();

  void zeroLoad();
  int addLoad(ElementalLoad* theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector& accel);

  const Vector& getResistingForce();
  const Vector& getResistingForceIncInertia();

  void Print(OPS_Stream& s, int flag);

  // public methods for element output
  int sendSelf(int commitTag, Channel&);
  int recvSelf(int commitTag, Channel&, FEM_ObjectBroker&);

  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information&);

private:
  void formInertiaTerms(int tangFlag);
  void formDampingTerms(int tangFlag);

  // Mixture mass density at integration point i
  double mixtureRho(int ipt);
  void computeBasis();
  void compuLocalShapeFunction();
  void Jacobian3d(int gaussPoint, double& xsj, int mode);
  const Matrix& getStiff(int flag);

  constexpr static int NEN = 20;
  // quadrature data
  static constexpr int nintu = 27;


  // static data
  static Matrix stiff;
  static Vector resid;
  static Matrix mass;
  static Matrix damp;


  // node information
  ID connectedExternalNodes;  // node tags
  Node* nodePointers[NEN];    // pointers to nodes

  std::array<NDMaterial*, nintu> materialPointers;

  //local nodal coordinates, three coordinates for each of twenty nodes
  //    static double xl[3][20] ;
  static double xl[3][NEN];
  double b[3]; // Body forces

  double appliedB[3]; // Body forces applied with load pattern, C.McGann, U.Washington
  int applyLoad;      // flag for body force in load, C.McGann, U.Washington

  static double shgu[4][NEN][27]; // Stores shape functions and derivatives (overwritten)
  static double shlu[4][NEN][27]; // Stores shape functions and derivatives
  static double wu[27];          // Stores quadrature weights
  static double dvolu[27];       // Stores detJacobian (overwritten)

  Vector* load;
  Matrix* Ki;

  // compute local shape functions
};
