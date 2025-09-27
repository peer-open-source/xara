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

#ifndef TwoNodeLinkSection_h
#define TwoNodeLinkSection_h

#include "Element.h"
#include <Matrix.h>
#include <cassert>

class Channel;
class SectionForceDeformation;
class Response;

class TwoNodeLinkSection : public Element
{
public:
  TwoNodeLinkSection(int tag, int dimension, int Nd1, int Nd2,
          SectionForceDeformation &section,
          const Vector y = 0, const Vector x = 0,
          const Vector Mratio = 0, const Vector shearDistI = 0,
          int addRayleigh = 0, double mass = 0.0);
  TwoNodeLinkSection();
  ~TwoNodeLinkSection();

  const char *getClassType() const override {return "TwoNodeLinkSection";}
  
  // public methods to obtain information about dof & connectivity
  int getNumExternalNodes() const;
  const ID &getExternalNodes();
  Node **getNodePtrs();
  int getNumDOF();
  void setDomain(Domain *);
  
  // public methods to set the state of the element
  int commitState() override;
  int revertToLastCommit() override;
  int revertToStart() override;
  int update() override;

  const Matrix &getTangentStiff() override;
  const Matrix &getInitialStiff() override;
  const Matrix &getDamp() override;
  const Matrix &getMass() override;

  void zeroLoad();
  int addLoad(ElementalLoad *, double loadFactor) override;
  int addInertiaLoadToUnbalance(const Vector &accel) override;
  
  const Vector &getResistingForce() override;
  const Vector &getResistingForceIncInertia();
  
  // public methods for element output
  int sendSelf(int commitTag, Channel &) override;
  int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &) override;
  void Print(OPS_Stream &s, int flag) override;
  
  // public methods for element recorder
  Response *setResponse(const char **argv, int argc, OPS_Stream &);
  int getResponse(int responseID, Information &);
  int setParameter(const char **argv, int argc, Parameter &);

private:
  // Type of dimension of element NxDy has dimension x=1,2,3 and
  // y=2,4,6,12 degrees-of-freedom for the element
  enum Etype { D1N2, D2N4, D2N6, D3N6, D3N12 };

  constexpr Etype SetType(int ndm, int ndf) {
    if (ndm == 1 && ndf == 1)
      return D1N2;
    else if (ndm == 2 && ndf == 2)
      return D2N4;
    else if (ndm == 2 && ndf == 3)
      return D2N6;
    else if (ndm == 3 && ndf == 3)
      return D3N6;
    else if (ndm == 3 && ndf == 6)
      return D3N12;

    assert(false); // Should never get here
    return D1N2;
  }
  Etype elemType;

  void setUp();
  void setTranGlobalLocal();
  void setTranLocalBasic();
  void addPDeltaForces(Vector &pLocal, const Vector &qBasic);
  void addPDeltaStiff(Matrix &kLocal, const Vector& qBasic);
  int getDirID(int sectionCode);

  // private attributes - a copy for each object of the class
  int numDIM;                         // 1, 2, or 3 dimensions
  int numDOF;                         // number of dof for TwoNodeLinkSection
  ID connectedExternalNodes;          // contains the tags of the end nodes
  Node *theNodes[2];                  // array of nodes
  SectionForceDeformation *theSection;

  // parameters
  Matrix trans;       // transformation matrix for orientation
  Vector x;           // local x direction
  Vector y;           // local y direction
  Vector Mratio;      // p-delta moment distribution ratios
  Vector shearDistI;  // shear distance from node I as fraction of length
  int addRayleigh;    // flag to add Rayleigh damping
  double mass;        // total mass
  double L;           // element length
  bool onP0;          // flag to indicate if the element is on P0
  
  Vector ub;          // trial displacements in basic system
  Vector ubdot;       // trial velocities in basic system
  Vector qb;          // resisting forces in basic system
  Vector ul;          // displacements in local system
  Matrix Tgl;         // transformation matrix from global to local system
  Matrix Tlb;         // transformation matrix from local to basic system
  
  Matrix *theMatrix;  // pointer to objects matrix (a class wide Matrix)
  Vector *theVector;  // pointer to objects vector (a class wide Vector)
  Vector *theLoad;    // pointer to the load vector
};

#endif
