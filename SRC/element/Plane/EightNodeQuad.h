/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// based on FourNodeQuad by MHS
// Written: Seweryn Kokot, Opole University of Technology, Poland
// Created: Aug 2020
//
// Description: This file contains the class definition for EightNodeQuad.

#ifndef EightNodeQuad_h
#define EightNodeQuad_h

#ifndef _bool_h
#include <stdbool.h>
#endif

#include <Element.h>
#include <ID.h>
#include <Matrix.h>
#include <Vector.h>
#include <quadrature/Plane/LegendreFixedQuadrilateral.h>

class Node;
class NDMaterial;
class Response;

class EightNodeQuad : public Element,
                   protected LegendreFixedQuadrilateral<9>
{
public:
  EightNodeQuad(int tag, 
               std::array<int,8> &nodes, 
               NDMaterial &m, 
               double t,
               double pressure = 0.0, double rho = 0.0, double b1 = 0.0,
               double b2 = 0.0);
  EightNodeQuad();
  ~EightNodeQuad();

  const char *getClassType(void) const { return "EightNodeQuad"; };

  int getNumExternalNodes(void) const;
  const ID &getExternalNodes();   
  Node **getNodePtrs();   

  int getNumDOF();   
  void setDomain(Domain *theDomain);

  // public methods to set the state of the element
  int commitState();   
  int revertToLastCommit();   
  int revertToStart();   
  int update();   

  // public methods to obtain stiffness, mass, damping and residual information
  const Matrix &getTangentStiff();
  const Matrix &getInitialStiff();
  const Matrix &getMass();

  void zeroLoad();
  int addLoad(ElementalLoad *theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector &accel);

  const Vector &getResistingForce();   
  const Vector &getResistingForceIncInertia();   

  // public methods for element output
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  void Print(OPS_Stream &s, int flag = 0);

  Response *setResponse(const char **argv, int argc, OPS_Stream &s);

  int getResponse(int responseID, Information &eleInformation);

  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);

  // RWB; PyLiq1 & TzLiq1 need to see the excess pore pressure and initial
  // stresses.
  friend class PyLiq1;
  friend class TzLiq1;
  friend class QzLiq1; // Sumeet

private:
  // private attributes - a copy for each object of the class

//static constexpr int nip = 9;    // number of integration/Gauss points
  static constexpr int NEN = 8;    // number of nodes
  static constexpr int NDM = 2;

  NDMaterial **theMaterial; // pointer to the ND material objects

  ID connectedExternalNodes; // Tags of quad nodes

  Node *theNodes[NEN];

  static double matrixData[(NEN*2)*(NEN*2)]; // array data for matrix
  static Matrix K;              // Element stiffness, damping, and mass Matrix
  static Vector P;              // Element resisting force vector
  Vector Q;                     // Applied nodal loads
  double b[2];                  // Body forces

  double appliedB[2]; // Body forces applied with load pattern, C.McGann,
                      // U.Washington
  int applyLoad;      // flag for body force in load, C.McGann, U.Washington

  Vector pressureLoad; // Pressure load at nodes

  double thickness; // Element thickness
  double pressure;  // Normal surface traction (pressure) over entire element
                   // Note: positive for outward normal
  double rho;
  static double shp[3][NEN]; // Stores shape functions and derivatives (overwritten)
//static double pts[nip][2]; // Stores quadrature points
//static double wts[nip];    // Stores quadrature weights

  // private member functions
  double shapeFunction(double xi, double eta);
  void setPressureLoadAtNodes();   

  Matrix *Ki;
};

#endif
