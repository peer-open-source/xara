/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// Description: This file contains the class definition for LagrangeQuad.
//
// Written: MHS
// Sensitivity by Quan Gu, Michele Barbato, Joel P. Conte @ UCSD.  2009 July.
//
#ifndef LagrangeQuad_h
#define LagrangeQuad_h

#include <array>
#include <Mate.h>
#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <quadrature/GaussLegendre2D.hpp>

#define ELE_TAG_LagrangeQuad 0

class Node;
class Response;

namespace OpenSees {

template<int NEN, bool enhanced=false>
class LagrangeQuad : public Element, protected GaussLegendre<2, 4> {
public:
  LagrangeQuad(int tag, const std::array<int, NEN>& nodes, 
               Mate<2>& m,  
               double thickness,
               double pressure = 0.0, 
               double rho = 0.0, 
               double b1 = 0.0, double b2 = 0.0);

  LagrangeQuad();
  ~LagrangeQuad();

  const char*
  getClassType() const
  {
    return "LagrangeQuad";
  }
  static constexpr const char* class_name = "LagrangeQuad";

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

  // Public methods for element output

  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information& eleInformation);

  // Inherited from TaggedObject
  void Print(OPS_Stream& s, int flag = 0);

  // Inherited from MovableObject
  int sendSelf(int tag, Channel&);
  int recvSelf(int tag, Channel&, FEM_ObjectBroker&);

  // Sensitivity
#if 0
    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);
    int            activateParameter           (int parameterID);
    const Vector & getResistingForceSensitivity(int gradNumber);
//  const Matrix & getKiSensitivity            (int gradNumber);
//  const Matrix & getMassSensitivity          (int gradNumber);
    int            commitSensitivity           (int gradNumber, int numGrads);
#endif

protected:
private:
  constexpr static int NDM = 2; // number of spatial dimensions
  constexpr static int NDF = 2; // number of DOFs per node
  constexpr static int NIP = 4; // number of integration points

  //
  // private member functions
  // - only objects of this class can call these
  double shapeFunction(double xi, double eta);
  void setPressureLoadAtNodes();

  //
  // private attributes
  // - a copy for each object of the class
  std::array<Mate<2>*, NIP> theMaterial; // pointer to the ND material objects

  ID connectedExternalNodes; // Tags of the nodes

  std::array<Node*, NEN> theNodes;

  Matrix* Ki;
//double matrixData[64]; // array data for matrix
//Matrix K;              // Element stiffness, damping, and mass Matrix
//Vector P;              // Element resisting force vector
  VectorND<NEN*NDF> Q;   // Applied nodal loads
  double b[2];           // Body forces

  double appliedB[2]; // Body forces applied with load pattern, C.McGann, U.Washington
  int applyLoad;      // flag for body force in load

  Vector pressureLoad; // Pressure load at nodes

  double thickness; // Element thickness
  double pressure;  // Normal surface traction (pressure) over entire element
                    // Note: positive for outward normal
  double rho;
  double shp[3][NEN]; // shape functions and derivatives

  int parameterID;
};
} // namespace OpenSees

#include "LagrangeQuad.tpp"
#endif
