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
                                                                       
#ifndef SSPquad_h
#define SSPquad_h

// Created: C.McGann, UW, 04.2011
//
// Description: This file contains the class definition for SSPquad (Stabilized Single-Point Quad)

#include <Element.h>
#include <Node.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <array>

// number of nodes per element
#define SSPQ_NUM_NODE 4
// number of dimensions
#define SSPQ_NUM_DIM  2
// degrees of freedom per element
#define SSPQ_NUM_DOF  8

class Domain;
class Node;
class Channel;
class NDMaterial;
class FEM_ObjectBroker;
class Response;

namespace OpenSees {

class PyLiq1;
class TzLiq1;
class QzLiq1;

class SSPquad : public Element
{
  public:
    SSPquad(int tag,
            const std::array<int, 4> &nodes,
            NDMaterial &theMat,
            double thick, 
            double density,
            double b1 = 0.0, double b2 = 0.0);
    SSPquad();
    ~SSPquad();

    const char* getClassType()  const { return "SSPquad"; }

    // public methods to obtain information about dof and connectivity
    int getNumExternalNodes() const; 
    const ID &getExternalNodes();
    Node **getNodePtrs();
    int getNumDOF();
    void setDomain(Domain *);

    // public methods to set the state of the element
    int commitState();
    int revertToLastCommit();
    int revertToStart();
    int update();

    // public methods to obtain stiffness, mass, damping, and residual info
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

    Response *setResponse(const char **argv, int argc, OPS_Stream &);
    int getResponse(int responseID, Information &);

    // public methods for material stage update
    int setParameter(const char **argv, int argc, Parameter &);
    int updateParameter(int parameterID, Information &);

    // allow PyLiq1 and TzLiq1 classes to get stresses from SSPquadUP class
    friend class PyLiq1;
    friend class TzLiq1;
    friend class QzLiq1; // Sumeet

  protected:

  private:

    // member functions
    void GetStab();                                     // compute stabilization stiffness matrix

    // objects
    NDMaterial *theMaterial;                            // pointer to NDMaterial object
    ID mExternalNodes;                                  // contains tags of the nodes
    Matrix mTangentStiffness;                           // tangent stiffness matrix
    Vector mInternalForces;                             // vector of internal forces
    Vector Q;                                           // vector of applied nodal forces
    Matrix mMass;                                       // mass matrix

    Node *theNodes[4];

    // input quantities
    double mThickness;                                  // thickness of element
    double b[2];                                        // body forces acting on element

    // load pattern variables
    double appliedB[2];                                 // body forces applied with load pattern
    int applyLoad;                                      // flag for body force in load pattern

    // calculation variables
    double J0;                                          // constant portion of jacobian
    double J1;                                          // linear (xi) portion of jacobian
    double J2;                                          // linear (eta) portion of jacobian
    
    Matrix Mmem;                                        // mapping matrix for membrane modes
    Matrix Kstab;                                       // stabilization stiffness matrix
    Matrix mNodeCrd;                                    // nodal coordinate array
    double density;                                     // mass density
};
} // namespace OpenSees
#endif
