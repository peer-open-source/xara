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
                                                                       
#ifndef SSPbrick_h
#define SSPbrick_h

// Created: C.McGann, UW, 10.2011
//
// Description: This file contains the class definition for SSPbrick (Stabilized Single-Point Brick)

#include <Element.h>
#include <Node.h>
#include <Vector.h>
#include <Matrix.h>
#include <VectorND.h>
#include <MatrixND.h>
#include <ID.h>

class Domain;
class Node;
class Channel;
class NDMaterial;
class FEM_ObjectBroker;
class Response;
using OpenSees::VectorND;
using OpenSees::MatrixND;

class SSPbrick : public Element
{
  public:
    SSPbrick(int tag, int Nd1, int Nd2, int Nd3, int Nd4, int Nd5, int Nd6, int Nd7, int Nd8,
                      NDMaterial &theMat, double b1 = 0.0, double b2 = 0.0, double b3 = 0.0);
    SSPbrick();
    ~SSPbrick();

	const char* getClassType()  const override {return "SSPbrick";}

    // public methods to obtain information about dof and connectivity
    int getNumExternalNodes(void) const; 
    const ID &getExternalNodes(void);
	Node **getNodePtrs(void);
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

	Response *setResponse(const char **argv, int argc, OPS_Stream &eleInfo);
	int getResponse(int responseID, Information &);

	// public methods for material stage update
	int setParameter(const char **argv, int argc, Parameter &);
    int updateParameter(int parameterID, Information &);

  protected:

  private:

	static constexpr int nen = 8,
	// number of nodes per element
	 SSPB_NUM_NODE = 8,
	// number of dimensions
	 SSPB_NUM_DIM  = 3,
	// degrees of freedom per element
	 SSPB_NUM_DOF  = 24;

    // member functions
	void GetStab();                                 // compute stabilization stiffness matrix
	Matrix Transpose(int d1, int d2, const Matrix &M);  // transpose operation

	// objects
	NDMaterial *theMaterial;                            // pointer to NDMaterial object
	ID mExternalNodes;                                  // contains tags of the nodes
	Matrix mTangentStiffness;                           // tangent stiffness matrix
	Vector mInternalForces;                             // vector of internal forces
	Vector Q;                                           // vector of applied nodal forces
	Matrix mMass;                                       // mass matrix

	Node *theNodes[8];

	// input quantities
	double b[3];                                        // body forces acting on element

	// load pattern variables
	double appliedB[3];                                 // body forces applied with load pattern
	int applyLoad;                                      // flag for body force in load pattern

	// calculation variables
	double J[20];                                       // jacobian determinant terms
	double mVol;                                        // element volume

	bool mInitialize;
	
	Matrix Bnot;                                        // mapping matrix for membrane modes
	Matrix Kstab;                                       // stabilization stiffness matrix
	MatrixND<SSPB_NUM_DIM,SSPB_NUM_NODE> mNodeCrd;      // nodal coordinate array
	
	VectorND<8> xi;                                          // xi evaluated at the nodes
	VectorND<8> et;                                          // eta evaluated at the nodes
	VectorND<8> ze;                                          // zeta evaluated at the nodes
	Vector hut;                                         // zeta*eta evaluated at the nodes
	Vector hus;                                         // zeta*xi evaluated at the nodes
	Vector hst;                                         // xi*eta evaluated at the nodes
	Vector hstu;                                        // xi*eta*zeta evaluated at the nodes
};

#endif
