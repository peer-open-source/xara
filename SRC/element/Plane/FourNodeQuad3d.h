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
// Description: This file contains the class definition for FourNodeQuad3d.
//
// Written: MHS
// Created: Feb 2000
// Revised: Dec 2000 for efficiency
//
#ifndef FourNodeQuad3d_h
#define FourNodeQuad3d_h

#ifndef _bool_h
#include <stdbool.h>
#endif

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class Node;
class NDMaterial;
class Response;

class FourNodeQuad3d : public Element
{
  public:
    FourNodeQuad3d(int tag, int nd1, int nd2, int nd3, int nd4,
		   NDMaterial &m, const char *type,
		   double t, double pressure = 0.0, 
		   double rho = 0.0,
		   double b1 = 0.0, double b2 = 0.0);
    FourNodeQuad3d();
    ~FourNodeQuad3d();

    const char *getClassType() const {return "FourNodeQuad3d";}
    static constexpr const char* class_name = "FourNodeQuad3d";

    int getNumExternalNodes() const;
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
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker 
		  &theBroker);
    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(const char **argv, int argc, 
			  OPS_Stream &s);

    int getResponse(int responseID, Information &eleInformation);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

    // RWB; PyLiq1 & TzLiq1 need to see the excess pore pressure and initial stresses.
    friend class PyLiq1;
    friend class TzLiq1;
    friend class QzLiq1; // Sumeet

  protected:
    
  private:
    constexpr static int NEN = 4; // number of nodes
    constexpr static int NDF = 3; // number of DOFs per node

    // private attributes - a copy for each object of the class

    NDMaterial **theMaterial; // pointer to the ND material objects
    
    ID connectedExternalNodes; // Tags of quad nodes

    Node *theNodes[NEN];

    static double matrixData[144];  // array data for matrix
    static Matrix K;		// Element stiffness, damping, and mass Matrix
    static Vector P;		// Element resisting force vector
    Vector Q;		        // Applied nodal loads
    double b[2];		// Body forces

    Vector pressureLoad;	// Pressure load at nodes

    double thickness;	        // Element thickness

    double appliedB[2]; // Body forces applied with load pattern, C.McGann, U.Washington
    int applyLoad;      // flag for body force in load, C.McGann, U.Washington
	
    double pressure;	        // Normal surface traction (pressure) over entire element
					 // Note: positive for outward normal
    double rho;
    static double shp[3][4];	// Stores shape functions and derivatives (overwritten)
    static double pts[4][2];	// Stores quadrature points
    static double wts[4];	// Stores quadrature weights

    // private member functions - only objects of this class can call these
    double shapeFunction(double xi, double eta);
    void setPressureLoadAtNodes();

    int dirn[2];
};

#endif

