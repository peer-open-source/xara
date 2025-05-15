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
// Description: This file contains the class definition for FourNodeQuad.
//
// Written: MHS
// Created: Feb 2000
// Revised: Dec 2000 for efficiency
// Sensitivity by Quan Gu, Michele Barbato, Joel P. Conte @ UCSD.  2009 July.
// 
#ifndef FourNodeQuad_h
#define FourNodeQuad_h

#include <array>
#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <quadrature/GaussLegendre2D.hpp>

class Node;
class NDMaterial;
class Response;

class FourNodeQuad : public Element,
                     protected GaussLegendre<2,4>
{
  public:
    FourNodeQuad(int tag,
                 std::array<int,4>& nodes,
                 NDMaterial &m, 
                 double thickness,
                 double p, double r, 
                 double b1, double b2);

    FourNodeQuad();
    ~FourNodeQuad();

    const char *getClassType() const {return "FourNodeQuad";}
    static constexpr const char* class_name = "FourNodeQuad";

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

    // Public methods for element output

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);

    // Inherited from TaggedObject
    void Print(OPS_Stream &s, int flag =0);

    // Inherited from MovableObject
    int sendSelf(int tag, Channel &);
    int recvSelf(int tag, Channel &, FEM_ObjectBroker &);

    // Sensitivity
    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);
    int            activateParameter           (int parameterID);
    const Vector & getResistingForceSensitivity(int gradNumber);
//  const Matrix & getKiSensitivity            (int gradNumber);
//  const Matrix & getMassSensitivity          (int gradNumber);
    int            commitSensitivity           (int gradNumber, int numGrads);

    // RWB; PyLiq1 & TzLiq1 need to see the excess pore pressure and initial stresses.
    friend class PyLiq1;
    friend class TzLiq1;
    friend class QzLiq1; // Sumeet

  protected:
    
  private:
    constexpr static int NDM = 2;    // number of spatial dimensions
    constexpr static int NEN = 4;    // number of nodes
    constexpr static int NDF = 2;    // number of DOFs per node
    constexpr static int NIP = 4;    // number of integration points
    constexpr static int NST = 3;    // number of stress components

    //
    // private member functions 
    //
    int stateDetermination(Matrix* Kptr, Vector* pptr, int flag);
    double shapeFunction(double xi, double eta);
    void setPressureLoadAtNodes();

    //
    // private attributes 
    //
    std::array<NDMaterial *, NIP> theMaterial; // pointers to materials
    
    ID connectedExternalNodes; // Tags of the nodes

    std::array<Node *, NEN> theNodes;

    Matrix *Ki;
    static double matrixData[64];   // array data for matrix
    static Matrix K;                // Element stiffness, damping, and mass Matrix
    static Vector P;                // Element resisting force vector
    Vector Q;                       // Applied nodal loads
    double b[2];                    // Body forces

    double appliedB[2];             // Body forces applied with load pattern, C.McGann, U.Washington
    int applyLoad;                  // flag for body force in load
        
    Vector pressureLoad;        // Pressure load at nodes

    double thickness;                // Element thickness
    double pressure;                 // Normal surface traction (pressure) over entire element
                                     // Note: positive for outward normal
    double rho;
    double shp[3][NEN];              // shape functions and derivatives

    int parameterID;

};

#endif

