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
// Written: Roozbeh Geraili Mikola (roozbehg@berkeley.edu)
// Created: Sep 2010
//
// Description: This file contains the class definition for Tri31.
//
#ifndef Tri31_h
#define Tri31_h

#ifndef _bool_h
#include <stdbool.h>
#endif

#include <array>
#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class Node;
class NDMaterial;
class Response;

class Tri31 : public Element
{
  public:
    Tri31(int tag, 
          std::array<int,3> &nodes,
          NDMaterial &m,
          double thickness, 
          double pressure,
          double rho,
          double b1, double b2);
    Tri31();
    ~Tri31();

    const char *getClassType() const {
      return "Tri31";
    }
    static constexpr const char* class_name = "Tri31";

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
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);

    int getResponse(int responseID, Information &eleInformation);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

    // RWB; PyLiq1 & TzLiq1 need to see the excess pore pressure and initial stresses.
    friend class PyLiq1;
    friend class TzLiq1;
    friend class QzLiq1; // Sumeet

  protected:
    
  private:

    constexpr static int NDM = 2;      // number of dimensions
    constexpr static int NEN = 3;      // number of nodes
    constexpr static int NDF = 2;      // number of nodes
    constexpr static int NIP = 1;      // number of gauss points


    std::array<NDMaterial *,NIP> theMaterial; // array of ND material objects
    
    ID connectedExternalNodes; // Tags of Tri31 nodes

    std::array<Node *, NEN> theNodes;

    static double matrixData[36];  // array data for matrix
    static Matrix K;               // Element stiffness, damping, and mass Matrix
    static Vector P;               // Element resisting force vector
    Vector Q;                   // Applied nodal loads
    double b[2];                // Body forces

    double appliedB[2];         // Body forces applied with load pattern
    int applyLoad;              // flag for body force in load

    Vector pressureLoad;        // Pressure load at nodes

    double thickness;           // Element thickness
    double pressure;            // Normal surface traction (pressure) over entire element
                                // Note: positive for outward normal
    double rho;
    static double shp[3][3];    // Stores shape functions and derivatives (overwritten)
    static double pts[1][2];    // Stores quadrature points
    static double wts[1];       // Stores quadrature weights

    // private member functions - only objects of this class can call these
    double shapeFunction(double xi, double eta);
    void setPressureLoadAtNodes();

    Matrix *Ki;
};

#endif

