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
#ifndef ZeroLength_h
#define ZeroLength_h

// File: ~/element/zeroLength/ZeroLength.h
// 
// Written: GLF
// Created: 12/99
// Revision: A
//
// Description: This file contains the class definition for ZeroLength.
// A ZeroLength element is defined by two nodes with the same coordinate.
// One or more material objects may be associated with the nodes to
// provide a force displacement relationship.
// ZeroLength element will work with 1d, 2d, or 3d material models.
//
// What: "@(#) ZeroLength.h, revA"

#include <Element.h>
#include <Matrix.h>
#include <Matrix3D.h>




class Node;
class Channel;
class UniaxialMaterial;
class Response;

namespace OpenSees {
class ZeroLength : public Element
{

  public:
  // Tolerance for zero length of element
  static constexpr double LENTOL = 1.0e-6;
  static constexpr double MaxLength = 1e-6;
  // Type of dimension of element NxDy has dimension x=1,2,3 and
  // y=2,4,6,12 degrees-of-freedom for the element
  enum Etype { D1N2, D2N4, D2N6, D3N6, D3N12 };

  // Constructor for a multiple 1d material models
  ZeroLength(int tag,                               
             int dimension,
             int Nd1, int Nd2, 
             const Matrix3D& T,
             int n1dMat,
             UniaxialMaterial** theMaterial,  
             const ID& direction,
             int doRaylieghDamping = 0);

  // Constructor for a multiple 1d material models
  ZeroLength(int tag,                               
             int dimension,
             int Nd1, int Nd2, 
             const Matrix3D& T,
             int n1dMat,
             UniaxialMaterial** theMaterial,  
             UniaxialMaterial** theDampMaterial,  
             const ID& direction,
             int doRaylieghDamping = 0);

    ZeroLength();    
    ~ZeroLength();

    const char *getClassType() const {return "ZeroLength";}

    // public methods to obtain information about dof & connectivity    
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

    // public methods to obtain stiffness, mass, damping and residual information    
    const Matrix &getTangentStiff();
    const Matrix &getInitialStiff();
    const Matrix &getDamp();
    const Matrix &getMass();

    void zeroLoad();        
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);    

    const Vector &getResistingForce();
    const Vector &getResistingForceIncInertia();            

    // public methods for element output
    int sendSelf(int commitTag, Channel &);
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);
    // TaggedObject
    void Print(OPS_Stream &s, int flag);    

    Response *setResponse(const char **argv, int argc, OPS_Stream &);
    int getResponse(int responseID, Information &eleInformation);

    int setParameter(const char **argv, int argc, Parameter &);

    // Sensitivity
    const Vector &getResistingForceSensitivity(int gradIndex);
    int commitSensitivity(int gradIndex, int numGrads);


    void onActivate();
    void onDeactivate();


  protected:
    
  private:
    Etype elemType;

    // private methods
    void   setUp(int Nd1, int Nd2);
    void   checkDirection (  ID& dir ) const;
    
    void   setTran1d (Etype e, int n );
    double computeCurrentStrain1d ( int mat, const Vector& diff ) const;    

    // private attributes - a copy for each object of the class
    ID  connectedExternalNodes;         // contains the tags of the end nodes
    int dimension;                      // = 1, 2, or 3 dimensions
    int numDOF;                                // number of dof for ZeroLength
    Matrix3D transformation;                // transformation matrix for orientation
    int useRayleighDamping;
        
    Node *theNodes[2];

    Matrix *theMatrix;                     // pointer to objects matrix (a class Matrix)
    Vector *theVector;              // pointer to objects vector (a class Vector)

    // Storage for uniaxial material models
    int numMaterials1d;                           // number of 1d materials
    UniaxialMaterial **theMaterial1d;      // array of pointers to 1d materials
    ID               *dir1d;                // array of directions 0-5 for 1d materials
    Matrix           *t1d;            // hold the transformation matrix

    // vector pointers to initial disp and vel if present
    Vector *d0;
    Vector *v0;

    // static data - single copy for all objects of the class        
    static Matrix ZeroLengthM2;   // class wide matrix for 2*2
    static Matrix ZeroLengthM4;   // class wide matrix for 4*4
    static Matrix ZeroLengthM6;   // class wide matrix for 6*6
    static Matrix ZeroLengthM12;  // class wide matrix for 12*12
    static Vector ZeroLengthV2;   // class wide Vector for size 2
    static Vector ZeroLengthV4;   // class wide Vector for size 4
    static Vector ZeroLengthV6;   // class wide Vector for size 6
    static Vector ZeroLengthV12;  // class wide Vector for size 12

    int mInitialize;  // tag to fix bug in recvSelf/setDomain when using database command
};
} // end of namespace OpenSees
#endif




