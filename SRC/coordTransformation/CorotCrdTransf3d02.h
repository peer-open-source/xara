//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//    See https://opensees.berkeley.edu/OpenSees/copyright.php for license.
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the class definition for
// CorotCrdTransf3d02.h. CorotCrdTransf3d02 provides the
// abstraction of a corotation transformation for a spatial frame element
//
// Written: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 04/2000
// Modified by : cmp, March 2024
//
#ifndef CorotCrdTransf3d02_h
#define CorotCrdTransf3d02_h

#include "CrdTransf.h"
#include <Vector.h>
#include <Matrix.h>
#include <Matrix3D.h>

class Triad;
namespace OpenSees {class Matrix3D;}

class CorotCrdTransf3d02: public CrdTransf
{
public:
    CorotCrdTransf3d02(int tag, const Vector &vecInLocXZPlane,
        const Vector &rigJntOffsetI, const Vector &rigJntOffsetJ);

    CorotCrdTransf3d02();
    ~CorotCrdTransf3d02();

    const char *getClassType() const {return "CorotCrdTransf3d02";}

    int initialize(Node *nodeIPointer, Node *nodeJPointer);
    int update();       // Set RI,RJ,Rbar, Ln, e and ul
    int commitState();
    int revertToLastCommit();        
    int revertToStart();

    double getInitialLength();
//  double getCurrentLength();
    double getDeformedLength();
    const Vector &getBasicTrialDisp();
    const Vector &getBasicIncrDisp();
    const Vector &getBasicIncrDeltaDisp();
    const Vector &getBasicTrialVel();
    const Vector &getBasicTrialAccel();

    const Vector &getGlobalResistingForce(const Vector &basicForce, const Vector &p0);
    const Matrix &getGlobalStiffMatrix(const Matrix &basicStiff, const Vector &basicForce);
    const Matrix &getInitialGlobalStiffMatrix(const Matrix &basicStiff);

    CrdTransf *getCopy3d() final;

    // method used to rotate consistent mass matrix
    const Matrix &getGlobalMatrixFromLocal(const Matrix &local);

    // methods used in post-processing only
    const Vector &getPointGlobalCoordFromLocal(const Vector &localCoords);
    const Vector &getPointGlobalDisplFromBasic(double xi, const Vector &basicDisps);
    const Vector &getPointLocalDisplFromBasic(double xi, const Vector &basicDisps);    
    
    int  getLocalAxes(Vector &xAxis, Vector &yAxis, Vector &zAxis);

    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag);
    
private:
    // compute the transformation matrix
    void compTransfMatrixBasicGlobal(const Triad&  r, const Triad&  E, const Triad&  rI, const Triad&  rJ);
    void compTransfMatrixLocalGlobal(Matrix &Tlg);
    void compTransfMatrixBasicLocal(Matrix &Tbl);


    //
    // internal data
    //
    Node *nodeIPtr, *nodeJPtr;  // pointers to the element two endnodes

    Vector xAxis;               // local x axis
    Vector vAxis;               // Vector that lies in local plane xz
    Vector nodeIOffset, 
           nodeJOffset;         // rigid joint offsets
    
    double L;                   // undeformed element length
    double Ln;                  // deformed element length
    
                                // (the columns of which are the element local axes)
    OpenSees::VectorND<4> alphaIq;             // quaternion for node I
    OpenSees::VectorND<4> alphaJq;             // quaternion for node I 
    OpenSees::VectorND<4> alphaIqcommit;       // commited quaternion for node I
    OpenSees::VectorND<4> alphaJqcommit;       // commited quaternion for node J
    Vector alphaI;              // last trial rotations end i
    Vector alphaJ;              // last trial rotatations end j

    Vector ul;                  // local displacements
    Vector ulcommit;            // commited local displacements
    Vector ulpr;                // previous local displacements

    // previously static variables
    Matrix T;                   // transformation matrix from basic to global system
    OpenSees::Matrix3D A;
    OpenSees::Matrix3D R0;      // rotation matrix from local to global coordinates
    OpenSees::Matrix3D e, RI, RJ, Rbar;
    
    double *nodeIInitialDisp, *nodeJInitialDisp;
    bool initialDispChecked;

    // static workspace variables
    static Matrix Tp;           // transformation matrix to renumber dofs
    static Matrix Tbl;          // transformation matrix from local to basic system
    static Matrix kg;           // global stiffness matrix
    static Matrix Lr2, Lr3;     // auxiliary matrices
//  static Matrix Tlg;          // transformation matrix from global to local system
//  static Matrix TlgInv;       // inverse of transformation matrix from global to local system
};

#endif
