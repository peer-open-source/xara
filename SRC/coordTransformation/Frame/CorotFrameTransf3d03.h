//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the class definition for
// CorotFrameTransf3d03.h. CorotFrameTransf3d03 provides the
// abstraction of a corotational transformation for a spatial frame element
//
// Written by : cmp, March 2024
//
// Adapted from work by: Remo Magalhaes de Souza
//
#ifndef CorotFrameTransf3d03_h
#define CorotFrameTransf3d03_h

#include <array>
#include "FrameTransform.h"
#include <Vector.h>
#include <Matrix.h>
#include <Versor.h>
#include <Matrix3D.h>
#include <Vector3D.h>
#include "Orient/CrisfieldTransform.h"

struct Triad;
using namespace OpenSees; // TODO

class CorotFrameTransf3d03: public FrameTransform3d
{
public:
    CorotFrameTransf3d03(int tag, const Vector &vecInLocXZPlane,
                       const Vector &rigJntOffsetI, 
                       const Vector &rigJntOffsetJ);

    CorotFrameTransf3d03();
    ~CorotFrameTransf3d03();

    const char *getClassType() const {
      return "CorotFrameTransf3d";
    }

    // NOTE: maybe add arg for rotation parameterization
    FrameTransform3d *getCopy();

    int initialize(Node *nodeI, Node *nodeJ);
    int update();
    int commitState();
    int revertToLastCommit();        
    int revertToStart();
    int getLocalAxes(Vector &xAxis, Vector &yAxis, Vector &zAxis);

    double getInitialLength();
    double getDeformedLength();
//  double getPresentLength();
//  double getCurrentLength();
    const Vector &getBasicTrialDisp();
    const Vector &getBasicIncrDisp();
    const Vector &getBasicIncrDeltaDisp();
    const Vector &getBasicTrialVel();
    const Vector &getBasicTrialAccel();

    virtual VectorND<12>    pushResponse(VectorND<12>&pl) override final;
    virtual VectorND<12>    pushConstant(const VectorND<12>&pl) const override final;

    virtual MatrixND<12,12> pushResponse(MatrixND<12,12>& kl, const VectorND<12>& pl) override final;
    virtual MatrixND<12,12> pushConstant(const MatrixND<12,12>& kl) override final;

    const Vector &getGlobalResistingForce(const Vector &qb, const Vector &p0);
    const Matrix &getGlobalStiffMatrix(const Matrix &Kb, const Vector &qb);
    const Matrix &getInitialGlobalStiffMatrix(const Matrix &Kb);

    // method used to rotate consistent mass matrix
    const Matrix &getGlobalMatrixFromLocal(const Matrix &Ml);

    // methods used in post-processing only
    const Vector &getPointGlobalCoordFromLocal(const Vector &localCoords);
    const Vector &getPointGlobalDisplFromBasic(double xi, const Vector &basicDisps);
    const Vector &getPointLocalDisplFromBasic(double xi, const Vector &basicDisps);

    // Sensitivity
    double getLengthGrad();

    // Movable Object
    int sendSelf(int tag, Channel &);
    int recvSelf(int tag, Channel &, FEM_ObjectBroker &);

    // Tagged Object
    void Print(OPS_Stream &s, int flag = 0);

protected:
    int addTangent(MatrixND<12,12>& M, const VectorND<12>& pl);

    VectorND<6>   pushResponse(const VectorND<6>& pa, int a, int b);
    MatrixND<6,6> pushResponse(const MatrixND<6,6>& K, const VectorND<12>& pl, int a, int b);
    int addTangent(MatrixND<6,6>& K, const VectorND<6>& p, int a, int b, int c);

protected:

private:
    // compute the transformation matrix
    void compTransfMatrixBasicGlobal(const Versor&, const Versor* Q);

    template<typename VecL, typename VecB>
    void
    LocalToBasic(const VecL& ul, VecB& ub)
    {
      ub[0] =  ul[jnx] - ul[inx];
      ub[1] =  ul[imz];
      ub[2] =  ul[jmz];
      ub[3] =  ul[imy];
      ub[4] =  ul[jmy];
      ub[5] =  ul[jmx] - ul[imx];
    }

    constexpr static int nn = 2;
    enum {
      inx= 0, // axial
      iny= 1, // Vy
      inz= 2, // Vz
      imx= 3, // torsion
      imy= 4, // rot y I
      imz= 5, // rot z I

      jnx= 6, // axial
      jny= 7,
      jnz= 8,
      jmx= 9, // torsion
      jmy=10, // rot y J
      jmz=11, // rot z J
    };

    //
    // Member data
    //
    std::array<Node*, nn> nodes;

    Vector3D xAxis;                              // local x axis
    Vector3D vz;                                 // Vector that lies in local plane xz
    Vector3D dX;

    // Rigid joint offsets
    enum {
      end_i = 0<<1,
      end_j = 0<<2,
    };
    int joint_offsets;

    Vector nodeIOffset, 
           nodeJOffset;
    Vector3D* offset[nn];

    double *nodeIInitialDisp, *nodeJInitialDisp;
    bool  initialDispChecked;                    

    double L;                        // initial element length
    double Ln;                       // current element length (at trial state)

    Versor Q_past[nn];                // commited rotations
    Versor Q_pres[nn];                // trial rotations

    Vector3D alphaI;                 // last trial rotations end i
    Vector3D alphaJ;                 // last trial rotatations end j

    VectorND<12> ul;                 // local displacements (size=7)
    Vector3D     vr[nn];             //
    VectorND<12> ulcommit;           // commited local displacements
    VectorND<12> ulpr;               // previous local displacements

    OpenSees::MatrixND<12,12> T;     // transformation from local to global system

    OpenSees::Matrix3D R0;           // rotation from local to global coordinates
    CrisfieldTransform crs;
    constexpr static Vector3D E1 {1, 0, 0},
                              E2 {0, 1, 0},
                              E3 {0, 0, 1};

    // Static workspace variables
    Matrix3D A;
    MatrixND<12,3> Lr2, Lr3;   // auxiliary matrices
};

#endif
