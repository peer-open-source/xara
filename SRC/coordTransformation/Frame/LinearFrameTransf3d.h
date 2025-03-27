//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the class definition for
// LinearFrameTransf3d.h. LinearFrameTransf3d provides the
// abstraction of a linear transformation for a spatial frame
// between the global and basic coordinate systems
//
// Written: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 04/2000
//
#ifndef LinearFrameTransf3d_h
#define LinearFrameTransf3d_h

#include <array>
#include <FrameTransform.h>
#include <Vector.h>
#include <Matrix.h>

class LinearFrameTransf3d: public FrameTransform3d
{
public:
    LinearFrameTransf3d(int tag, const Vector &vecInLocXZPlane);
    LinearFrameTransf3d(int tag, const Vector &vecInLocXZPlane,
        const Vector &rigJntOffsetI,
        const Vector &rigJntOffsetJ);
    
    LinearFrameTransf3d();
    ~LinearFrameTransf3d();
    
    const char *getClassType() const {return "LinearFrameTransf3d";}
    
    virtual int getLocalAxes(Vector &xAxis, Vector &yAxis, Vector &zAxis);
    
    virtual FrameTransform3d *getCopy();

    virtual double getInitialLength();
    virtual double getDeformedLength();
    
    virtual int initialize(Node *node1Pointer, Node *node2Pointer) override final;
    virtual int update() override final;
    virtual int commitState() override final;
    virtual int revertToLastCommit() override final;
    virtual int revertToStart() override final;

    virtual const Vector &getBasicTrialDisp() override;
    virtual const Vector &getBasicIncrDisp();
    virtual const Vector &getBasicIncrDeltaDisp();
    virtual const Vector &getBasicTrialVel();
    virtual const Vector &getBasicTrialAccel();

    virtual VectorND<12>    pushResponse(VectorND<12>&pl) override final;
    virtual VectorND<12>    pushConstant(const VectorND<12>&pl) const override final;

    virtual MatrixND<12,12> pushResponse(MatrixND<12,12>& kl, const VectorND<12>& pl) override final;
    virtual MatrixND<12,12> pushConstant(const MatrixND<12,12>& kl) override final;

    virtual const Vector &getGlobalResistingForce(const Vector &basicForce, const Vector &p0);
    virtual const Matrix &getGlobalStiffMatrix(const Matrix &basicStiff, const Vector &basicForce);
    virtual const Matrix &getInitialGlobalStiffMatrix(const Matrix &basicStiff);
    
    // method used to rotate consistent mass matrix
    const Matrix &getGlobalMatrixFromLocal(const Matrix &local);
    
    // methods used in post-processing only
    const Vector &getPointGlobalCoordFromLocal(const Vector &localCoords);
    const Vector &getPointGlobalDisplFromBasic(double xi, const Vector &basicDisps);
    const Vector &getPointLocalDisplFromBasic(double xi, const Vector &basicDisps);    
    

    // Sensitivity
    
    const Vector & getBasicDisplFixedGrad();
    const Vector & getBasicDisplTotalGrad(int gradNumber);
    const Vector &getGlobalResistingForceShapeSensitivity (const Vector &basicForce, const Vector &p0, int grad);
    bool isShapeSensitivity();
    double getLengthGrad();
    double getd1overLdh();

    
    // MovableObject
    virtual int sendSelf(int cTag, Channel &theChannel);
    virtual int recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    // TaggedObject
    void Print(OPS_Stream &s, int flag = 0);

protected:

private:

    int computeElemtLengthAndOrient();
    void compTransfMatrixLocalGlobal(Matrix &Tlg);
    
    // internal data
    Node *nodeIPtr, *nodeJPtr;  // pointers to the element two endnodes
//  std::array<Node*, 2> nodes;                // pointers to the element two endnodes
    
    double *nodeIOffset, *nodeJOffset;    // rigid joint offsets

    Vector3D xi, xj, vz;

    Matrix3D R;         // rotation matrix
    double RWI[3][3];
    double RWJ[3][3];

    double L;           // undeformed element length

//  static Matrix Tlg;  // matrix that transforms from global to local coordinates
    static Matrix kg;   // global stiffness matrix

    double *nodeIInitialDisp, *nodeJInitialDisp;
    bool initialDispChecked;
};

#endif

