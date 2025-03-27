//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
// 
// The purpose of this class is to wrap the more general FrameTransform<>
// templates to reproduce the legacy CrdTransf classes that were derived
// for elements in the "basic" coordinate system.
//
// cmp
//
#ifndef BasicFrameTransf3d_h
#define BasicFrameTransf3d_h

#include <array>
#include <FrameTransform.h>
#include <Vector.h>
#include <Matrix.h>

namespace OpenSees {
class BasicFrameTransf3d: public FrameTransform3d
{
public:
    BasicFrameTransf3d(FrameTransform<2,6> *t);

    ~BasicFrameTransf3d();

    virtual int getLocalAxes(Vector &x, Vector &y, Vector &z);

    virtual FrameTransform3d *getCopy();

    virtual double getInitialLength();
    virtual double getDeformedLength();

    virtual int initialize(Node *ni, Node *nj) override final;
    virtual int update() override final;
    virtual int commitState() override final;
    virtual int revertToLastCommit() override final;
    virtual int revertToStart() override final;

    virtual const Vector &getBasicTrialDisp() override final;
    virtual const Vector &getBasicIncrDisp() override final;
    virtual const Vector &getBasicIncrDeltaDisp() override final;
    virtual const Vector &getBasicTrialVel() override final;

    virtual VectorND<12>    pushResponse(VectorND<12>&pl) override final;
    virtual VectorND<12>    pushConstant(const VectorND<12>&pl) const override final;

    virtual MatrixND<12,12> pushResponse(MatrixND<12,12>& kl, const VectorND<12>& pl) override final;
    virtual MatrixND<12,12> pushConstant(const MatrixND<12,12>& kl) override final;

    virtual const Vector &getGlobalResistingForce(const Vector &basicForce, const Vector &p0) final;
    virtual const Matrix &getGlobalStiffMatrix(const Matrix &basicStiff, const Vector &basicForce) final;
    virtual const Matrix &getInitialGlobalStiffMatrix(const Matrix &basicStiff) final;

    // method used to rotate consistent mass matrix
    const Matrix &getGlobalMatrixFromLocal(const Matrix &local);
    
    // methods used in post-processing only
    const Vector &getPointGlobalCoordFromLocal(const Vector &localCoords);
    const Vector &getPointGlobalDisplFromBasic(double xi, const Vector &basicDisps);
    const Vector &getPointLocalDisplFromBasic( double xi, const Vector &basicDisps);    

    //
    // Sensitivity
    //
    const Vector & getBasicDisplFixedGrad();
    const Vector & getBasicDisplTotalGrad(int gradNumber);
    const Vector &getGlobalResistingForceShapeSensitivity (const Vector &basicForce, const Vector &p0, int grad);
    bool isShapeSensitivity();
    double getLengthGrad();
    double getd1overLdh();


    // MovableObject
    virtual int sendSelf(int tag, Channel &);
    virtual int recvSelf(int tag, Channel &, FEM_ObjectBroker &);
    const char *getClassType() const {return "BasicFrameTransf3d";}
    
    // TaggedObject
    void Print(OPS_Stream &s, int flag = 0);

protected:
private:
    FrameTransform<2,6> &t;
};
} // namespace OpenSees
#endif

