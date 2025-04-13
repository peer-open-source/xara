//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the class definition for
// RigidFrameTransf.h. RigidFrameTransf provides the
// abstraction of a linear transformation for a spatial frame
// between the global and basic coordinate systems
//
// Written: cmp
// Created: 04/2025
//
#ifndef RigidFrameTransf_hpp
#define RigidFrameTransf_hpp

#include <array>
#include <FrameTransform.h>
#include <Vector3D.h>
#include <MatrixND.h>

template <int nn, int ndf, typename BasisT>
class RigidFrameTransf: public FrameTransform<nn,ndf>
{
public:
    constexpr static int n = nn*ndf;

    RigidFrameTransf(int tag, 
                      const Vector3D &vecxz,
                      const std::array<Vector3D, nn> *offset=nullptr,
                      int offset_flags = 0);

    ~RigidFrameTransf();
    
    const char *getClassType() const {return "RigidFrameTransf";}
    
    virtual int getLocalAxes(Vector3D &x, Vector3D &y, Vector3D &z) const;
    
    virtual FrameTransform<nn,ndf> *getCopy() const;

    virtual double getInitialLength();
    virtual double getDeformedLength();
    virtual const std::array<Vector3D,nn> *getRigidOffsets() const {return offsets;}
    
    virtual int initialize(std::array<Node*, nn>& new_nodes) override final;
    virtual int update() override final;
    virtual int commit() override final;
    virtual int revertToLastCommit() override final;
    virtual int revertToStart() override final;

    virtual VectorND<nn*ndf> getStateVariation() final;
    virtual Vector3D getNodePosition(int tag) final;
    virtual Vector3D getNodeRotationLogarithm(int tag) final;

    virtual VectorND<nn*ndf>        pushResponse(VectorND<nn*ndf>&pl) override final;
    virtual MatrixND<nn*ndf,nn*ndf> pushResponse(MatrixND<nn*ndf,nn*ndf>& kl, const VectorND<nn*ndf>& pl) override final;

    // // method used to rotate consistent mass matrix
    // const Matrix &getGlobalMatrixFromLocal(const Matrix &local);
    

    // Sensitivity
    //
    // const Vector & getBasicDisplFixedGrad();
    // const Vector & getBasicDisplTotalGrad(int gradNumber);
    // const Vector &getGlobalResistingForceShapeSensitivity (const Vector &basicForce, const Vector &p0, int grad);
    bool isShapeSensitivity();
    double getLengthGrad();
    double getd1overLdh();

    // TaggedObject
    void Print(OPS_Stream &s, int flag = 0);

private:

    int computeElemtLengthAndOrient();

    inline VectorND<nn*ndf> 
    pullConstant(const VectorND<nn*ndf>& ug, 
                const Matrix3D& R, 
                const std::array<Vector3D, nn> *offset = nullptr,
                int offset_flags = 0);

    template<const Vector& (Node::*Getter)()>
    const Vector3D
    pullPosition(int node)
    {
        const Vector &u = (nodes[node]->*Getter)();

        Vector3D v;
        for (int i=0; i<3; i++)
          v[i] = u[i];

        // 1) Offsets
        if (offsets) {
          if (!(offset_flags&OffsetLocal)) {
            Vector3D w {u[3], u[4], u[5]};
            v -= offsets->at(node).cross(w);
          }
        }

        // 2) Constant Rotation
        return R^v;
    }

    std::array<Node*, nn> nodes;
    std::array<Vector3D, nn> ur; // rotation vector
    std::array<Vector3D, nn> ux; // displacement vector

    std::array<Vector3D, nn> *offsets;
    int offset_flags;

    Vector3D xi, xj, vz;
    Matrix3D R0;        // rotation matrix
    double L;           // undeformed element length

    BasisT basis;
};

#include "RigidFrameTransf.tpp"
#endif

