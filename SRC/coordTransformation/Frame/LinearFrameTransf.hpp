//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the class definition for
// LinearFrameTransf.h. LinearFrameTransf provides the
// abstraction of a linear transformation for a spatial frame
// between the global and basic coordinate systems
//
// Written: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 04/2000
//
#ifndef LinearFrameTransf_hpp
#define LinearFrameTransf_hpp

#include <array>
#include <FrameTransform.h>
#include <Vector.h>
#include <Matrix.h>

template <int nn, int ndf>
class LinearFrameTransf: public FrameTransform<nn,ndf>
{
public:
    constexpr static int n = nn*ndf;

    LinearFrameTransf(int tag, 
                      const Vector3D &vecxz,
                      const std::array<Vector3D, nn> *offset=nullptr,
                      int offset_flags = 0);

    ~LinearFrameTransf();
    
    const char *getClassType() const {return "LinearFrameTransf";}
    
    virtual int getLocalAxes(Vector3D &x, Vector3D &y, Vector3D &z);
    
    virtual FrameTransform<nn,ndf> *getCopy() const;

    virtual double getInitialLength();
    virtual double getDeformedLength();
    
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

    // method used to rotate consistent mass matrix
    const Matrix &getGlobalMatrixFromLocal(const Matrix &local);
    

    // Sensitivity
    
    const Vector & getBasicDisplFixedGrad();
    const Vector & getBasicDisplTotalGrad(int gradNumber);
    const Vector &getGlobalResistingForceShapeSensitivity (const Vector &basicForce, const Vector &p0, int grad);
    bool isShapeSensitivity();
    double getLengthGrad();
    double getd1overLdh();

    // TaggedObject
    void Print(OPS_Stream &s, int flag = 0);

    static inline VectorND<nn*ndf> 
    pullConstant(const VectorND<nn*ndf>& ug, 
                const Matrix3D& R, 
                const std::array<Vector3D, nn> *offset = nullptr);

private:

    int computeElemtLengthAndOrient();

    std::array<Node*, nn> nodes;

    // Rigid joint offsets
    std::array<Vector3D, nn> *offsets;
    int offset_flags;

    Vector3D xi, xj, vz;

    Matrix3D R;         // rotation matrix

    double L;           // undeformed element length

//  static Matrix Tlg;  // matrix that transforms from global to local coordinates
//  static Matrix kg;   // global stiffness matrix

    std::array<VectorND<ndf>*, nn> u_init;
    bool initialDispChecked;
};

#include "LinearFrameTransf.tpp"
#endif

