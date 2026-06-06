//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Please cite the following resource in any derivative works:
//
// [1] Perez, C.M., and Filippou F.C.. "On Nonlinear Geometric Transformations
//     of Finite Elements" Int. J. Numer. Meth. Engrg. 2024; 
//     https://doi.org/10.1002/nme.7506
//
//===----------------------------------------------------------------------===//
//
// Written: cmp
//
#ifndef FrameTransform_h
#define FrameTransform_h

#include <VectorND.h>
#include <MatrixND.h>
#include <Matrix3D.h>
#include <Vector3D.h>
#include <TaggedObject.h>

class Information;
class Response;
class Node;

enum {
 CRDTR_TAG_CorotFrameTransfWarping3d,
 CRDTR_TAG_CorotFrameTransf3d,
 CRDTR_TAG_LinearFrameTransf3d,
 CRDTR_TAG_PDeltaFrameTransf3d
};

enum {
  OffsetGlobal     = 0, // 1<<0,
  OffsetLocal      = 1, // 1<<1,
  OffsetNormalized = 2, // 1<<2,

  LogIter          = 1<<3,
  LogIncr          = 1<<4,
  LogInit          = 1<<5,
  LogDefault       = 1<<6
};

namespace OpenSees {

namespace Transform {
enum Action : int {
  None        = 0,

  Logarithm   = 1u << 0,
  // LocalOffset = 1u << 1,
  Isometry    = 1u << 2,
  Rotation    = 1u << 3,
  Offset      = 1u << 4,
  Bubnov      = 1u << 6,
  Adjoint     = 1u << 7,
  Tangent     = 1u << 8,
  Total  = (1u<<0)
         | (1u<<1)
         | (1u<<2) | (1u<<3) | (1u<<4) | (1u<<5) 
         | (1u<<6)| (1u<<7)
         | (1u<<8)
};
}

template <int nn, int ndf>
class FrameTransform : public TaggedObject
{
public:
  explicit FrameTransform(int tag) : TaggedObject(tag) {}

  virtual FrameTransform<nn,ndf> *getCopy() const =0;


  virtual Vector3D  getNodePosition(int tag) =0;
  virtual Vector3D  getNodeRotationLogarithm(int tag) =0;
  virtual Vector3D getNodeRotationUpdateLogarithm(int tag) {
    return Vector3D{};
  }
#if 0
  virtual Versor         getNodeRotation(int tag);
  virtual Vector3D       getNodeRotationVariation(int tag);
  virtual VectorND<ndf>  getNodeRotationIncrement(int tag);

  virtual VectorND<ndf>  getNodeLogarithm(int tag) =0;
  virtual VectorND<ndf>  getNodeVariation(int tag) =0;
  virtual VectorND<ndf>  getNodeVelocity(int tag);
#endif
  // \bar{x}
  virtual Vector3D  getNodeLocation(int tag) {return Vector3D{};}

  virtual int initialize(std::array<Node*, nn>& nodes)=0;
  virtual int update() noexcept =0;
  virtual int commit() =0;
  virtual int restart() {return 0;} // TODO
  virtual int recover() {return 0;} // TODO

  virtual int revertToLastCommit() =0;
  virtual int revertToStart() =0;

//virtual VectorND<nn*ndf> getStateLogarithm() =0; //
  virtual VectorND<nn*ndf> getStateVariation() =0; // pull

  virtual int push(VectorND<nn*ndf>&pl, int=Transform::Total) =0;
  virtual int push(MatrixND<nn*ndf,nn*ndf>& kl, 
                   const VectorND<nn*ndf>& pl, 
                   int=Transform::Total) =0;

  virtual VectorND<nn*ndf>    pushResponse(VectorND<nn*ndf>&pl) final {
    VectorND<nn*ndf> pg{pl};
    push(pg, Transform::Total);
    return pg;
  }

  virtual MatrixND<nn*ndf,nn*ndf> pushResponse(MatrixND<nn*ndf,nn*ndf>& kl,
                                               const VectorND<nn*ndf>& pl) final {
    MatrixND<nn*ndf,nn*ndf> kg{kl};
    push(kg, pl, Transform::Total);
    return kg;              
  }

  virtual double getInitialLength() =0;
  virtual double getDeformedLength() =0;
  Matrix3D getInitialRotation() const {
    Vector3D x, y, z;
    getLocalAxes(x, y, z);
    
    Matrix3D R;
    for (int i=0; i<3; i++) {
      R(i,0) = x[i];
      R(i,1) = y[i];
      R(i,2) = z[i];
    }
    return R;
  }

  virtual Matrix3D getRotation() const noexcept =0;

  virtual MatrixND<3,nn*ndf> 
  getRotationTangent() {
    MatrixND<3,nn*ndf> dR{};
    return dR;
  }
  virtual const std::array<Vector3D,nn> *getRigidOffsets() const =0;

  //
  virtual int getLocalAxes(Vector3D &x, Vector3D &y, Vector3D &z) const =0;

  Vector3D getNormalVector() const {
    Vector3D x, y, z;
    if (getLocalAxes(x, y, z) < 0)
      return Vector3D{{0.0, 0.0, 1.0}};
    return z;
  }

  // Recorders
  virtual Response *setResponse(const char **argv, int argc, OPS_Stream &) {
    return nullptr;
  }
  virtual int getResponse(int responseID, Information &) {
    return -1;
  }

  // Sensitivity

  virtual void   pushGrad(VectorND<nn*ndf>& dp, VectorND<nn*ndf>& pl) {}
  virtual void   pullFixedGrad(VectorND<nn*ndf>&) {}
  virtual void   pullTotalGrad(VectorND<nn*ndf>&, int) {}
  virtual bool   isShapeSensitivity() {return false;}
  virtual double getLengthGrad() {return 0.0;}
  virtual double getd1overLdh() {return 0.0;}

  virtual Matrix3D getRotationSensitivity() {
    return Matrix3D{};
  }

protected:
  constexpr static int ndm = 3;
  static inline constexpr void
  pushRotation(MatrixND<nn*ndf,nn*ndf>& Kg, const Matrix3D& R);

  static inline constexpr void
  pushOffsets(MatrixND<nn*ndf,nn*ndf>& Kg, const std::array<Vector3D,nn>& offsets);

  static int
  Orient(const Vector3D& dx, const Vector3D& vz, Matrix3D &R) {

    // calculate the element local x axis components wrt to the global coordinates

    Vector3D e1 = dx/dx.norm();

    //
    Vector3D e2 = vz.cross(e1);

    const double ynorm = e2.norm();

    if (ynorm == 0.0)
        return -1;

    e2 /= ynorm;

    Vector3D e3 = e1.cross(e2);

    for (int i = 0; i < 3; i++) {
      R(i,0) = e1[i];
      R(i,1) = e2[i];
      R(i,2) = e3[i];
    }
    return 0;
  }
#if 0
  VectorND<nn*ndf>    pushConstant(const VectorND<nn*ndf>&pl);
  MatrixND<nn*ndf,nn*ndf> pushConstant(const MatrixND<nn*ndf,nn*ndf>& kl);
#endif

  void pushRotationOffset(VectorND<nn*ndf>&pl, const Matrix3D& R);
};

}
#include "FrameTransform.tpp"

#endif // include guard
