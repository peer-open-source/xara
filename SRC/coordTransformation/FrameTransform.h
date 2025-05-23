//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Written: cmp
//
#ifndef FrameTransform_h
#define FrameTransform_h

#include <vector>
#include <Versor.h>
#include <VectorND.h>
#include <MatrixND.h>
#include <Matrix3D.h>
#include <CrdTransf.h>

#define MAYBE_STATIC static

using OpenSees::VectorND;
using OpenSees::MatrixND;
using OpenSees::Matrix3D;


enum {
 CRDTR_TAG_CorotFrameTransfWarping3d,
 CRDTR_TAG_CorotFrameTransf3d,
 CRDTR_TAG_LinearFrameTransf3d,
 CRDTR_TAG_PDeltaFrameTransf3d
};

enum {
  OffsetGlobal = 0,
  OffsetLocal  = 1,
  OffsetNormalized = 2
};

//
// Generalized 
//
template <int nn, int ndf>
class FrameTransform : public TaggedObject
{
public:
  constexpr static int ndm = 3;

public:
  FrameTransform(int tag) : TaggedObject(tag) {}

  virtual FrameTransform<nn,ndf> *getCopy() const =0;

  virtual VectorND<nn*ndf> getStateVariation() =0;

  virtual Vector3D  getNodePosition(int tag) =0;
  virtual Vector3D  getNodeRotationLogarithm(int tag) =0;
  // virtual Versor         getNodeRotation(int tag);
  // virtual Vector3D       getNodeRotationVariation(int tag);
  // virtual VectorND<ndf>  getNodeRotationIncrement(int tag);

  // virtual VectorND<ndf>  getNodeLogarithm(int tag) =0;
  // virtual VectorND<ndf>  getNodeVariation(int tag) =0;
  // virtual VectorND<ndf>  getNodeVelocity(int tag);
  // virtual VectorND<ndm>  getNodeLocation(int tag);

  virtual int initialize(std::array<Node*, nn>& nodes)=0;
  virtual int update() =0;
  virtual int commit() =0;
  virtual int revertToLastCommit() =0;
  virtual int revertToStart() =0;

  virtual double getInitialLength() =0;
  virtual double getDeformedLength() =0;
  virtual const std::array<Vector3D,nn> *getRigidOffsets() const =0;

  virtual VectorND<nn*ndf>    pushResponse(VectorND<nn*ndf>&pl) =0;
  virtual MatrixND<nn*ndf,nn*ndf> pushResponse(MatrixND<nn*ndf,nn*ndf>& kl, const VectorND<nn*ndf>& pl) =0;

  // // template<
  // VectorND<nn*ndf-6> pushResponse(VectorND<nn*ndf-6>&pb);

  VectorND<nn*ndf>    pushConstant(const VectorND<nn*ndf>&pl);
  MatrixND<nn*ndf,nn*ndf> pushConstant(const MatrixND<nn*ndf,nn*ndf>& kl);

  //
  virtual int getLocalAxes(Vector3D &x, Vector3D &y, Vector3D &z) const =0;

  // Recorders
  virtual Response *setResponse(const char **argv, int argc, 
                                OPS_Stream &theHandler) {
    return nullptr;
  };
  virtual int getResponse(int responseID, Information &) {
    return -1;
  };

  // Sensitivity
  // virtual const Vector &getBasicDisplTotalGrad(int grad)=0;
  // virtual const Vector &getBasicDisplFixedGrad()=0;
  // virtual const Vector &getGlobalResistingForceShapeSensitivity(const Vector &pb, const Vector &p0, int grad)=0;
  virtual bool   isShapeSensitivity() {return false;}
  virtual double getLengthGrad() {return 0.0;}
  virtual double getd1overLdh() {return 0.0;}

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
};

#include "FrameTransform.tpp"

//
// 3D
//
class FrameTransform3d : public CrdTransf {

public:
  FrameTransform3d(int tag, int classTag) : CrdTransf(tag, classTag) {}

  virtual FrameTransform3d *getCopy() {
    return nullptr;
  }

  virtual CrdTransf *getCopy3d() {
    return getCopy();
  }

};

#endif // include guard
