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

typedef std::vector<int> Layout;

enum {
 CRDTR_TAG_CorotFrameTransfWarping3d,
 CRDTR_TAG_CorotFrameTransf3d,
 CRDTR_TAG_LinearFrameTransf3d,
 CRDTR_TAG_PDeltaFrameTransf3d
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
  FrameTransform<nn,ndf>(int tag) : TaggedObject(tag) {}

  // TODO(cmp) : make (almost?) everything pure virtual
  virtual FrameTransform<nn,ndf> *getCopy() =0;

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
  virtual int update() = 0;
  virtual int commit() = 0;
  virtual int revertToLastCommit() = 0;
  virtual int revertToStart() = 0;

  virtual double getInitialLength() = 0;
  virtual double getDeformedLength() = 0;

  virtual VectorND<nn*ndf>    pushResponse(VectorND<nn*ndf>&pl) =0;
  virtual MatrixND<nn*ndf,nn*ndf> pushResponse(MatrixND<nn*ndf,nn*ndf>& kl, const VectorND<nn*ndf>& pl) =0;

  VectorND<nn*ndf>    pushConstant(const VectorND<nn*ndf>&pl);
  MatrixND<nn*ndf,nn*ndf> pushConstant(const MatrixND<nn*ndf,nn*ndf>& kl);

  //
  virtual int getLocalAxes(Vector3D &x, Vector3D &y, Vector3D &z) = 0;

  // Recorders
  virtual Response *setResponse(const char **argv, int argc, 
                                OPS_Stream &theHandler) {
    return nullptr;
  };
  virtual int getResponse(int responseID, Information &eleInformation) {
    return -1;
  };

  // Sensitivity
  virtual const Vector &getBasicDisplTotalGrad(int grad)=0;
  virtual const Vector &getBasicDisplFixedGrad()=0;
  virtual const Vector &getGlobalResistingForceShapeSensitivity(const Vector &pb, const Vector &p0, int grad)=0;
  virtual bool   isShapeSensitivity() {return false;}
  virtual double getLengthGrad() {return 0.0;}
  virtual double getd1overLdh() {return 0.0;}
  //
};

#include "FrameTransform.tpp"

//
// 2D
//
class FrameTransform2d : public CrdTransf {
  public:
  FrameTransform2d(int tag, int classTag) : CrdTransf(tag, classTag) {};
};

//
// 3D
//
class FrameTransform3d : public CrdTransf {
public:
  enum {
    N, Vy, Vz, T, My, Mz,
  };

public:
  FrameTransform3d(int tag, int classTag) : CrdTransf(tag, classTag) {}

  // TODO(cmp) : make almost everything pure virtual
  virtual FrameTransform3d *getCopy() {
    return nullptr;
  }

  virtual CrdTransf *getCopy3d() {
    return getCopy();
  }

  //
  //
  //


  virtual VectorND<12>    pushResponse(VectorND<12>&pl) {
    static VectorND<12> empty{};
    return empty;
  }

  virtual VectorND<12>    pushConstant(const VectorND<12>&pl) const {
    static VectorND<12> empty{};
    return empty;
  }
  virtual MatrixND<12,12> pushResponse(MatrixND<12,12>& kl, const VectorND<12>& pl) {
    static MatrixND<12,12> empty{};
    return empty;
  }
  virtual MatrixND<12,12> pushConstant(const MatrixND<12,12>& kl) {
    static MatrixND<12,12> empty{};
    return empty;
  }

};

template<int nn, int ndf, typename VecT>
static inline VectorND<nn*ndf>
getLocal(VecT ug, const Matrix3D& R, double nodeIOffset[], double nodeJOffset[])
{
  VectorND<nn*ndf> ul;

  for (int i=0; i<4; i++)
    for (int j=0; j<3; j++)
      ul[i*3+j] = R(0,j)*ug[3*i] + R(1,j)*ug[3*i+1] + R(2,j)*ug[3*i+2];

  double Wu[3];
  if (nodeIOffset) {
    Wu[0] =  nodeIOffset[2] * ug[4] - nodeIOffset[1] * ug[5];
    Wu[1] = -nodeIOffset[2] * ug[3] + nodeIOffset[0] * ug[5];
    Wu[2] =  nodeIOffset[1] * ug[3] - nodeIOffset[0] * ug[4];

    ul[0] += R(0,0) * Wu[0] + R(1,0) * Wu[1] + R(2,0) * Wu[2];
    ul[1] += R(0,1) * Wu[0] + R(1,1) * Wu[1] + R(2,1) * Wu[2];
    ul[2] += R(0,2) * Wu[0] + R(1,2) * Wu[1] + R(2,2) * Wu[2];
  }

  if (nodeJOffset) {
    Wu[0] =  nodeJOffset[2] * ug[10] - nodeJOffset[1] * ug[11];
    Wu[1] = -nodeJOffset[2] * ug[ 9] + nodeJOffset[0] * ug[11];
    Wu[2] =  nodeJOffset[1] * ug[ 9] - nodeJOffset[0] * ug[10];

    ul[6] += R(0,0) * Wu[0] + R(1,0) * Wu[1] + R(2,0) * Wu[2];
    ul[7] += R(0,1) * Wu[0] + R(1,1) * Wu[1] + R(2,1) * Wu[2];
    ul[8] += R(0,2) * Wu[0] + R(1,2) * Wu[1] + R(2,2) * Wu[2];
  }

  return ul;
}

template<int nn, int ndf, typename VecT>
static inline VectorND<nn*ndf>
getLocal(VectorND<nn*ndf>& ug, const Matrix3D& R, std::array<Vector3D,nn>* offset)
{
  VectorND<nn*ndf> ul = ug;

  for (int i=0; i<nn; i++)
    for (int j=0; j<6; j++)
      ul[i*ndf+j%3] = R(0,j%3)*ug[i*ndf] + R(1,j%3)*ug[i*ndf+1] + R(2,j%3)*ug[3*i+2];

  if (offset) {
    double Wu[3];
    std::array<Vector3D, nn>& offsets = *offset;

    Wu[0] =  offsets[0][2] * ug[4] - offsets[0][1] * ug[5];
    Wu[1] = -offsets[0][2] * ug[3] + offsets[0][0] * ug[5];
    Wu[2] =  offsets[0][1] * ug[3] - offsets[0][0] * ug[4];

    ul[0] += R(0,0) * Wu[0] + R(1,0) * Wu[1] + R(2,0) * Wu[2];
    ul[1] += R(0,1) * Wu[0] + R(1,1) * Wu[1] + R(2,1) * Wu[2];
    ul[2] += R(0,2) * Wu[0] + R(1,2) * Wu[1] + R(2,2) * Wu[2];

    Wu[0] =  offsets[1][2] * ug[10] - offsets[1][1] * ug[11];
    Wu[1] = -offsets[1][2] * ug[ 9] + offsets[1][0] * ug[11];
    Wu[2] =  offsets[1][1] * ug[ 9] - offsets[1][0] * ug[10];

    ul[6] += R(0,0) * Wu[0] + R(1,0) * Wu[1] + R(2,0) * Wu[2];
    ul[7] += R(0,1) * Wu[0] + R(1,1) * Wu[1] + R(2,1) * Wu[2];
    ul[8] += R(0,2) * Wu[0] + R(1,2) * Wu[1] + R(2,2) * Wu[2];
  }

  return ul;
}

static inline VectorND<6>
getBasic(VectorND<12>& ul, double oneOverL)
{
  VectorND<6> ub;
  // VectorND<12> ul = getLocal<2,6>(ug, R, nodeIOffset, nodeJOffset);

#if 0
  static double ul[12];

  for (int i=0; i<4; i++)
    for (int j=0; j<3; j++)
      ul[i*3+j] = R(0,j)*ug[3*i] + R(1,j)*ug[3*i+1] + R(2,j)*ug[3*i+2];

  double Wu[3];
  if (nodeIOffset) {
    Wu[0] =  nodeIOffset[2] * ug[4] - nodeIOffset[1] * ug[5];
    Wu[1] = -nodeIOffset[2] * ug[3] + nodeIOffset[0] * ug[5];
    Wu[2] =  nodeIOffset[1] * ug[3] - nodeIOffset[0] * ug[4];

    ul[0] += R(0,0) * Wu[0] + R(1,0) * Wu[1] + R(2,0) * Wu[2];
    ul[1] += R(0,1) * Wu[0] + R(1,1) * Wu[1] + R(2,1) * Wu[2];
    ul[2] += R(0,2) * Wu[0] + R(1,2) * Wu[1] + R(2,2) * Wu[2];
  }

  if (nodeJOffset) {
    Wu[0] =  nodeJOffset[2] * ug[10] - nodeJOffset[1] * ug[11];
    Wu[1] = -nodeJOffset[2] * ug[ 9] + nodeJOffset[0] * ug[11];
    Wu[2] =  nodeJOffset[1] * ug[ 9] - nodeJOffset[0] * ug[10];

    ul[6] += R(0,0) * Wu[0] + R(1,0) * Wu[1] + R(2,0) * Wu[2];
    ul[7] += R(0,1) * Wu[0] + R(1,1) * Wu[1] + R(2,1) * Wu[2];
    ul[8] += R(0,2) * Wu[0] + R(1,2) * Wu[1] + R(2,2) * Wu[2];
  }
#endif

  ub[0] = ul[6] - ul[0];

  double tmp;
  tmp   = oneOverL * (ul[1] - ul[7]);
  ub[1] = ul[ 5] + tmp;
  ub[2] = ul[11] + tmp;
  tmp   = oneOverL * (ul[8] - ul[2]);
  ub[3] = ul[ 4] + tmp;
  ub[4] = ul[10] + tmp;
  ub[5] = ul[ 9] - ul[3];

  return ub;
}

template <int nen=2, int ndf=6>
static VectorND<nen*ndf>
pushLocal(const Vector& q, double L)
{
  
  VectorND<12> pl;

  double q0 = q(0);
  double q1 = q(1);
  double q2 = q(2);
  double q3 = q(3);
  double q4 = q(4);
  double q5 = q(5);

  double oneOverL = 1.0 / L;

  pl[0]  = -q0;                    // Ni
  pl[1]  =  oneOverL * (q1 + q2);  // Viy
  pl[2]  = -oneOverL * (q3 + q4);  // Viz
  pl[3]  = -q5;                    // Ti
  pl[4]  =  q3; 
  pl[5]  =  q1;                    // Mzi
  pl[6]  =  q0;                    // Nj
  pl[7]  = -pl[1];                 // Vjy
  pl[8]  = -pl[2];                 // Vjz
  pl[9]  =  q5;                    // Tj
  pl[10] =  q4;
  pl[11] =  q2;                    // Mzj

  return pl;
}


static inline void 
formOffsets(const Matrix3D& R, 
            const double nodeIOffset[3], 
            const double nodeJOffset[3], 
            double RWI[3][3], double RWJ[3][3])
{
  if (nodeIOffset) {
    // Compute RWI
    RWI[0][0] = -R(1,0) * nodeIOffset[2] + R(2,0) * nodeIOffset[1];
    RWI[1][0] = -R(1,1) * nodeIOffset[2] + R(2,1) * nodeIOffset[1];
    RWI[2][0] = -R(1,2) * nodeIOffset[2] + R(2,2) * nodeIOffset[1];

    RWI[0][1] =  R(0,0) * nodeIOffset[2] - R(2,0) * nodeIOffset[0];
    RWI[1][1] =  R(0,1) * nodeIOffset[2] - R(2,1) * nodeIOffset[0];
    RWI[2][1] =  R(0,2) * nodeIOffset[2] - R(2,2) * nodeIOffset[0];

    RWI[0][2] = -R(0,0) * nodeIOffset[1] + R(1,0) * nodeIOffset[0];
    RWI[1][2] = -R(0,1) * nodeIOffset[1] + R(1,1) * nodeIOffset[0];
    RWI[2][2] = -R(0,2) * nodeIOffset[1] + R(1,2) * nodeIOffset[0];
  }

  if (nodeJOffset) {
    // Compute RWJ
    RWJ[0][0] = -R(1,0) * nodeJOffset[2] + R(2,0) * nodeJOffset[1];
    RWJ[1][0] = -R(1,1) * nodeJOffset[2] + R(2,1) * nodeJOffset[1];
    RWJ[2][0] = -R(1,2) * nodeJOffset[2] + R(2,2) * nodeJOffset[1];

    RWJ[0][1] =  R(0,0) * nodeJOffset[2] - R(2,0) * nodeJOffset[0];
    RWJ[1][1] =  R(0,1) * nodeJOffset[2] - R(2,1) * nodeJOffset[0];
    RWJ[2][1] =  R(0,2) * nodeJOffset[2] - R(2,2) * nodeJOffset[0];

    RWJ[0][2] = -R(0,0) * nodeJOffset[1] + R(1,0) * nodeJOffset[0];
    RWJ[1][2] = -R(0,1) * nodeJOffset[1] + R(1,1) * nodeJOffset[0];
    RWJ[2][2] = -R(0,2) * nodeJOffset[1] + R(1,2) * nodeJOffset[0];
  }
}



#endif // include guard
