//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the implementation for the
// RigidFrameTransf class. RigidFrameTransf is a nonlinear
// transformation for a space frame between the global
// and basic coordinate systems
//
// Written: cmp
// Created: 04/2025
//
#pragma once
#include <Vector.h>
#include <Matrix.h>
#include <Matrix3D.h>
#include <Node.h>
#include <Channel.h>
#include <Logging.h>
#include <Rotations.hpp>
#include <RigidFrameTransf.hpp>

using namespace OpenSees;


template <int nn, int ndf, typename BasisT>
RigidFrameTransf<nn,ndf,BasisT>::RigidFrameTransf(int tag, 
                                           const Vector3D &vecxz, 
                                           const std::array<Vector3D, nn> *offset,
                                           int offset_flags)
  : FrameTransform<nn,ndf>(tag),
    Du{0},
    L(0),
    offsets{nullptr},
    offset_flags(offset_flags),
    basis{nodes}
{
  R.zero();

  for (int i=0; i<3; i++)
    vz[i] = vecxz[i];

  R(0,2) = vz(0);
  R(1,2) = vz(1);
  R(2,2) = vz(2);

  // Rigid joint offsets
  if (offset != nullptr) {
    offsets = new std::array<Vector3D, nn>{};
    *offsets = *offset;
  }
}



template <int nn, int ndf, typename BasisT>
RigidFrameTransf<nn,ndf,BasisT>::~RigidFrameTransf()
{
  if (offsets != nullptr)
    delete offsets;
}

template <int nn, int ndf, typename BasisT>
int
RigidFrameTransf<nn,ndf,BasisT>::commit()
{
  return 0;
}

template <int nn, int ndf, typename BasisT>
int
RigidFrameTransf<nn,ndf,BasisT>::revertToLastCommit()
{
  return 0;
}

template <int nn, int ndf, typename BasisT>
int
RigidFrameTransf<nn,ndf,BasisT>::revertToStart()
{
  return 0;
}


template <int nn, int ndf, typename BasisT>
int
RigidFrameTransf<nn,ndf,BasisT>::initialize(std::array<Node*, nn>& new_nodes)
{

  for (int i=0; i<nn; i++) {
    nodes[i] = new_nodes[i];
    if (nodes[i] == nullptr) {
      opserr << "invalid pointers to the element nodes\n";
      return -1;
    }
  }

  int error;
  // get element length and orientation
  if ((error = this->computeElemtLengthAndOrient()))
    return error;

  return 0;
}


template <int nn, int ndf, typename BasisT>
int
RigidFrameTransf<nn,ndf,BasisT>::computeElemtLengthAndOrient()
{

  const Vector &XI = nodes[   0]->getCrds();
  const Vector &XJ = nodes[nn-1]->getCrds();

  for (int i=0; i<3; i++) {
    xi[i] = XI[i];
    xj[i] = XJ[i];
  }
  
  Vector3D dx = xj - xi;

  if (offsets != nullptr) {
    for (int i=0; i<3; i++)
      dx(i) -= (*offsets)[   0][i];
    for (int i=0; i<3; i++)
      dx(i) += (*offsets)[nn-1][i];
  }


  if (u_init[0] != 0) {
    for (int i=0; i<3; i++)
      dx(i) -= (*u_init[0])[i];
  }

  if (u_init[nn-1] != 0) {
    for (int i=0; i<3; i++)
      dx(i) += (*u_init[nn-1])[i];
  }

  // calculate the element length
  L = dx.norm();

  if (L == 0.0)
    return -2;

  return FrameTransform<nn,ndf>::Orient(dx, vz, R);
}


template <int nn, int ndf, typename BasisT>
int
RigidFrameTransf<nn,ndf,BasisT>::getLocalAxes(Vector3D &e1, Vector3D &e2, Vector3D &e3) const
{
  for (int i = 0; i < 3; i++) {
    e1[i] = R(i,0);
    e2[i] = R(i,1);
    e3[i] = R(i,2);
  }
  return 0;
}

template <int nn, int ndf, typename BasisT>
double
RigidFrameTransf<nn,ndf,BasisT>::getInitialLength()
{
  return L;
}

template <int nn, int ndf, typename BasisT>
double
RigidFrameTransf<nn,ndf,BasisT>::getDeformedLength()
{
  return L;
}


//
// Pull
//
template <int nn, int ndf, typename BasisT>
int
RigidFrameTransf<nn,ndf,BasisT>::update()
{
  basis.update();

  Versor   R = basis.getRotation();
  Vector3D c = basis.getTranslation();


  for (int i=0; i<nn; i++)
    ur[i] = LogSO3(R^nodes[i]->getRotation())


  return 0;
}

template <int nn, int ndf, typename BasisT>
VectorND<nn*ndf> 
RigidFrameTransf<nn,ndf,BasisT>::pullConstant(const VectorND<nn*ndf>& ug, 
             const Matrix3D& R, 
             const std::array<Vector3D, nn> *offset,
             int offset_flags) 
{

  constexpr static int N = nn * ndf;

  // Initialize ul = ug
  VectorND<N> ul = ug;

  // (1)
  // Do ui -= ri x wi
  if constexpr (ndf >= 6)
    if (offset && !(offset_flags&OffsetLocal)) {
      const std::array<Vector3D, nn>& offsets = *offset;
      for (int i=0; i<nn; i++) {

        const int j = i * ndf;
        Vector3D w {ul[j+3],ul[j+4],ul[j+5]};

        ul.assemble(j, offsets[i].cross(w), -1.0);
      }
    }

  // (2) Rotations and translations
  for (int i=0; i<nn; i++) {
    const int j = i * ndf;
    ul.insert(j  , R^Vector3D{ul[j+0], ul[j+1], ul[j+2]}, 1.0);
    ul.insert(j+3, R^Vector3D{ul[j+3], ul[j+4], ul[j+5]}, 1.0);
  }

  // 3)
  if constexpr (ndf >= 6)
    if (offset && (offset_flags&OffsetLocal)) {
      const std::array<Vector3D, nn>& offsets = *offset;
      for (int i=0; i<nn; i++) {

        const int j = i * ndf;
        Vector3D w {ul[j+3],ul[j+4],ul[j+5]};

        ul.assemble(j, offsets[i].cross(w), -1.0);
      }
    }

  // (4)
  // TODO (nn>2)
  constexpr static Vector3D iv {1,0,0};
  Vector3D uI = ul.template extract<3>(0);
  Vector3D Du = ul.template extract<3>((nn-1)*ndf) - uI;
  Vector3D ixDu = iv.cross(Du);
  for (int i=0; i<nn; i++) {
    // Translation
    ul.assemble(i*ndf, uI, -1.0);
    for (int j=1; j<3; j++)
      ul[i*ndf+j] -= double(i)/(nn-1.0)*Du[j];

    // Rotation
    ul.assemble(i*ndf+3, ixDu, -1.0/L);
  }

  return ul;
}

template <int nn, int ndf, typename BasisT>
VectorND<nn*ndf>
RigidFrameTransf<nn,ndf,BasisT>::getStateVariation()
{

  static VectorND<nn*ndf> ug;
  for (int i=0; i<nn; i++) {
    const Vector &ddu = nodes[i]->getIncrDeltaDisp();
    for (int j = 0; j < ndf; j++) {
      ug[i*ndf+j] = ddu(j);
    }
  }
  return RigidFrameTransf<nn,ndf,BasisT>::pullConstant(ug, R, offsets, offset_flags);
}

template <int nn, int ndf, typename BasisT>
Vector3D
RigidFrameTransf<nn,ndf,BasisT>::getNodePosition(int node)
{

  Vector3D v = this->pullPosition<&Node::getTrialDisp>(node) 
             - this->pullPosition<&Node::getTrialDisp>(0);
  return v;
}


template <int nn, int ndf, typename BasisT>
Vector3D
RigidFrameTransf<nn,ndf,BasisT>::getNodeRotationLogarithm(int node)
{
  return ur[node];
}


//
// Push
//
template <int nn, int ndf, typename BasisT>
VectorND<nn*ndf>
RigidFrameTransf<nn,ndf,BasisT>::pushResponse(VectorND<nn*ndf>&p)
{
  VectorND<nn*ndf> pa = p;
  constexpr Vector3D iv{1, 0, 0};
  constexpr Matrix3D ix = Hat(iv);

  // 1) Sum of moments: m = sum_i mi + sum_i (xi x ni)
  Vector3D m{};
  for (int i=0; i<nn; i++) {
    // m += mi
    for (int j=0; j<3; j++)
      m[j] += p[i*ndf+3+j];

    const Vector3D n = Vector3D{p[i*ndf+0], p[i*ndf+1], p[i*ndf+2]};
    m.addVector(1, iv.cross(n), double(i)/double(nn-1)*L);
  }
  const Vector3D ixm = ix*m;

  // 2) Adjust force part
  for (int i=0; i<nn; i++) {
    pa.assemble(i*ndf,  ixm,  (i? 1.0:-1.0)/L);
    pa[i*ndf+3] += m[0]*(i? -1:1)*0.5;
  }

  // 3) Rotate and do joint offsets
  auto pg = this->FrameTransform<nn,ndf>::pushConstant(pa);
  return pg;
}

template <int nn, int ndf, typename BasisT>
MatrixND<nn*ndf,nn*ndf>
RigidFrameTransf<nn,ndf,BasisT>::pushResponse(MatrixND<nn*ndf,nn*ndf>&kb, const VectorND<nn*ndf>&)
{

  MatrixND<nn*ndf,nn*ndf> A{};
  A.addDiagonal(1.0);
  constexpr Vector3D axis{1, 0, 0};
  constexpr Matrix3D ix = Hat(axis);
  constexpr Matrix3D ioi = axis.bun(axis);

  MatrixND<3,ndf> Gb{};
  Gb.template insert<0, 3>(ioi, 0.5);
  for (int a = 0; a<nn; a++) {
    for (int b = 0; b<nn; b++) {
      // TODO(nn>2): Interpolate coordinate?
      if (b == 0)
        Gb.template insert<0,0>(ix, -1/L);
      else if (b == nn-1)
        Gb.template insert<0,0>(ix,  1/L);
      // TODO(nn>2): Interpolate coordinate?
      A.assemble(ix*Gb, a*ndf  , b*ndf,  double(a)/double(nn-1)*L);
      A.assemble(   Gb, a*ndf+3, b*ndf, -1.0);
    }
  }

  MatrixND<12,12> kl;
  kl.addMatrixTripleProduct(0, A, kb, 1);
  return this->FrameTransform<nn,ndf>::pushConstant(kl);
}


template <int nn, int ndf, typename BasisT>
FrameTransform<nn,ndf> *
RigidFrameTransf<nn,ndf,BasisT>::getCopy() const
{

  Vector3D xz;
  xz(0) = R(0,2);
  xz(1) = R(1,2);
  xz(2) = R(2,2);


  RigidFrameTransf *theCopy = new RigidFrameTransf<nn,ndf>(this->getTag(), xz, offsets);

  theCopy->nodes = nodes;
  theCopy->L     = L;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      theCopy->R(j,i) = R(j,i);

  return theCopy;
}


//
// Sensitivity
//
template <int nn, int ndf, typename BasisT>
bool
RigidFrameTransf<nn,ndf,BasisT>::isShapeSensitivity()
{
  int nodeParameterI = nodes[   0]->getCrdsSensitivity();
  int nodeParameterJ = nodes[nn-1]->getCrdsSensitivity();
  // TODO(sensitivity): implement dvz

  return (nodeParameterI != 0 || nodeParameterJ != 0);
}


template <int nn, int ndf, typename BasisT>
double
RigidFrameTransf<nn,ndf,BasisT>::getLengthGrad()
{
  const int di = nodes[0]->getCrdsSensitivity();
  const int dj = nodes[1]->getCrdsSensitivity();

  Vector3D dxi{0.0};
  Vector3D dxj{0.0};

  if (di != 0)
    dxi(di-1) = 1.0;
  if (dj != 0)
    dxj(dj-1) = 1.0;

  return 1/L*(xj - xi).dot(dxj - dxi);
}

template <int nn, int ndf, typename BasisT>
double
RigidFrameTransf<nn,ndf,BasisT>::getd1overLdh()
{
  return -getLengthGrad()/(L*L);
}


template <int nn, int ndf, typename BasisT>
void
RigidFrameTransf<nn,ndf,BasisT>::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"RigidFrameTransf\"";
    s << ", \"vecxz\": [" 
      << R(0,2) << ", " 
      << R(1,2) << ", "
      << R(2,2) << "]";
    if (offsets != nullptr) {
      s << ", \"offsets\": [";
      for (int i=0; i<nn; i++) {
        s << "["
          << (*offsets)[i][0] << ", " 
          << (*offsets)[i][1] << ", "
          << (*offsets)[i][2] << "]";
        if (i < nn-1)
          s << ", ";
      }
      s << "]";
    }

    s << "}";

    return;
  }

  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "\nFrameTransform: " << this->getTag() << " Type: RigidFrameTransf\n";
    s << "\tOrientation: " << Matrix(&R(0,0), 3,3) << "\n";
  }
}

