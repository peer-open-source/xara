//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//----------------------------------------------------------------------------//
//
// Please cite the following resource in any derivative works:
//
// [1] Perez, C.M., and Filippou F.C.. "On Nonlinear Geometric Transformations
//     of Finite Elements" Int. J. Numer. Meth. Engrg. 2024; 
//     https://doi.org/10.1002/nme.7506
//
//===----------------------------------------------------------------------===//

//
// Written: Claudio M. Perez
// Created: 04/2025
//

#pragma once
#include <Vector.h>
#include <Matrix.h>
#include <Matrix3D.h>
#include <Vector3D.h>
#include <Node.h>
#include <Logging.h>
#include <AxisAngle.h>
#include <GroupSO3.h>
#include "EuclidFrameTransf.h"
#include <utility/Unroll.h>

namespace OpenSees {

template <int nn, int ndf, typename IsoT>
EuclidFrameTransf<nn,ndf,IsoT>::EuclidFrameTransf(int tag, 
                                           const Vector3D &vecxz, 
                                           const std::array<Vector3D, nn> *offset,
                                           int offset_flags)
  : FrameTransform<nn,ndf>(tag),
    L(0),
    nodes{},
    ur{},
    offsets{nullptr},
    offset_flags(offset_flags),
    basis{vecxz}
{
  R0.zero();
  R0.addDiagonal(1.0);
  double nz = vecxz.norm();
  for (int i=0; i<3; i++)
    vz[i] = vecxz[i]/nz;

  // Rigid joint offsets
  if (offset != nullptr) {
    offsets = new std::array<Vector3D, nn>{};
    *offsets = *offset;
    basis.setOffsets(offsets);
  }
}



template <int nn, int ndf, typename IsoT>
EuclidFrameTransf<nn,ndf,IsoT>::~EuclidFrameTransf()
{
  if (offsets != nullptr)
    delete offsets;
}

template <int nn, int ndf, typename IsoT>
int
EuclidFrameTransf<nn,ndf,IsoT>::commit()
{
  return 0;
}

template <int nn, int ndf, typename IsoT>
int
EuclidFrameTransf<nn,ndf,IsoT>::revertToLastCommit()
{
  return 0;
}

template <int nn, int ndf, typename IsoT>
int
EuclidFrameTransf<nn,ndf,IsoT>::revertToStart()
{
  return 0;
}


template <int nn, int ndf, typename IsoT>
int
EuclidFrameTransf<nn,ndf,IsoT>::initialize(std::array<Node*, nn>& new_nodes)
{
  for (int i=0; i<nn; i++) {
    nodes[i] = new_nodes[i];
    if (nodes[i] == nullptr) {
      opserr << "invalid pointers to element nodes\n";
      return -1;
    }
    // ensure the node is initialized
    nodes[i]->getTrialRotation();

    ur[i] = AxisAngle(Eye3);
  }
  
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

  // calculate the element length
  L = dx.norm();

  if (L == 0.0)
    return -2;

  int error = basis.initialize(nodes);

  R0 = basis.getRotation();
  // R0 = basis.getInitialRotation();
  return error;
}


template <int nn, int ndf, typename IsoT>
FrameTransform<nn,ndf> *
EuclidFrameTransf<nn,ndf,IsoT>::getCopy() const
{
  return new EuclidFrameTransf<nn,ndf,IsoT>(this->getTag(), vz, offsets);
}


template <int nn, int ndf, typename IsoT>
int
EuclidFrameTransf<nn,ndf,IsoT>::getLocalAxes(Vector3D &e1, Vector3D &e2, Vector3D &e3) const
{
  const Matrix3D R = basis.getInitialRotation();
  for (int i = 0; i < 3; i++) {
    e1[i] = R(i,0);
    e2[i] = R(i,1);
    e3[i] = R(i,2);
  }
  return 0;
}

template <int nn, int ndf, typename IsoT>
double
EuclidFrameTransf<nn,ndf,IsoT>::getInitialLength()
{
  return L;
}


template <int nn, int ndf, typename IsoT>
double
EuclidFrameTransf<nn,ndf,IsoT>::getDeformedLength()
{
  return basis.getLength();
}


//
// Pull
//
template <int nn, int ndf, typename IsoT>
int
EuclidFrameTransf<nn,ndf,IsoT>::update() noexcept
{
  if ((offset_flags & LogIter) && !getenv("XARA_LOG_ITER")) {
    skip_log_iter = true;
  }

  if (basis.update(nodes) < 0) 
    return -1;

  const Matrix3D R = basis.getRotation();

  for (int i=0; i<nn; i++) {
    Versor q = nodes[i]->getTrialRotation();
    const Matrix3D Ri = R^(MatrixFromVersor(q)*R0);
    ddur[i] = AxisAngle(Ri*ur[i].matrix().transpose());
    ur[i] = AxisAngle(Ri);
    // if (offset_flags & LogIter) 
    //   ur[i] = ddur[i];
  }

  return 0;
}

template <int nn, int ndf, typename IsoT>
Versor
EuclidFrameTransf<nn,ndf,IsoT>::getNodeRotation(int tag)
{
  return nodes[tag]->getTrialRotation();
}


template <int nn, int ndf, typename IsoT>
Vector3D 
EuclidFrameTransf<nn,ndf,IsoT>::getNodeLocation(int node)
{
  const Vector &Xn = nodes[node]->getCrds();
  Vector3D X3 {Xn[0], Xn[1], Xn[2]};
  if (offsets != nullptr) {
    X3.addVector(1.0, (*offsets)[node], 1.0);
  }
  Vector3D xn = basis.getRotation()^X3;

  xn += this->pullPosition<&Node::getTrialDisp>(node);
  
  return xn - basis.getLocation();
}


template <int nn, int ndf, typename IsoT>
Vector3D
EuclidFrameTransf<nn,ndf,IsoT>::getNodePosition(int node)
{
  Vector3D u = this->pullPosition<&Node::getTrialDisp>(node);
  u -= basis.getPosition();
  const Vector& Xn = nodes[node]->getCrds();
  Vector3D X3 {Xn[0], Xn[1], Xn[2]};
  if (offsets != nullptr) {
    X3.addVector(1.0, (*offsets)[node], 1.0);
  }
  u += basis.getRotationDelta()^X3;
  return u;
}


template <int nn, int ndf, typename IsoT>
Vector3D
EuclidFrameTransf<nn,ndf,IsoT>::getNodeRotationLogarithm(int node)
{
  return ur[node].vector;
}


template <int nn, int ndf, typename IsoT>
VectorND<nn*ndf>
EuclidFrameTransf<nn,ndf,IsoT>::getStateVariation()
{

  static VectorND<nn*ndf> ul;
  for (int i=0; i<nn; i++) {
    const Vector &ddu = nodes[i]->getIncrDeltaDisp();
    for (int j = 0; j < ndf; j++) {
      ul[i*ndf+j] = ddu(j);
    }
  }


  // -) Global Offsets
  // Do ui -= ri x wi
  if constexpr (ndf >= 6)
    if (offsets && !(offset_flags&OffsetLocal)) [[unlikely]] {
      const std::array<Vector3D, nn> roffsets = this->getCurrentOffsets();
      for (int i=0; i<nn; i++) {

        const int j = i * ndf;
        Vector3D w {ul[j+3],ul[j+4],ul[j+5]};

        ul.assemble(j, roffsets[i].cross(w), -1.0);
      }
    }

  // 2) Isometry
  // 2.1) Rotation
  const Matrix3D& R = basis.getRotation();
  for (int i=0; i<nn; i++) {
    const int j = i * ndf;
    ul.insert(j  , R^Vector3D{ul[j+0], ul[j+1], ul[j+2]}, 1.0);
    ul.insert(j+3, R^Vector3D{ul[j+3], ul[j+4], ul[j+5]}, 1.0);
  }

  // 2.2) Projection
  {
    const Vector3D wr = basis.getRotationVariation(ndf, &ul[0]);
    const Vector3D dc = basis.getPositionVariation(ndf, &ul[0]);

    for (int i=0; i<nn; i++) {
      Vector3D ui = this->getNodePosition(i);
      ul.assemble(i*ndf+0, dc, -1.0);
      ul.assemble(i*ndf+0, ui.cross(wr), 1.0);
      ul.assemble(i*ndf+3, wr, -1.0);
    }
  }

  // -) Offsets
  if constexpr (ndf >= 6)
    if (offsets && (offset_flags&OffsetLocal)) [[unlikely]] {
      const std::array<Vector3D, nn>& offset = *offsets;
      for (int i=0; i<nn; i++) {

        const int j = i * ndf;
        Vector3D w {ul[j+3],ul[j+4],ul[j+5]};

        ul.assemble(j, offset[i].cross(w), -1.0);
      }
    }

  // 3) Element parameters
  if (!(offset_flags & LogIter)) {
    for (int i=0; i<nn; i++) {
      const int j = i * ndf+3;
      Vector3D v {ul[j+0], ul[j+1], ul[j+2]};
      ul.insert(i*ndf+3, ur[i].dLog(v), 1.0);
    }
  }

  return ul;
}


template <int nn, int ndf, typename IsoT>
// VectorND<nn*ndf>
void
EuclidFrameTransf<nn,ndf,IsoT>::pull(VectorND<nn*ndf>& ul, const Matrix3D& R, int op)
{

  // static VectorND<nn*ndf> ul;
  // for (int i=0; i<nn; i++) {
  //   const Vector &ddu = nodes[i]->getIncrDeltaDisp();
  //   for (int j = 0; j < ndf; j++) {
  //     ul[i*ndf+j] = ddu(j);
  //   }
  // }


  // 2) Isometry
  // 2.1) Rotation
  if (op & Transform::Rotation) {
    for (int i=0; i<nn; i++) {
      const int j = i * ndf;
      ul.insert(j  , R^Vector3D{ul[j+0], ul[j+1], ul[j+2]}, 1.0);
      ul.insert(j+3, R^Vector3D{ul[j+3], ul[j+4], ul[j+5]}, 1.0);
    }
  }

  // 2.2) Projection
  if (op & Transform::Adjoint)
  {
    const Vector3D wr = basis.getRotationVariation(ndf, &ul[0]);
    const Vector3D dc = basis.getPositionVariation(ndf, &ul[0]);

    for (int i=0; i<nn; i++) {
      Vector3D ui = this->getNodePosition(i);
      ul.assemble(i*ndf+0, dc, -1.0);
      ul.assemble(i*ndf+0, ui.cross(wr), 1.0);
      ul.assemble(i*ndf+3, wr, -1.0);
    }
  }

  // -) Offsets
  if constexpr (ndf >= 6)
    if (offsets && (offset_flags&OffsetLocal)) [[unlikely]] {
      const std::array<Vector3D, nn>& offset = *offsets;
      for (int i=0; i<nn; i++) {

        const int j = i * ndf;
        Vector3D w {ul[j+3],ul[j+4],ul[j+5]};

        ul.assemble(j, offset[i].cross(w), -1.0);
      }
    }

  // 3) Element parameters
  if (1) { // !(offset_flags & LogIter)) {
    for (int i=0; i<nn; i++) {
      const int j = i * ndf+3;
      Vector3D v {ul[j+0], ul[j+1], ul[j+2]};
      ul.insert(i*ndf+3, ur[i].dLog(v), 1.0);
    }
  }

  // return ul;
}

//
// Push
//


template <int nn, int ndf, typename IsoT>
int
EuclidFrameTransf<nn,ndf,IsoT>::push(VectorND<nn*ndf>&p, int op)
{
  VectorND<nn*ndf>& pa = p;

  // 3) Element Parameters
  if ((op & Transform::Logarithm) && !skip_log_iter ) { // (offset_flags & LogIter)) {
    for (int i=0; i<nn; i++) {
      Vector3D m {pa[i*ndf + 3], pa[i*ndf + 4], pa[i*ndf + 5]};
      pa.insert(i*ndf + 3, ur[i].dLog()^m, 1.0);
    }
  }


  // 2.1) Isometry
  if (op & Transform::Adjoint) {

    // a) Sum of moments: m = sum_i mi + sum_i (xi x ni)
    Vector3D M{}, N{};
    for (int i=0; i<nn; i++) {
      const int l = i*ndf;
      Vector3D ni {pa[l+0], pa[l+1], pa[l+2]};
      // m += mi
      for (int j=0; j<3; j++)
        M[j] += pa[i*ndf+3+j];
      // m += xi x ni
      M += this->getNodeLocation(i).cross(ni);
      N += ni;
    }
    // b) Adjust
    for (int i=0; i<nn; i++) {
      pa.template assemble<6>(i*ndf, basis.getRotationGradient(i)^M, -1.0);
      pa.template assemble<6>(i*ndf, basis.getTranslationGradient(i)^N, -1.0);
    }
  }

  // if (op == Transform::Isometry)
  //   return 0;

  // 2.2) Rotation
  const Matrix3D& R = basis.getRotation();
  if (op & Transform::Rotation) {
    for (int i=0; i<nn; i++) {
      const int base = i * ndf;
      pa.insert(base,   R*Vector3D{pa[base  ], pa[base+1], pa[base+2]}, 1.0);
      pa.insert(base+3, R*Vector3D{pa[base+3], pa[base+4], pa[base+5]}, 1.0);
    }
  }

  // 1) Offset
  if (offsets != nullptr) [[unlikely]] {
    const std::array<Vector3D, nn> roffsets = this->getCurrentOffsets();
    for (int i=0; i<nn; i++) {
      const int j = i * ndf;
      const Vector3D ni {pa[j+0], pa[j+1], pa[j+2]};
      pa.assemble(j+3, roffsets[i].cross(ni), 1.0);
    }
  }
  return 0;
}


template <int nn, int ndf, typename IsoT>
VectorND<6>
EuclidFrameTransf<nn,ndf,IsoT>::getWrench(const VectorND<nn*ndf>& p)
{
  Vector3D M{}, N{};
  for (int i=0; i<nn; i++) {
    const int l = i*ndf;
    Vector3D ni {p[l+0], p[l+1], p[l+2]};
    Vector3D mi = //(offset_flags & LogIter) ?
      Vector3D{p[l+3], p[l+4], p[l+5]} ;//:
      // ur[i].dLog()^Vector3D{p[l+3], p[l+4], p[l+5]};
    // m += mi
    M += mi;
    // m += xi x ni
    M += this->getNodeLocation(i).cross(ni);
    N += ni;
  }
  VectorND<6> w;
  for (int i=0; i<3; i++) {
    w[i]   = N[i];
    w[i+3] = M[i];
  }
  return w;
}


template <int nn, int ndf, typename IsoT>
int
EuclidFrameTransf<nn,ndf,IsoT>::push(MatrixND<nn*ndf,nn*ndf>&kb,
                                     const VectorND<nn*ndf>& pb, 
                                     int op)
{

  const double flip = getenv("XARA_FLIP")? 1.0 : -1.0;

  VectorND<nn*ndf> p = pb;
  for (int i=0; i<nn; i++) {
    for (int j=0; j<ndf-6; j++) {
      p(i*ndf +6 + j) = 0.0;
    }
  }

  MatrixND<nn*ndf,nn*ndf> Kb = kb;
  // 1) Element parameters
  if ((op & Transform::Logarithm)  && !skip_log_iter) { // (offset_flags & LogIter)) {
  // if (op != Transform::Rotation) {// && (op != Transform::Bubnov)) {//!(offset_flags & LogIter)) {
    for (int i=0; i<nn; i++) {
      Vector3D m{pb[i*ndf+3], pb[i*ndf+4], pb[i*ndf+5]};
      AxisAngle& ur_i = (offset_flags & LogIter) ? ddur[i] : ur[i];
      const Matrix3D Ai = ur_i.dLog();
      p.insert(i*ndf+3, Ai^m, 1.0);

      const Matrix3D kg = ur_i.ddLog(m);

      for (int j=0; j<nn; j++) {
        AxisAngle& ur_j = (offset_flags & LogIter) ? ddur[j] : ur[j];
        const Matrix3D Aj = ur_j.dLog();
        // loop over 3x3 blocks for n and m
        for (int k=0; k<2; k++) {
          for (int l=0; l<2; l++) {
            if (k == 0 && l == 0)
              continue;

            Matrix3D Kab {{
              Kb(i*ndf+3*k+0, j*ndf+3*l  ), Kb(i*ndf+3*k+1, j*ndf+3*l  ), Kb(i*ndf+3*k+2, j*ndf+3*l  ),
              Kb(i*ndf+3*k+0, j*ndf+3*l+1), Kb(i*ndf+3*k+1, j*ndf+3*l+1), Kb(i*ndf+3*k+2, j*ndf+3*l+1),
              Kb(i*ndf+3*k+0, j*ndf+3*l+2), Kb(i*ndf+3*k+1, j*ndf+3*l+2), Kb(i*ndf+3*k+2, j*ndf+3*l+2)
            }};
            if (k == 1)
              Kab = Ai^Kab; // row rotation block
            if (l == 1)
              Kab = Kab*Aj; // column rotation block

            Kb.insert(Kab, i*ndf+3*k, j*ndf+3*l, 1.0);
            if (i == j && k == 1 && l == 1)
              Kb.assemble(kg, i*ndf+3*k, j*ndf+3*l, 1.0);
          }
        }
      }
    }
  }

  //
  // 2) Isometry
  //
  // 2.2) Projection
  // 2.2.1) Kl = A ^ k * A ...
  // if ((op & Transform::Adjoint) || (op & Transform::Tangent))
  {
    MatrixND<nn*ndf,nn*ndf>& Kl = kb;
    const MatrixND<nn*ndf,nn*ndf> A = getProjection();
    if (op == Transform::Bubnov)
      Kl = A^kb;
    else
      Kl.addMatrixTripleProduct(0, A, Kb, 1);

    const VectorND<nn*ndf> Ap = A^p;
#if 0
    // 2.2.2) Kl += Kw * A
    {
      VectorND<nn*6> qwx{};
      for (int i=0; i<nn; i++)
        for (int j=0; j<6; j++)
          qwx[i*6+j] = p[i*ndf+j] - Ap[i*ndf+j];


      if (op != Transform::Bubnov) [[likely]] {
        MatrixND<nn*ndf,nn*ndf> Kb{};
        Kb.zero();
        if constexpr (ndf == 6) {
          Kl.addMatrixProduct(basis.getRotationJacobian(qwx), A, 1.0);
        }
        else {
          const MatrixND<6*nn,6*nn> Kw = basis.getRotationJacobian(qwx);
          Kb.assemble(Kw.template extract<0, 6,  0, 6>(),   0,   0, 1.0);
          Kb.assemble(Kw.template extract<0, 6,  6,12>(),   0, ndf, 1.0);
          Kb.assemble(Kw.template extract<6,12,  0, 6>(), ndf,   0, 1.0);
          Kb.assemble(Kw.template extract<6,12,  6,12>(), ndf, ndf, 1.0);
          Kl.addMatrixProduct(Kb, A, 1.0);
        }
      }
    }
#else 
    {
      const VectorND<6> pw = this->getWrench(p);
      // VectorND<nn*ndf> Ap = A^p;
      if (op != Transform::Bubnov) [[likely]] {

        if constexpr (ndf == 6) {
          Kl -= basis.getHessian(pw);
        }
        else {
          const MatrixND<6*nn,6*nn> Kw = basis.getHessian(pw);
          Unroll<0,nn>([&](auto aa) {
            constexpr int a = decltype(aa)::value;
            Unroll<0,nn>([&](auto bb) {
              constexpr int b = decltype(bb)::value;
              Kl.assemble(
                Kw.template extract<6*a, 6*(a+1), 6*b, 6*(b+1)>(),
                a*ndf,
                b*ndf,
                -1.0
              );
            });
          });
        }
      }
    }
#endif
    //
    // Kl += -W'*Pn'*A  - Pnm * W
    //
    {
      MatrixND<nn*ndf,nn*ndf> Kb{};
      Kb.zero();
      for (int j=0; j<nn; j++) {
        const MatrixND<3,6> Gj = basis.getRotationGradient(j);
        for (int i=0; i<nn; i++) {
          Kb.assemble(Hat(&p[i*ndf+0])*Gj,  i*ndf+0, j*ndf, -1.0);

          // Kl += -Pnm*W
          Kl.assemble(Hat(&Ap[i*ndf+0])*Gj, i*ndf+0, j*ndf, -1.0);
          Kl.assemble(Hat(&Ap[i*ndf+3])*Gj, i*ndf+3, j*ndf, -1.0);
        }
      }
      if (op != Transform::Bubnov)
        Kl.addMatrixTransposeProduct(1.0, Kb, A,  -1.0*flip);
    }
  } // Projection


  //
  // 2.1) Rotation
  //
  // Kl = diag(R) * Kl * diag(R)^T
  FrameTransform<nn,ndf>::pushRotation(kb, basis.getRotation());

  // 1) Offsets
  if (op & Transform::Offset && offsets != nullptr) [[unlikely]] {
    //
    // a) a'*k*a
    //
    const std::array<Vector3D, nn> roffsets = this->getCurrentOffsets();
    this->pushOffsets(kb, roffsets);
    //
    // b) += kg
    //
    p = pb;
    this->push(p, Transform::Total);// Transform::Adjoint|Transform::Rotation|Transform::Offset);
    for (int i=0; i<nn; i++) {
      kb.assemble(Hat(&p[i*ndf])*Hat(roffsets[i]), i*ndf+3, i*ndf+3, 1.0);
    }
  }

  return 0;
}



//
// Sensitivity
//
template <int nn, int ndf, typename IsoT>
bool
EuclidFrameTransf<nn,ndf,IsoT>::isShapeSensitivity()
{
  int nodeParameterI = nodes[   0]->getCrdsSensitivity();
  int nodeParameterJ = nodes[nn-1]->getCrdsSensitivity();
  // TODO(sensitivity): implement dvz

  return (nodeParameterI != 0 || nodeParameterJ != 0);
}


template <int nn, int ndf, typename IsoT>
double
EuclidFrameTransf<nn,ndf,IsoT>::getLengthGrad()
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

template <int nn, int ndf, typename IsoT>
double
EuclidFrameTransf<nn,ndf,IsoT>::getd1overLdh()
{
  return -getLengthGrad()/(L*L);
}


template <int nn, int ndf, typename IsoT>
void
EuclidFrameTransf<nn,ndf,IsoT>::pushGrad(VectorND<nn*ndf>& dp,
                                    VectorND<nn*ndf>& pl)
{  
  //
  // dp += T_{gl} dpl
  //
  double dL = this->getLengthGrad();
  double doneOverL = -dL/(L*L);

  constexpr Vector3D iv{1, 0, 0};
#if 1
  // 1.1) Sum of moments: m = sum_i mi + sum_i (xi x ni)
  Vector3D m{};
  for (int i=0; i<nn; i++) {
    // m += mi
    for (int j=0; j<3; j++)
      m[j] += pl[i*ndf+3+j];

    const Vector3D n = Vector3D{pl[i*ndf+0], pl[i*ndf+1], pl[i*ndf+2]};
    m.addVector(1, iv.cross(n), double(i)/double(nn-1)*L);
  }
  const Vector3D ixm = iv.cross(m);

  // 1.2) Adjust force part
  for (int i=0; i<nn; i++)
    dp.assemble(i*ndf,  ixm,  (i? 1.0:-1.0)*doneOverL);
#else 
#endif

  // 2) Rotate and do joint offsets

  //
  // dp = T_{lg}' pl
  //
  // int dv = 0; // TODO
  // int di = nodes[0]->getCrdsSensitivity();
  // int dj = nodes[1]->getCrdsSensitivity();

  Matrix3D R = basis.getRotation();
  for (int i=0; i<nn; i++) {
    const int base = i * ndf;
    dp.assemble(base,    R*Vector3D{dp[base  ], dp[base+1], dp[base+2]}, 1.0);
    dp.assemble(base+3,  R*Vector3D{dp[base+3], dp[base+4], dp[base+5]}, 1.0);
  }

  this->push(pl, Transform::Adjoint);
  Matrix3D dR = basis.getRotationSensitivity(nodes);
  
  for (int i=0; i<nn; i++) {
    const int base = i * ndf;
    dp.assemble(base,   dR*Vector3D{pl[base  ], pl[base+1], pl[base+2]}, 1.0);
    dp.assemble(base+3, dR*Vector3D{pl[base+3], pl[base+4], pl[base+5]}, 1.0);
  }


  return;
}


template <int nn, int ndf, typename IsoT>
void
EuclidFrameTransf<nn,ndf,IsoT>::pullFixedGrad(VectorND<nn*ndf>& du)
{
  //
  // dub += (T_{bl}' T_{lg}  +   T_{bl} T_{lg}') * ug
  //

  //
  // Form ug
  //
  VectorND<nn*ndf> ug;
  for (int i = 0; i < nn; i++) {
    const Vector& u = nodes[i]->getTrialDisp();
    for (int j = 0; j < ndf; j++) {
      ug[i*ndf+j] = u(j);
    }
  }


  // TODO: Sensitivity

  //
  // du = Tbl dR^ug
  {
    Matrix3D dR = basis.getRotationSensitivity(nodes);
  
    VectorND<nn*ndf> u1 = ug;
    EuclidFrameTransf<nn,ndf,IsoT>::pull(u1, dR, Transform::Total);
    du += u1;
  }

  {
    // double dL = this->getLengthGrad();
    // double doneOverL = -dL/(L*L);
    // double length = L;
    // L = 1/doneOverL;
    // VectorND<nn*ndf> u2 = ug;
    // EuclidFrameTransf<nn,ndf,IsoT>::pull(u2, R, offsets, offset_flags);
    // L = length; // restore
    // du += u2;
  }
  return;
}


template <int nn, int ndf, typename IsoT>
void
EuclidFrameTransf<nn,ndf,IsoT>::pullTotalGrad(VectorND<nn*ndf>& du, int gradNumber)
{
  for (int n=0; n<nn; n++)
    for (int i = 0; i < ndf; i++) {
      du[i + n*ndf] = nodes[n]->getDispSensitivity((i + 1), gradNumber);
    }

  // dub = T_{bl} T_{lg} * ug'
  const Matrix3D R = basis.getRotation();
  this->pull(du, R, Transform::Total);

  // dub += (T_{bl}' T_{lg}  +   T_{bl} T_{lg}') * ug
  this->pullFixedGrad(du);

  return;
}


template <int nn, int ndf, typename IsoT>
void
EuclidFrameTransf<nn,ndf,IsoT>::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"EuclidFrameTransf\"";
    s << ", \"vecxz\": [" 
      << vz[0] << ", " 
      << vz[1] << ", "
      << vz[2] << "]";

    //
    // Rotation parameters
    //
    s << ", \"internal_rotation\": ";
    if (offset_flags & LogIter)
      s << "\"Iter\"";
    else if (offset_flags & LogIncr)
      s << "\"Incr\"";
    else
      s << "\"Init\"";
    
    s << ", \"external_rotations\": [";
    for (int i=0; i<nn; i++) {
      switch (external_rotation_type[i]) {
        case Rotations::Parameters::None:
          s << "\"None\"";
          break;
        case Rotations::Parameters::Init:
          s << "\"Init\"";
          break;
        case Rotations::Parameters::Incr:
          s << "\"Incr\"";
          break;
        case Rotations::Parameters::Iter:
          s << "\"Iter\"";
          break;
      }
      if (i < nn-1)
        s << ", ";
    }
    s << "]";

    //
    if (offsets != nullptr) {
      s << ", \"offset_basis\": " << (offset_flags & OffsetLocal ? "\"local\"" : "\"global\"");
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
}

} // namespace OpenSees