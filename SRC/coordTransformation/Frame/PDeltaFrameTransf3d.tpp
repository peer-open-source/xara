//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//

//
// Adapted: Remo Magalhaes de Souza
//          04/2000
//
// Written: cmp
//
#include <Vector.h>
#include <Matrix.h>
#include <Matrix3D.h>
#include <Node.h>

#include "PDeltaFrameTransf3d.h"

namespace OpenSees {

template <int nn, int ndf>
PDeltaFrameTransf<nn,ndf>::PDeltaFrameTransf(int tag, 
                            const Vector3D &vecxz,
                            const std::array<Vector3D, nn> *offset,
                            int offset_flags,
                            bool ctan
                          )
  : FrameTransform<nn,ndf>(tag),
    offset_flags(offset_flags),
    linear(tag, vecxz, offset, offset_flags),
    consistent_tangent(ctan),
    vz(vecxz)
{

}


template <int nn, int ndf>
PDeltaFrameTransf<nn,ndf>::~PDeltaFrameTransf()
{

}


template <int nn, int ndf>
int
PDeltaFrameTransf<nn,ndf>::commit()
{
  return linear.commit();
}


template <int nn, int ndf>
int
PDeltaFrameTransf<nn,ndf>::revertToLastCommit()
{
  return linear.revertToLastCommit();
}

template <int nn, int ndf>
int
PDeltaFrameTransf<nn,ndf>::revertToStart()
{
  return linear.revertToStart();
}

template <int nn, int ndf>
int
PDeltaFrameTransf<nn,ndf>::initialize(std::array<Node*, nn>& new_nodes)
{
  return linear.initialize(new_nodes);
}

template <int nn, int ndf>
int
PDeltaFrameTransf<nn,ndf>::update() noexcept
{
  return linear.update();
}

template <int nn, int ndf>
int
PDeltaFrameTransf<nn,ndf>::getLocalAxes(Vector3D &XAxis, Vector3D &YAxis, Vector3D &ZAxis) const
{
  return linear.getLocalAxes(XAxis, YAxis, ZAxis);
}

template <int nn, int ndf>
double
PDeltaFrameTransf<nn,ndf>::getInitialLength()
{
  return linear.getInitialLength();
}

template <int nn, int ndf>
double
PDeltaFrameTransf<nn,ndf>::getDeformedLength()
{
  return linear.getDeformedLength();
}


template <int nn, int ndf>
int
PDeltaFrameTransf<nn,ndf>::push(VectorND<nn*ndf>&pl, int op)
{
  // Axial force
  const double N = pl[1*ndf+0];
  //
  // 1) let linear.push do redistribution (isometry only)
  //
  if (op & Transform::Adjoint) {

    linear.push(pl, Transform::Adjoint);

    // Include leaning column effects (P-Delta)
    const Vector3D Du = linear.getDelta()/linear.getInitialLength();

    pl[0*ndf+1] -= Du[1] * N;
    pl[1*ndf+1] += Du[1] * N;

    pl[0*ndf+2] -= Du[2] * N;
    pl[1*ndf+2] += Du[2] * N;
  }

  linear.push(pl, op&~Transform::Adjoint);

  return 0;
}


// template <int nn, int ndf>
// int
// PDeltaFrameTransf<nn,ndf>::push(MatrixND<nn*ndf,nn*ndf>& kl, 
//                                 const VectorND<nn*ndf> &pl, int op)
// {
//   double NoverL = pl[6] / linear.getInitialLength();

//   //
//   if (op & Transform::Adjoint) {
//     linear.push(kl, pl, Transform::Adjoint);

//     //
//     // Include geometric stiffness effects in local system;
//     //
//     // Kl += [ ]
//     kl(1, 1) += NoverL;
//     kl(2, 2) += NoverL;
//     kl(7, 7) += NoverL;
//     kl(8, 8) += NoverL;

//     kl(1, 7) -= NoverL;
//     kl(7, 1) -= NoverL;
//     kl(2, 8) -= NoverL;
//     kl(8, 2) -= NoverL;
//   }

//   return linear.push(kl, pl, op&~Transform::Adjoint);
// }
template <int nn, int ndf>
int
PDeltaFrameTransf<nn,ndf>::push(MatrixND<nn*ndf,nn*ndf>& kl,
                                const VectorND<nn*ndf> &pl,
                                int op)
{
  static_assert(nn == 2, "PDelta equivalence here assumes a two-node frame.");
  static_assert(ndf >= 6, "PDeltaFrameTransf requires at least 6 dofs per node.");

  if (op & Transform::Adjoint) {
    constexpr int i  = 0;
    constexpr int j  = 1;

    // const int ix = i*ndf + 0;
    constexpr int iy = i*ndf + 1;
    constexpr int iz = i*ndf + 2;

    constexpr int jx = j*ndf + 0;
    constexpr int jy = j*ndf + 1;
    constexpr int jz = j*ndf + 2;

    const double L = linear.getInitialLength();
    const double N = pl[jx];

    // First produce B^T Kb B in the local system.
    linear.push(kl, pl, Transform::Adjoint);

    // Row jx is dN/du_l, since local axial end force at node J is +N.
    VectorND<nn*ndf> dNdu{};
    for (int a = 0; a < nn*ndf; ++a)
      dNdu[a] = kl(jx, a);

    // Classical geometric stiffness:
    // (N/L) * (hy hy^T + hz hz^T).
    const double NoverL = N / L;

    kl(iy, iy) += NoverL;
    kl(iz, iz) += NoverL;
    kl(jy, jy) += NoverL;
    kl(jz, jz) += NoverL;

    kl(iy, jy) -= NoverL;
    kl(jy, iy) -= NoverL;
    kl(iz, jz) -= NoverL;
    kl(jz, iz) -= NoverL;

    if (consistent_tangent) {
      // ((Delta_y/L) hy + (Delta_z/L) hz) * dN/du_l.
      //
      // linear.getDelta()/L is (u_j - u_i)/L, so
      // Delta_y/L = (u_i_y - u_j_y)/L = -Du[1],
      // Delta_z/L = (u_i_z - u_j_z)/L = -Du[2].
      const Vector3D Du = linear.getDelta() / L;
      const double dyOverL = -Du[1];
      const double dzOverL = -Du[2];

      for (int a = 0; a < nn*ndf; ++a) {
        kl(iy, a) += dyOverL * dNdu[a];
        kl(jy, a) -= dyOverL * dNdu[a];

        kl(iz, a) += dzOverL * dNdu[a];
        kl(jz, a) -= dzOverL * dNdu[a];
      }
    }
  }

  return linear.push(kl, pl, op & ~Transform::Adjoint);
}

template <int nn, int ndf>
VectorND<nn*ndf>
PDeltaFrameTransf<nn,ndf>::getStateVariation()
{
  return linear.getStateVariation();
}

template <int nn, int ndf>
Vector3D
PDeltaFrameTransf<nn,ndf>::getNodePosition(int tag)
{
  return linear.getNodePosition(tag);
}

template <int nn, int ndf>
Vector3D
PDeltaFrameTransf<nn,ndf>::getNodeRotationLogarithm(int tag)
{
  return linear.getNodeRotationLogarithm(tag);
}

template <int nn, int ndf>
FrameTransform<nn,ndf> *
PDeltaFrameTransf<nn,ndf>::getCopy() const
{
  Vector3D e1, e2, e3;
  linear.getLocalAxes(e1, e2, e3);

  return new PDeltaFrameTransf(this->getTag(), 
                               e3, 
                               linear.getRigidOffsets(),
                               offset_flags,
                               consistent_tangent);
}


template <int nn, int ndf>
double
PDeltaFrameTransf<nn,ndf>::getLengthGrad()
{
  return linear.getLengthGrad();
}


template <int nn, int ndf>
void
PDeltaFrameTransf<nn,ndf>::Print(OPS_Stream &s, int flag)
{
  const std::array<Vector3D,nn> *offsets = linear.getRigidOffsets();
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"" << this->getClassType() << "\"";
    s << ", \"vecxz\": [" 
      << vz[0] << ", " 
      << vz[1] << ", "
      << vz[2] 
      << "]";
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
}

} // namespace OpenSees