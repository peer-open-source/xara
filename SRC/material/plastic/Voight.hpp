//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
#pragma once
#include <MatrixND.h>
#include <VectorND.h>
#include <cmath>

namespace Voight {
enum class Basis {
  Stress, Strain, Mixed
};

using namespace OpenSees;

static inline double 
Trace(const VectorND<6> &v)
{
  return v(0) + v(1) + v(2);
}

static inline VectorND<6> 
Dev(const VectorND<6> &v)
{
  VectorND<6> dev = v;
  double mean = (v(0) + v(1) + v(2)) / 3.0;
  dev(0) -= mean;
  dev(1) -= mean;
  dev(2) -= mean;
  return dev;
}


static inline double 
Dot(const VectorND<6> &a, const VectorND<6> &b)
{
  // the last 3 components contain the symmetric terms.
  return       a(0)*b(0) + a(1)*b(1) + a(2)*b(2)
        + 2.0*(a(3)*b(3) + a(4)*b(4) + a(5)*b(5));
}

// static inline VectorND<6>
// Mul(const MatrixND<6,6> &A, const VectorND<6> &v)
// {
//   VectorND<6> result{};
//   for (int i=0;i<6;i++)
//     for (int j=0;j<6;j++)
//       result(i) += A(i,j)*v(j);
//   return result;
// }

static inline double 
Oct(const VectorND<6> &v)
{
  return 0;
}

static inline void 
Bun(MatrixND<6,6>& P,
    const VectorND<3>& n1,
    const VectorND<3>& n2,
    const VectorND<3>& n3,
    const VectorND<3>& n4,
    const double scale)
{
    const VectorND<6> pi = {
        n1[0]*n2[0], n1[1]*n2[1], n1[2]*n2[2],
        n1[0]*n2[1], n1[1]*n2[2], n1[2]*n2[0]
    };
    const VectorND<6> pj = {
        n3[0]*n4[0], n3[1]*n4[1], n3[2]*n4[2],
        n3[0]*n4[1], n3[1]*n4[2], n3[2]*n4[0]
    };
    P.addTensorProduct(pi, pj, scale);
}


static inline double
J2(const VectorND<6> &v)
{
  //J2 = Dot(Dev(v), Dev(v))/2.0 == NormDev(v)/2
  return ( std::pow((v(0) - v(1)),2) 
        +  std::pow((v(0) - v(2)),2)
        +  std::pow((v(1) - v(2)),2))/6.0
        +  std::pow(v(3),2) + std::pow(v(4),2) + std::pow(v(5),2);
}


// 2nd order Identity Tensor
static constexpr MatrixND<3,3> I1 {{
      1.0, 0.0, 0.0 ,
      0.0, 1.0, 0.0 ,
      0.0, 0.0, 1.0 
}};

// 4th order mixed variant identity tensor (51b)
static constexpr MatrixND<6,6> IImix {{
      1, 0, 0,  0, 0, 0 ,
      0, 1, 0,  0, 0, 0 ,
      0, 0, 1,  0, 0, 0 ,
      0, 0, 0,  1, 0, 0 ,
      0, 0, 0,  0, 1, 0 ,
      0, 0, 0,  0, 0, 1 ,
}};

// 4th order covariant identity tensor (51a)
static constexpr MatrixND<6,6> IIco {{
      1, 0, 0,  0, 0, 0 ,
      0, 1, 0,  0, 0, 0 ,
      0, 0, 1,  0, 0, 0 ,
      0, 0, 0,  2, 0, 0 ,
      0, 0, 0,  0, 2, 0 ,
      0, 0, 0,  0, 0, 2 ,
}};

// 4th order contravariant identity tensor (51a)
static constexpr MatrixND<6,6> IIcon {{
      1, 0, 0,   0 ,  0 ,  0  ,
      0, 1, 0,   0 ,  0 ,  0  ,
      0, 0, 1,   0 ,  0 ,  0  ,
      0, 0, 0,  0.5,  0 ,  0  ,
      0, 0, 0,   0 , 0.5,  0  ,
      0, 0, 0,   0 ,  0 , 0.5 ,
}};

// 4th order Volumetric Tensor (57), or Ivp = ivol ⊗ ivol
// IIvol = I1 tensor I1
static constexpr MatrixND<6,6> IoI {{
      1, 1, 1,  0, 0, 0 ,
      1, 1, 1,  0, 0, 0 ,
      1, 1, 1,  0, 0, 0 ,
      0, 0, 0,  0, 0, 0 ,
      0, 0, 0,  0, 0, 0 ,
      0, 0, 0,  0, 0, 0 ,
}};
static constexpr MatrixND<6,6> IIvol {{
      1, 1, 1,  0, 0, 0 ,
      1, 1, 1,  0, 0, 0 ,
      1, 1, 1,  0, 0, 0 ,
      0, 0, 0,  0, 0, 0 ,
      0, 0, 0,  0, 0, 0 ,
      0, 0, 0,  0, 0, 0 ,
}};

static constexpr VectorND<6> ivol {{1, 1, 1, 0, 0, 0}};

// 4th order Deviatoric Tensor
//
// Note:  this is the contravariant form!
//        usable for s^a = 2G * IIdev^ab * epsilon_b
// (Need a different form for s^a = IIdev ^a_b * sigma^a)
static constexpr MatrixND<6,6> IIdev {{
       2./3.,  -1./3.,  -1./3.,   0 ,  0 ,  0  ,
      -1./3.,   2./3.,  -1./3.,   0 ,  0 ,  0  ,
      -1./3.,  -1./3.,   2./3.,   0 ,  0 ,  0  ,
          0.,      0.,      0.,  0.5,  0 ,  0  ,
          0.,      0.,      0.,   0 , 0.5,  0  ,
          0.,      0.,      0.,   0 ,  0 , 0.5 
}};

// 4th order contravariant deviatoric tensor (Id)
static constexpr MatrixND<6,6> IIdevCon = IIcon - 1./3.*IIvol;
//    0.6667   -0.3333   -0.3333         0         0         0
//   -0.3333    0.6667   -0.3333         0         0         0
//   -0.3333   -0.3333    0.6667         0         0         0
//         0         0         0    0.5000         0         0
//         0         0         0         0    0.5000         0
//         0         0         0         0         0    0.5000

// 4th order covariant deviatoric tensor
static constexpr MatrixND<6,6> IIdevCo  = IIco  - 1./3.*IIvol;

// 4th order mixed variant deviatoric tensor (Idp)
static constexpr MatrixND<6,6> IIdevMix = IImix - 1./3.*IIvol;
//     0.6667   -0.3333   -0.3333         0         0         0
//    -0.3333    0.6667   -0.3333         0         0         0
//    -0.3333   -0.3333    0.6667         0         0         0
//          0         0         0    1.0000         0         0
//          0         0         0         0    1.0000         0
//          0         0         0         0         0    1.0000

} // namespace Voight

namespace Mises {
}