#pragma once
#include "FrameTransform.h"
#include "Logging.h"
#include <MatrixND.h>

#if 1
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
#endif

static inline VectorND<6>
getBasic(VectorND<12>& ul, double oneOverL)
{
  VectorND<6> ub;

  ub[0] = ul[6] - ul[0]; // N

  double tmp;
  //                    ujy  -  uiy
  tmp   = -oneOverL * (ul[7] - ul[1]);
  ub[1] = ul[ 5] + tmp;   // Mzi
  ub[2] = ul[11] + tmp;   // Mzj

  //                   ujz  -  uiz
  tmp   = oneOverL * (ul[8] - ul[2]);
  ub[3] = ul[ 4] + tmp;   // Myi
  ub[4] = ul[10] + tmp;   // Myj
  ub[5] = ul[ 9] - ul[3]; // T

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



template <int nn, int ndf>
VectorND<nn*ndf> 
FrameTransform<nn,ndf>::pushConstant(const VectorND<nn*ndf>& pl) 
{
    Matrix3D R;
    Vector3D x, y, z;
    getLocalAxes(x, y, z);
    for (int i=0; i<3; i++) {
      R(i,0) = x[i];
      R(i,1) = y[i];
      R(i,2) = z[i];
    }
    // TODO: Rigid offsets
    constexpr int N = nn * ndf;
    const std::array<Vector3D,nn> *offset = this->getRigidOffsets();

    //
    // Initialize
    //
    VectorND<N> pg = pl;

    // (A) First pass: just do the direct transformations
    for (int i=0; i<nn; i++) {
        const int base = i * ndf;

        pg.insert(base,   R*Vector3D{pg[base  ], pg[base+1], pg[base+2]}, 1.0);
        pg.insert(base+3, R*Vector3D{pg[base+3], pg[base+4], pg[base+5]}, 1.0);
    }

    // (B) Second pass: add offset cross product ( r x F ) into the
    //     rotational DOFs [3..5]
    if constexpr (ndf >= 6)
        if (offset) {
            const std::array<Vector3D, nn>& offsets = *offset;
            for (int i=0; i<nn; i++) {

                const int base = i * ndf;
                const Vector3D ni {
                    pg[base+0],
                    pg[base+1],
                    pg[base+2]
                };

                // Add M = r x F
                pg.assemble(base+3, offsets[i].cross(ni), 1.0);
            }
        }

    return pg;
}

template<int nn, int ndf>
MatrixND<nn*ndf,nn*ndf>
FrameTransform<nn,ndf>::pushConstant(const MatrixND<nn*ndf,nn*ndf>& kl)
{
    //
    // Do diag(R)*M*diag(R)'
    //
    static MatrixND<nn*ndf,nn*ndf> Kg;
    Kg = kl;

    Matrix3D R;
    Vector3D x, y, z;
    getLocalAxes(x, y, z);
    for (int i=0; i<3; i++) {
      R(i,0) = x[i];
      R(i,1) = y[i];
      R(i,2) = z[i];
    }
    
    const Matrix3D RT = R.transpose();
    for (int i=0; i<nn; i++) {
      for (int j=0; j<nn; j++) {
        for (int k=0; k<2; k++) {
          for (int l=0; l<2; l++) {
            Matrix3D Kab {{
              {Kg(i*ndf+3*k+0, j*ndf+3*l  ), Kg(i*ndf+3*k+1, j*ndf+3*l  ), Kg(i*ndf+3*k+2, j*ndf+3*l  )},
              {Kg(i*ndf+3*k+0, j*ndf+3*l+1), Kg(i*ndf+3*k+1, j*ndf+3*l+1), Kg(i*ndf+3*k+2, j*ndf+3*l+1)},
              {Kg(i*ndf+3*k+0, j*ndf+3*l+2), Kg(i*ndf+3*k+1, j*ndf+3*l+2), Kg(i*ndf+3*k+2, j*ndf+3*l+2)}
            }};
            Kg.insert(R*(Kab*RT), i*ndf+3*k, j*ndf+3*l, 1.0);
          }
        }
      }
    }

    //
    // TODO: Rigid offsets
    //
    return Kg;
}
