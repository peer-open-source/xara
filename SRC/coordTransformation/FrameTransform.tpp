#pragma once
#include "FrameTransform.h"
#include <MatrixND.h>

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
    const std::array<Vector3D,nn> *offset = nullptr;
    constexpr int N = nn * ndf;

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
