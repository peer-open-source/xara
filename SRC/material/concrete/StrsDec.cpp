#include <array>
#include <Vector3D.h>
#include <Matrix3D.h>
#include <VectorND.h>
#include <MatrixND.h>
#include <cmath>
#include <eigen/SymEigDirect3D.h>


using namespace OpenSees;

inline double
Heaviside(double X) { 
  return X > 1e-16 ? 1.0 : (X < 0.0 ? 0.0 : 0.0); // 0.5 
}


static inline void 
addVoightTensorProduct(MatrixND<6,6>& P,
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


int 
StrsDecA(const VectorND<6> &sig, 
        VectorND<6> &sigpos,
        MatrixND<6,6>* Qpos)
{
  // convert to matrix form
  MatrixND<3,3> sigM{}; // = Vector2Tensor(sig);

  // sigM(0, 0) = sig[0];
  // sigM(1, 1) = sig[1];
  // sigM(2, 2) = sig[2];
  // sigM(0, 1) = sigM(1, 0) = sig[3];
  // sigM(1, 2) = sigM(2, 1) = sig[4];
  // sigM(0, 2) = sigM(2, 0) = sig[5];

  // find eigenvalues sigI and eigenvectors n
  Vector3D sigI{};
  std::array<Vector3D,3> V{};
  SymEigDirect3D<double,-1> eigen;
  eigen(sig[0], sig[3], sig[5], 
        sig[1], sig[4], 
        sig[2],
        sigI, V);

  // eigenvectors are n = [n1 n2 n3]
  Vector3D  & n1 = V[0], 
            & n2 = V[1], 
            & n3 = V[2];

  // Heaviside function
  double  H1 = Heaviside(sigI[0]), 
          H2 = Heaviside(sigI[1]), 
          H3 = Heaviside(sigI[2]);

  // positive stress

  // McCaulay of principal stresses, M1 = <sigI1>...
  double  M1 = H1*sigI[0],
          M2 = H2*sigI[1],
          M3 = H3*sigI[2];

  MatrixND<3,3> sigposM{};
  sigposM.addTensorProduct(n1, n1, M1);
  sigposM.addTensorProduct(n2, n2, M2);
  sigposM.addTensorProduct(n3, n3, M3);
  sigpos[0] = sigposM(0, 0);
  sigpos[1] = sigposM(1, 1);
  sigpos[2] = sigposM(2, 2);
  sigpos[3] = sigposM(0, 1);
  sigpos[4] = sigposM(1, 2);
  sigpos[5] = sigposM(2, 0);


  // Qpos, Qneg : linearization of stress projections, i.e.
  // d(sig_pos) = Qpos*d(sig)


  if (Qpos != nullptr) {    
    MatrixND<6,6>& Ppos = *Qpos;
    Ppos.zero();

    // Ppos : positive projection operator;
    addVoightTensorProduct(Ppos, n1,n1,  n1,n1,  H1);
    addVoightTensorProduct(Ppos, n2,n2,  n2,n2,  H2);
    addVoightTensorProduct(Ppos, n3,n3,  n3,n3,  H3);

    // p12(x)p12
    double term12 = (fabs(sigI[0] - sigI[1]) <= 1e-16) ? H1 : (M1 - M2)/(sigI[0] - sigI[1]);
    addVoightTensorProduct(Ppos, n1,n2,  n1,n2,  0.25*2*term12); // 1212
    addVoightTensorProduct(Ppos, n1,n2,  n2,n1,  0.25*2*term12); // 1221
    addVoightTensorProduct(Ppos, n2,n1,  n1,n2,  0.25*2*term12); // 2112
    addVoightTensorProduct(Ppos, n2,n1,  n2,n1,  0.25*2*term12); // 2121
    // p13(x)p13
    double term13 = (fabs(sigI[0] - sigI[2]) <= 1e-16) ? H1 : (M1 - M3)/(sigI[0] - sigI[2]);
    addVoightTensorProduct(Ppos, n1,n3,  n1,n3,  0.25*2*term13); // 1313
    addVoightTensorProduct(Ppos, n1,n3,  n3,n1,  0.25*2*term13); // 1321
    addVoightTensorProduct(Ppos, n3,n1,  n1,n3,  0.25*2*term13); // 3113
    addVoightTensorProduct(Ppos, n3,n1,  n3,n1,  0.25*2*term13); // 3121
    // p23(x)p23
    double term23 = (fabs(sigI[1] - sigI[2]) <= 1e-16) ? H2 : (M2 - M3)/(sigI[1] - sigI[2]);
    addVoightTensorProduct(Ppos, n2,n3,  n2,n3,  0.25*2*term23); // 2323
    addVoightTensorProduct(Ppos, n2,n3,  n3,n2,  0.25*2*term23); // 2332
    addVoightTensorProduct(Ppos, n3,n2,  n2,n3,  0.25*2*term23); // 3223
    addVoightTensorProduct(Ppos, n3,n2,  n3,n2,  0.25*2*term23); // 3232

    for (int i=0; i<6; i++)
      for (int j=3; j<6; j++)
        Ppos(i,j) = 2*Ppos(i,j); // 0.5*Ppos(i,j) + 0.5*Ppos(j,i);
  }
  return 0;
  // Qpos = Ppos + 2*(term12*(p12*p12') + term13*(p13*p13') + term23*(p23*p23'));
  // Qpos(:,4:6) = 2*Qpos(:,4:6);
}