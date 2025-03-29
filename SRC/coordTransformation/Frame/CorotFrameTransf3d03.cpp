//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the implementation for the
// CorotFrameTransf3d03 class. CorotFrameTransf3d03 is a Corotational
// transformation for a spatial frame element between the global
// and basic coordinate systems. The formulation is derived from
// Crisfield (1991) and employs a heuristic approximation to the
// logarithm on SO(3).
//
// Written: Claudio Perez
// Created: 05/2024
//
// Adapted from: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
//
#if __cplusplus >= 202302L
    // C++23 (and later) code
#  include <utility>
   using std::unreachable;
#else  
# define GCC_VERSION (__GNUC__ * 10000 \
                               + __GNUC_MINOR__ * 100 \
                               + __GNUC_PATCHLEVEL__)
# if GCC_VERSION >= 40500
# define unreachable()  __builtin_unreachable()
# else
# define unreachable() do {;} while(0)
# endif
#endif

#include <math.h>

#include <Node.h>
#include <Channel.h>
#include <Logging.h>
#include <CorotFrameTransf3d03.h>

#include <Triad.h>
#include <Vector.h>
#include <Matrix.h>
#include <VectorND.h>
#include <Vector3D.h>
#include <MatrixND.h>
#include <Matrix3D.h>
#include <Rotations.hpp>
#include "blk3x12x3.h"
#include "Orient/CrisfieldTransform.h"
using namespace OpenSees;

// MatrixND<12,3> CorotFrameTransf3d03::Lr2{};
// MatrixND<12,3> CorotFrameTransf3d03::Lr3{};

#undef OPS_STATIC
#define OPS_STATIC static 
#ifndef THREAD_LOCAL
# define THREAD_LOCAL static
#endif


// Permutation matrix (to renumber basic dof's)

// v = Tp * ul
//
//      |   thI  |   thJ  |
//      | x z -y | x 
//      | 0 1  2 | 3 4  5 | 6
// Tp=  [ 0 0  0   0 0  0   1;  0 // axial
//        0 1  0   0 0  0   0;  1 // rot z I
//        0 0  0   0 1  0   0;  2 // rot z J
//        0 0 -1   0 0  0   0;  3 // rot y I
//        0 0  0   0 0 -1   0;  4 // rot y J
//       -1 0  0   1 0  0   0]; 5 // torsion
//
//
// constexpr MatrixND<6,7> T = {{
//     {0,    0,    0,    0,    0,   -1 },
//     {0,    1,    0,    0,    0,    0 },
//     {0,    0,    0,   -1,    0,    0 },
//     {0,    0,    0,    0,    0,    1 },
//     {0,    0,    1,    0,    0,    0 },
//     {0,    0,    0,    0,   -1,    0 },
//     {1,    0,    0,    0,    0,    0 }}};

// So-called "Matrix of rigid body modes"
static MatrixND<6,12> T12_6 {
//   N Mz Mz My My  T
   {{0, 0, 0, 0, 0, 0}, // Ni
    {0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0}, // Ti
    {0, 0, 0, 1, 0, 0}, // My
    {0, 1, 0, 0, 0, 0}, // Mz

    {1, 0, 0, 0, 0, 0}, // Nj
    {0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 0},
    {0, 0, 0, 0, 0, 2}, // Tj
    {0, 0, 0, 0, 1, 0}, // My
    {0, 0, 1, 0, 0, 0}} // Mz
};

CorotFrameTransf3d03::CorotFrameTransf3d03(int tag, const Vector &vecInLocXZPlane,
                                       const Vector &rigJntOffsetI,
                                       const Vector &rigJntOffsetJ)
  : FrameTransform3d(tag, CRDTR_TAG_CorotFrameTransf3d),
    nodeIOffset(3), nodeJOffset(3),
    offset{nullptr, nullptr},
    L(0), Ln(0),
    nodeIInitialDisp(0), nodeJInitialDisp(0),
    initialDispChecked(false)
{
    alphaI.zero();
    alphaJ.zero();

    // Check vector that defines local xz plane
    if (vecInLocXZPlane.Size() != 3 ) {
        vz[0] = 0;
        vz[1] = 0;
        vz[2] = 1;
    }
    else {
      vz[0] = vecInLocXZPlane[0];
      vz[1] = vecInLocXZPlane[1];
      vz[2] = vecInLocXZPlane[2];
    }

    // check rigid joint offset for node I
    if (rigJntOffsetI.Size() != 3 ) {
        opserr << "CorotFrameTransf3d: Invalid rigid joint offset vector for node I\n";
        opserr << "Size must be 3\n";
        nodeIOffset.Zero();
    }
    else
        nodeIOffset = rigJntOffsetI;

    // Check rigid joint offset for node J
    if (rigJntOffsetJ.Size() != 3 ) {
        opserr << "CorotFrameTransf3d:  Invalid rigid joint offset vector for node J\n";
        opserr << "Size must be 3\n";
        nodeJOffset.Zero();
    }
    else
        nodeJOffset = rigJntOffsetJ;

    if (nodeIOffset.Norm() != 0)
      joint_offsets |= end_i;
    if (nodeJOffset.Norm() != 0)
      joint_offsets |= end_j;

    if ((joint_offsets & end_i) || (joint_offsets & end_j)) {
        opserr << "Rigid joint zones not implemented; ";
        opserr << "Using zero values\n"; 
        nodeIOffset.Zero();
        nodeJOffset.Zero();
    }
}



// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
CorotFrameTransf3d03::CorotFrameTransf3d03():
  FrameTransform3d(0, CRDTR_TAG_CorotFrameTransf3d),
  nodeIOffset(3), nodeJOffset(3),
  L(0.0), Ln(0.0),
  nodeIInitialDisp(0), nodeJInitialDisp(0), 
  initialDispChecked(false)
{

}


CorotFrameTransf3d03::~CorotFrameTransf3d03()
{
  if (nodeIInitialDisp != nullptr)
    delete [] nodeIInitialDisp;

  if (nodeJInitialDisp != nullptr)
    delete [] nodeJInitialDisp;
}


double
CorotFrameTransf3d03::getInitialLength()
{
  return L;
}


double
CorotFrameTransf3d03::getDeformedLength()
{
  return Ln;
}


FrameTransform3d *
CorotFrameTransf3d03::getCopy()
{

  CorotFrameTransf3d03 *theCopy =
    new CorotFrameTransf3d03(this->getTag(), Vector(vz), nodeIOffset, nodeJOffset);

  if (!theCopy) {
    opserr << "CorotFrameTransf3d03::getCopy() - out of memory creating copy\n";
    return 0;
  }

  theCopy->nodes[0]  = nodes[0];
  theCopy->nodes[1]  = nodes[1];
  theCopy->xAxis     = xAxis;
  theCopy->L         = L;
  theCopy->Ln        = Ln;
  theCopy->R0        = R0;
  theCopy->Q_pres[0] = Q_pres[0];
  theCopy->Q_pres[1] = Q_pres[1];
  theCopy->Q_past[0] = Q_past[0];
  theCopy->Q_past[1] = Q_past[1];
  theCopy->ul = ul;
  theCopy->ulcommit = ulcommit;

  return theCopy;
}


int
CorotFrameTransf3d03::revertToStart()
{
  ul.zero();
  Q_pres[0] = VersorFromMatrix(R0);
  Q_pres[1] = VersorFromMatrix(R0);

  alphaI.zero();
  alphaJ.zero();

  this->update();
  return 0;
}


int
CorotFrameTransf3d03::commitState()
{
  ulcommit  = ul;
  Q_past[0] = Q_pres[0];
  Q_past[1] = Q_pres[1];
  return 0;
}


int
CorotFrameTransf3d03::revertToLastCommit()
{
    // determine global displacement increments from last iteration
    const Vector &dispI = nodes[0]->getTrialDisp();
    const Vector &dispJ = nodes[1]->getTrialDisp();

    for (int k = 0; k < 3; k++) {
      alphaI(k) =  dispI(k+3);
      alphaJ(k) =  dispJ(k+3);
    }

    ul        = ulcommit;
    Q_pres[0] = Q_past[0];
    Q_pres[1] = Q_past[1];

    this->update();

    return 0;
}


int
CorotFrameTransf3d03::initialize(Node *nodeIPointer, Node *nodeJPointer)
{

    nodes[0] = nodeIPointer;
    nodes[1] = nodeJPointer;

    if ((!nodes[0]) || (!nodes[1])) {
      opserr << "\nCorotFrameTransf3d03::initialize";
      opserr << "\ninvalid pointers to the element nodes\n";
      return -1;
    }

    dX  = nodes[1]->getCrds() - nodes[0]->getCrds();

    // Add initial displacements at nodes
    if (initialDispChecked == false) {
      const Vector &nodeIDisp = nodes[0]->getDisp();
      const Vector &nodeJDisp = nodes[1]->getDisp();
      for (int i = 0; i<6; i++)
      if (nodeIDisp[i] != 0.0) {
        nodeIInitialDisp = new double [6];
        for (int j = 0; j<6; j++)
          nodeIInitialDisp[j] = nodeIDisp[j];
        i = 6;
      }

      for (int j = 0; j<6; j++)
        if (nodeJDisp[j] != 0.0) {
          nodeJInitialDisp = new double [6];
          for (int i = 0; i<6; i++)
            nodeJInitialDisp[i] = nodeJDisp[i];
          j = 6;
        }
      initialDispChecked = true;
    }

    // if (nodeIInitialDisp != nullptr) {
    //   dX[0] -= nodeIInitialDisp[0];
    //   dX[1] -= nodeIInitialDisp[1];
    //   dX[2] -= nodeIInitialDisp[2];
    // }

    // if (nodeJInitialDisp != nullptr) {
    //   dX[0] += nodeJInitialDisp[0];
    //   dX[1] += nodeJInitialDisp[1];
    //   dX[2] += nodeJInitialDisp[2];
    // }

    int error;
    static Vector XAxis(3);
    static Vector YAxis(3);
    static Vector ZAxis(3);

    // Set rotation matrix
    if ((error = this->getLocalAxes(XAxis, YAxis, ZAxis)))
      return error;

    // Compute initial pseudo-vectors for nodal triads
    Q_pres[0] = Q_pres[1] = VersorFromMatrix(R0);

    this->commitState();

    return 0;
}

void inline
CorotFrameTransf3d03::compTransfMatrixBasicGlobal(
                                                const Versor& Qbar,
                                                const Versor* Q_pres)
{

    // extract columns of rotation matrices
    const Triad r {MatrixFromVersor(Qbar)},
                rI{MatrixFromVersor(Q_pres[0])},
                rJ{MatrixFromVersor(Q_pres[1])};
    const Vector3D 
      &e1  =  crs.getBasisE1(), // E[1],
      &e2  =  crs.getBasisE2(), // E[2],
      &e3  =  crs.getBasisE3(), // E[3],
      &r1  =  r[1], // .rotate(E1), 
      &r2  =  r[2], // .rotate(E2), 
      &r3  =  r[3], // .rotate(E3),
      &rI1 = rI[1], // .rotate(E1),
      &rI2 = rI[2], // .rotate(E2), 
      &rI3 = rI[3], // .rotate(E3),
      &rJ1 = rJ[1], // .rotate(E1), 
      &rJ2 = rJ[2], // .rotate(E2), 
      &rJ3 = rJ[3]; // .rotate(E3);

    // Compute the transformation matrix from the basic to the
    // global system
    //   A = (1/Ln)*(I - e1*e1');
    // Matrix3D A;
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        A(i,j) = (double(i==j) - e1[i]*e1[j])/Ln;

    // This must be called up here
    CrisfieldTransform::getLMatrix(A, e1, r1, r2, Lr2);
    CrisfieldTransform::getLMatrix(A, e1, r1, r3, Lr3);

    //               3 |             3            |     3    |           3              |
    //   T1 = [      O', (-S(rI3)*e2 + S(rI2)*e3)',        O',                        O']'; imx
    //   T2 = [(A*rI2)', (-S(rI2)*e1 + S(rI1)*e2)', -(A*rI2)',                        O']'; imz
    //   T3 = [(A*rI3)', (-S(rI3)*e1 + S(rI1)*e3)', -(A*rI3)',                        O']'; imy
    //
    //   T4 = [      O',                        O',        O', ( S(rJ2)*e3 - S(rJ3)*e2 )']'; jmx
    //   T5 = [(A*rJ2)',                        O', -(A*rJ2)', ( S(rJ1)*e2 - S(rJ2)*e1 )']'; jmz
    //   T6 = [(A*rJ3)',                        O', -(A*rJ3)', ( S(rJ1)*e3 - S(rJ3)*e1 )']'; jmy

    T.zero();

    //   T1 = [      O', (-S(rI3)*e2 + S(rI2)*e3)',        O', O']';

    // (-S(rI3)*e2 + S(rI2)*e3)
    Vector3D Se  = rI2.cross(e3);
    Se -= rI3.cross(e2);
    for (int i = 0; i < 3; i++)
      // T(jmx,i+3) =  -Se[i];
      T(imx,i+3) =  Se[i];

    //   T2 = [(A*rI2)', (-S(rI2)*e1 + S(rI1)*e2)', -(A*rI2)', O']';

    Vector3D At = A*rI2;

    // (-S(rI2)*e1 + S(rI1)*e2)'
    Se  = rI1.cross(e2);
    Se -= rI2.cross(e1);
    for (int i = 0; i < 3; i++) {
        T(imz,i  ) =  At[i];
        T(imz,i+3) =  Se[i];
        T(imz,i+6) = -At[i];
    }

    //   T3 = [(A*rI3)', (-S(rI3)*e1 + S(rI1)*e3)', -(A*rI3)', O']';

    At = A*rI3;
    
    // -S(rI3)*e1 + S(rI1)*e3
    Se  = rI1.cross(e3);
    Se -= rI3.cross(e1);
    for (int i = 0; i < 3; i++) {
        T(imy,i  ) =  At[i]*-1; // TODO: Check
        T(imy,i+3) =  Se[i]*-1;
        T(imy,i+6) = -At[i]*-1;
    }

    //   T4 = [      O', O',        O', (-S(rJ3)*e2 + S(rJ2)*e3)']';
    Se  = rJ2.cross(e3);
    Se -= rJ3.cross(e2);
    for (int i = 0; i < 3; i++)
      T(jmx, i+9) =  Se[i];   // S(rJ2)*e3 - S(rJ3)*e2

    // T5 = [(A*rJ2)', O', -(A*rJ2)', (-S(rJ2)*e1 + S(rJ1)*e2)']';
    At = A*rJ2;
    Se  = rJ1.cross(e2); 
    Se -= rJ2.cross(e1);
    for (int i = 0; i < 3; i++) {
        T(jmz, i  ) =  At[i];
        T(jmz, i+6) = -At[i];
        T(jmz, i+9) =  Se[i]; // (-S(rJ2)*e1 + S(rJ1)*e2)
    }

    // T6 = [(A*rJ3)', O', -(A*rJ3)', (-S(rJ3)*e1 + S(rJ1)*e3)']'
    At  = A*rJ3;
    Se  = rJ1.cross(e3);  // (-S(rJ3)*e1 + S(rJ1)*e3)
    Se -= rJ3.cross(e1);
    for (int i = 0; i < 3; i++) {
        T(jmy,i  ) =  At[i]*-1; // TODO: Check
        T(jmy,i+6) = -At[i]*-1;
        T(jmy,i+9) =  Se[i]*-1;
    }

    //
    // Second part
    //

    // T(:,1) += Lr3*rI2 - Lr2*rI3;
    // T(:,2) +=           Lr2*rI1; z
    // T(:,3) += Lr3*rI1          ; y

    // T(:,4) += Lr3*rJ2 - Lr2*rJ3;
    // T(:,5) += Lr2*rJ1          ; z    // ?????? check sign
    // T(:,6) += Lr3*rJ1          ; y    // ?????? check sign

    // Bending Z
    for (int i = 0; i < 12; i++) {
      double T1i = 0;
      for (int k=0; k<3; k++)
        T1i += Lr2(i,k)*rI1[k];
      T(imz,i) += T1i;
    }

    for (int i = 0; i < 12; i++) {
      double T4i = 0;
      for (int k=0; k<3; k++)
        T4i += Lr2(i,k)*rJ1[k]; // Lr[i];
      T(jmz,i) += T4i;
    }

    // Torsion
    for (int i = 0; i < 12; i++) {
      double T0i = 0;
      for (int k=0; k<3; k++)
        T0i += Lr3(i,k)*rI2[k] - Lr2(i,k)*rI3[k];
      // T(jmx,i) += -T0i;
      T(imx,i) += T0i;
    }
    for (int i = 0; i < 12; i++) {
      double T3i = 0;
      for (int k=0; k<3; k++)
        T3i += Lr3(i,k)*rJ2[k] - Lr2(i,k)*rJ3[k];
      T(jmx,i) += T3i;
    }
    // Bending Y
    for (int i = 0; i < 12; i++) {
      double T2i = 0;
      for (int k=0; k<3; k++)
        T2i += Lr3(i,k)*rI1[k]; // Lr[i];
      T(imy,i) += T2i*-1; // TODO: Check
    }
    for (int i = 0; i < 12; i++) {
      double T5i = 0;
      for (int k=0; k<3; k++)
        T5i += Lr3(i,k)*rJ1[k]; // Lr[i];
      T(jmy,i) += T5i*-1; // TODO: Check
    }

    //
    //
    //
    for (int node=0; node < 2; node++)
      for (int j = 0; j < 3; j++) {
        const double c = 0.5 / std::cos(ul[(node? jmx : imx) + j]);
        for (int i = 0; i < 12; i++)
          T((node? jmx : imx) + j, i) *= c;
      }

    // Axial
    // T(:,7) = [-e1' O' e1' O']';
    for (int i = 0; i < 3; i++) {
        T(jnx,i  ) = -e1[i];
        T(jnx,i+6) =  e1[i];
    }

    // Combine torsion
    for (int i=0; i<12; i++) {
      T(jmx,i) -= T(imx,i);
      T(imx,i) = 0;
    }
}

//
// Set RI,RJ,Rbar, Ln, e and ul
//
int
CorotFrameTransf3d03::update()
{
    // determine global displacement increments from last iteration

    const Vector& dispI = nodes[0]->getTrialDisp();
    const Vector& dispJ = nodes[1]->getTrialDisp();

    //
    // Update state
    //
    // 1.1 Relative translation
    Vector3D dx = dX;// dx = dX + dJI;
    {
      // Relative translation
      for (int k = 0; k < 3; k++)
        dx[k] += dispJ(k) - dispI(k);


      // Calculate the deformed length
      Ln = dx.norm();

      if (Ln == 0.0) {
        opserr << "\nCorotFrameTransf3d03: deformed length is 0.0\n";
        return -2;
      }
    }

    // 1.2 Rotational displacement increments
    {
      Vector3D dAlphaI, dAlphaJ;

      for (int k = 0; k < 3; k++) {
        dAlphaI[k] =  dispI(k+3) - alphaI[k];
        dAlphaJ[k] =  dispJ(k+3) - alphaJ[k];
        alphaI[k]  =  dispI(k+3);
        alphaJ[k]  =  dispJ(k+3);
      }

      // Update the nodal rotations
      Q_pres[0] = VersorProduct(Q_pres[0],  Versor::from_vector(dAlphaI));
      Q_pres[1] = VersorProduct(Q_pres[1],  Versor::from_vector(dAlphaJ));
    }

    //
    // 2) Form transformation
    //

    crs.update(Q_pres[0], Q_pres[1], dx);
    // Form the transformation tangent
    this->compTransfMatrixBasicGlobal(crs.getReference(), Q_pres);

    //
    // 3) Local deformations
    //

    // Save previous state
    ulpr = ul;

    // Rotations
    {
      Matrix3D e = crs.getRotation();
      vr[0] = LogC90(e^MatrixFromVersor(Q_pres[0]));
      for (int i=0; i<3; i++)
        ul[imx+i] = vr[0][i];

      vr[1] = LogC90(e^MatrixFromVersor(Q_pres[1]));
      for (int i=0; i<3; i++)
        ul[jmx+i] = vr[1][i];
    }

    // Axial
    ul(inx) = 0;
    ul(jnx) = Ln - L;

    return 0;
}


inline VectorND<12>
CorotFrameTransf3d03::pushResponse(VectorND<12>&pl)
{
  // return T^pl;
  VectorND<12> pg{};
  for (int a = 0; a<2; a++) {
    VectorND<6> pa {pl(a*6+0), pl(a*6+1), pl(a*6+2), 
                    pl(a*6+3), pl(a*6+4), pl(a*6+5)};

    for (int b = 0; b<2; b++) {
      VectorND<6> pab = pushResponse(pa, a, b);
      pg.assemble(b*6, pab, 1.0);
    }
  }
  return pg;
}

VectorND<6>
CorotFrameTransf3d03::pushResponse(const VectorND<6>&pa, int a, int b)
{
  VectorND<6> pg{};
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
      pg[j] += T(a*6 + i, b*6+j) * pa[i];

  return pg;
}


MatrixND<12,12>
CorotFrameTransf3d03::pushConstant(const MatrixND<12,12>&kl)
{

  MatrixND<12,12> kg;
  static double RWI[3][3]{};
  static double RWJ[3][3]{};

#if 0

  if (joint_offsets & end_i) {
    // Compute RWI
    RWI[0][0] = -R0(1, 0) * nodeIOffset[2] + R0(2, 0) * nodeIOffset[1];
    RWI[1][0] = -R0(1, 1) * nodeIOffset[2] + R0(2, 1) * nodeIOffset[1];
    RWI[2][0] = -R0(1, 2) * nodeIOffset[2] + R0(2, 2) * nodeIOffset[1];

    RWI[0][1] = R0(0, 0) * nodeIOffset[2] - R0(2, 0) * nodeIOffset[0];
    RWI[1][1] = R0(0, 1) * nodeIOffset[2] - R0(2, 1) * nodeIOffset[0];
    RWI[2][1] = R0(0, 2) * nodeIOffset[2] - R0(2, 2) * nodeIOffset[0];

    RWI[0][2] = -R0(0, 0) * nodeIOffset[1] + R0(1, 0) * nodeIOffset[0];
    RWI[1][2] = -R0(0, 1) * nodeIOffset[1] + R0(1, 1) * nodeIOffset[0];
    RWI[2][2] = -R0(0, 2) * nodeIOffset[1] + R0(1, 2) * nodeIOffset[0];
  }


  if (joint_offsets & end_j) {
    // Compute RWJ
    RWJ[0][0] = -R0(1, 0) * nodeJOffset[2] + R0(2, 0) * nodeJOffset[1];
    RWJ[1][0] = -R0(1, 1) * nodeJOffset[2] + R0(2, 1) * nodeJOffset[1];
    RWJ[2][0] = -R0(1, 2) * nodeJOffset[2] + R0(2, 2) * nodeJOffset[1];

    RWJ[0][1] = R0(0, 0) * nodeJOffset[2] - R0(2, 0) * nodeJOffset[0];
    RWJ[1][1] = R0(0, 1) * nodeJOffset[2] - R0(2, 1) * nodeJOffset[0];
    RWJ[2][1] = R0(0, 2) * nodeJOffset[2] - R0(2, 2) * nodeJOffset[0];

    RWJ[0][2] = -R0(0, 0) * nodeJOffset[1] + R0(1, 0) * nodeJOffset[0];
    RWJ[1][2] = -R0(0, 1) * nodeJOffset[1] + R0(1, 1) * nodeJOffset[0];
    RWJ[2][2] = -R0(0, 2) * nodeJOffset[1] + R0(1, 2) * nodeJOffset[0];
  }
#endif

  // Transform local stiffness to global system
  // First compute kl*T_{lg}
  static double tmp[12][12];  // Temporary storage
  for (int m = 0; m < 12; m++) {
    tmp[m][0] = kl(m, 0) * R0(0, 0) + kl(m, 1) * R0(0, 1) + kl(m, 2) * R0(0, 2);
    tmp[m][1] = kl(m, 0) * R0(1, 0) + kl(m, 1) * R0(1, 1) + kl(m, 2) * R0(1, 2);
    tmp[m][2] = kl(m, 0) * R0(2, 0) + kl(m, 1) * R0(2, 1) + kl(m, 2) * R0(2, 2);

    tmp[m][3] = kl(m, 3) * R0(0, 0) + kl(m, 4) * R0(0, 1) + kl(m, 5) * R0(0, 2);
    tmp[m][4] = kl(m, 3) * R0(1, 0) + kl(m, 4) * R0(1, 1) + kl(m, 5) * R0(1, 2);
    tmp[m][5] = kl(m, 3) * R0(2, 0) + kl(m, 4) * R0(2, 1) + kl(m, 5) * R0(2, 2);

    if (joint_offsets & end_i) {
      tmp[m][3] += kl(m, 0) * RWI[0][0] + kl(m, 1) * RWI[1][0] + kl(m, 2) * RWI[2][0];
      tmp[m][4] += kl(m, 0) * RWI[0][1] + kl(m, 1) * RWI[1][1] + kl(m, 2) * RWI[2][1];
      tmp[m][5] += kl(m, 0) * RWI[0][2] + kl(m, 1) * RWI[1][2] + kl(m, 2) * RWI[2][2];
    }

    tmp[m][6] = kl(m, 6) * R0(0, 0) + kl(m, 7) * R0(0, 1) + kl(m, 8) * R0(0, 2);
    tmp[m][7] = kl(m, 6) * R0(1, 0) + kl(m, 7) * R0(1, 1) + kl(m, 8) * R0(1, 2);
    tmp[m][8] = kl(m, 6) * R0(2, 0) + kl(m, 7) * R0(2, 1) + kl(m, 8) * R0(2, 2);

    tmp[m][9]  = kl(m, 9) * R0(0, 0) + kl(m, 10) * R0(0, 1) + kl(m, 11) * R0(0, 2);
    tmp[m][10] = kl(m, 9) * R0(1, 0) + kl(m, 10) * R0(1, 1) + kl(m, 11) * R0(1, 2);
    tmp[m][11] = kl(m, 9) * R0(2, 0) + kl(m, 10) * R0(2, 1) + kl(m, 11) * R0(2, 2);

    if (joint_offsets & end_j) {
      tmp[m][ 9] += kl(m, 6) * RWJ[0][0] + kl(m, 7) * RWJ[1][0] + kl(m, 8) * RWJ[2][0];
      tmp[m][10] += kl(m, 6) * RWJ[0][1] + kl(m, 7) * RWJ[1][1] + kl(m, 8) * RWJ[2][1];
      tmp[m][11] += kl(m, 6) * RWJ[0][2] + kl(m, 7) * RWJ[1][2] + kl(m, 8) * RWJ[2][2];
    }
  }

  // Now compute T'_{lg}*(kl*T_{lg})
  for (int m = 0; m < 12; m++) {
    kg(0, m) = R0(0, 0) * tmp[0][m] + R0(0, 1) * tmp[1][m] + R0(0, 2) * tmp[2][m];
    kg(1, m) = R0(1, 0) * tmp[0][m] + R0(1, 1) * tmp[1][m] + R0(1, 2) * tmp[2][m];
    kg(2, m) = R0(2, 0) * tmp[0][m] + R0(2, 1) * tmp[1][m] + R0(2, 2) * tmp[2][m];

    kg(3, m) = R0(0, 0) * tmp[3][m] + R0(0, 1) * tmp[4][m] + R0(0, 2) * tmp[5][m];
    kg(4, m) = R0(1, 0) * tmp[3][m] + R0(1, 1) * tmp[4][m] + R0(1, 2) * tmp[5][m];
    kg(5, m) = R0(2, 0) * tmp[3][m] + R0(2, 1) * tmp[4][m] + R0(2, 2) * tmp[5][m];

    if (joint_offsets & end_i) {
      kg(3, m) += RWI[0][0] * tmp[0][m] + RWI[1][0] * tmp[1][m] + RWI[2][0] * tmp[2][m];
      kg(4, m) += RWI[0][1] * tmp[0][m] + RWI[1][1] * tmp[1][m] + RWI[2][1] * tmp[2][m];
      kg(5, m) += RWI[0][2] * tmp[0][m] + RWI[1][2] * tmp[1][m] + RWI[2][2] * tmp[2][m];
    }

    kg( 6, m) = R0(0, 0) * tmp[6][m] + R0(0, 1) * tmp[7][m] + R0(0, 2) * tmp[8][m];
    kg( 7, m) = R0(1, 0) * tmp[6][m] + R0(1, 1) * tmp[7][m] + R0(1, 2) * tmp[8][m];
    kg( 8, m) = R0(2, 0) * tmp[6][m] + R0(2, 1) * tmp[7][m] + R0(2, 2) * tmp[8][m];

    kg( 9, m) = R0(0, 0) * tmp[9][m] + R0(0, 1) * tmp[10][m] + R0(0, 2) * tmp[11][m];
    kg(10, m) = R0(1, 0) * tmp[9][m] + R0(1, 1) * tmp[10][m] + R0(1, 2) * tmp[11][m];
    kg(11, m) = R0(2, 0) * tmp[9][m] + R0(2, 1) * tmp[10][m] + R0(2, 2) * tmp[11][m];

    if (joint_offsets & end_j) {
      kg( 9, m) += RWJ[0][0] * tmp[6][m] + RWJ[1][0] * tmp[7][m] + RWJ[2][0] * tmp[8][m];
      kg(10, m) += RWJ[0][1] * tmp[6][m] + RWJ[1][1] * tmp[7][m] + RWJ[2][1] * tmp[8][m];
      kg(11, m) += RWJ[0][2] * tmp[6][m] + RWJ[1][2] * tmp[7][m] + RWJ[2][2] * tmp[8][m];
    }
  }

  return kg;
}


VectorND<12>
CorotFrameTransf3d03::pushConstant(const VectorND<12>&pl) const
{
  // transform vector from local to global coordinates

  VectorND<12> pg;
  for (int i=0; i<4; i++)
    for (int j=0; j<3; j++)
      pg[i*3+j] = R0(j,0)*pl[3*i] + R0(j,1)*pl[3*i+1] + R0(j,2)*pl[3*i+2];

//
// TODO(cmp): joint offsets
//
//if (joint_offsets & end_i) {
//  pg[3] += -nodeIOffset[2] * pg[1] + nodeIOffset[1] * pg[2];
//  pg[4] +=  nodeIOffset[2] * pg[0] - nodeIOffset[0] * pg[2];
//  pg[5] += -nodeIOffset[1] * pg[0] + nodeIOffset[0] * pg[1];
//}

//if (joint_offsets & end_j) {
//  pg[ 9] += -nodeJOffset[2] * pg[7] + nodeJOffset[1] * pg[8];
//  pg[10] +=  nodeJOffset[2] * pg[6] - nodeJOffset[0] * pg[8];
//  pg[11] += -nodeJOffset[1] * pg[6] + nodeJOffset[0] * pg[7];
//}

  return pg;
}



const Vector &
CorotFrameTransf3d03::getGlobalResistingForce(const Vector &pb, const Vector &p0)
{
    static VectorND<12> pg;
    static Vector wrapper(pg);

    // Transform resisting forces from the basic to local
    static VectorND<12> pl;
    pl = pushLocal(pb, L);

    // Transform from local to global
    pg = pushResponse(pl);

    // If there are no element loads present, just return
    if (p0 == 0.0)
      return wrapper;

    // Add end forces due to element p0 loads
    // assuming member loads are in local system
    static VectorND<12> pl0{0};
    pl0[0] = p0[0]; // N
    pl0[1] = p0[1]; // Viy
    pl0[7] = p0[2]; // Vjy
    pl0[2] = p0[3]; // Viz
    pl0[8] = p0[4]; // Vjz

    pg += pushConstant(pl0);

    return wrapper;
}



const Matrix &
CorotFrameTransf3d03::getInitialGlobalStiffMatrix(const Matrix &kb)
{
  // transform tangent from the basic system to local coordinates
  static Matrix kl(12,12);
  kl.addMatrixTripleProduct(0.0, T12_6, kb, 1.0);   // kl = Tp ^ kb * Tp;

  // transform tangent from local to global coordinates
  static Matrix kg(12,12);

  // compute the tangent matrix in global coordinates
  kg.addMatrixTripleProduct(0.0, T, kl, 1.0);

  return kg;
}


const Matrix &
CorotFrameTransf3d03::getGlobalStiffMatrix(const Matrix &kb, const Vector &pb)
{
  // transform kb from the basic system to local
  static MatrixND<12,12> kl;
  static Matrix Kl(kl);

  // transform basic stiffness to 7x7 Remo layout
  // Kl.addMatrixTripleProduct(0.0, Tp, kb, 1.0); // kl = Tp ^ kb * Tp;
  Kl.addMatrixTripleProduct(0.0, T12_6, kb, 1.0); //


  // transform resisting forces from the basic system to local
  VectorND<12> pl = pushLocal(pb, L);

  //
  // Transform kl from local to global system
  //
  THREAD_LOCAL MatrixND<12,12> kg;
  static Matrix Wrapper(kg);
  kg.addMatrixTripleProduct(0.0, T, kl, 1.0);
  
  this->addTangent(kg, pl);
  return Wrapper;
}


// do 
//    K = ag'*Km*ag + Kp
MatrixND<12,12>
CorotFrameTransf3d03::pushResponse(MatrixND<12,12>& kl, const VectorND<12>& pl)
{    
  MatrixND<12,12> K;
  K.addMatrixTripleProduct(0.0, T, kl, 1.0);

  // Add geometric part kg
  this->addTangent(K, pl);

  return K;
}

//
// Add geometric part of the transformation tangent
//
//  kg += T'*kl*T + ks1 + T * diag(m.*tan(thetal))*T' + ...
//         m(4)*(ks2r2t3_u3 + ks2r3u2_t2) + ...
//         m(2)*ks2r2t1 + m(3)*ks2r3t1 + ...
//         m(5)*ks2r2u1 + m(6)*ks2r3u1 + ...
//         ks3 + ks3' + ks4 + ks5;
int
CorotFrameTransf3d03::addTangent(MatrixND<12,12>& kg, const VectorND<12>& pl)
{    
    const Triad r {MatrixFromVersor(crs.getReference())},
                rI{MatrixFromVersor(Q_pres[0])},
                rJ{MatrixFromVersor(Q_pres[1])};
    const Vector3D 
      &e1  =  crs.getBasisE1(), // E[1],
      &e2  =  crs.getBasisE2(), // E[2],
      &e3  =  crs.getBasisE3(), // E[3],
      &r1  =  r[1], // .rotate(E1), 
      &r2  =  r[2], // .rotate(E2), 
      &r3  =  r[3], // .rotate(E3),
      &rI1 = rI[1], // .rotate(E1), 
      &rI2 = rI[2], // .rotate(E2), 
      &rI3 = rI[3], // .rotate(E3),
      &rJ1 = rJ[1], // .rotate(E1), 
      &rJ2 = rJ[2], // .rotate(E2), 
      &rJ3 = rJ[3]; // .rotate(E3);
    // NOTE[cmp] 
    // CorotFrameTransf3d03::compTransfMatrixBasicGlobal must be 
    // called first to set Lr1, Lr2 and T

    // Matrix3D A;
    // for (int i = 0; i < 3; i++)
    //   for (int j = 0; j < 3; j++)
    //     A(i,j) = (double(i==j) - e1[i]*e1[j])/Ln;
    // getLMatrix(A, e1, r1, r2, Lr2);
    // getLMatrix(A, e1, r1, r3, Lr3);

    //
    // Ksigma1
    //
    {
      const double N = -pl[0]; // Axial force
      // a=0
      kg.assemble(A, 0, 0,  N);
      kg.assemble(A, 0, 6, -N);
      // a=1
      kg.assemble(A, 6, 0, -N);
      kg.assemble(A, 6, 6,  N);
    }

    //
    // Ksigma3
    //
    //  ks3 = [o kbar2  |  o kbar4];
    //
    //  where
    //
    //    kbar2 = -Lr2*(m(3)*S(rI3) + m(1)*S(rI1)) + Lr3*(m(3)*S(rI2) - m(2)*S(rI1)) ;
    //
    //    kbar4 =  Lr2*(m(3)*S(rJ3) - m(4)*S(rJ1)) - Lr3*(m(3)*S(rJ2) + m(5)*S(rJ1));
    //
    // or
    //
    //  ks3 = [o ka+kb  |  o kc+kd];
    //      = [o ka     |  o kc] + [o kb  |  o kd];
    //
    //  where
    //
    //    ka = -Lr2*S(rI3)*m(3)  
    //         +Lr2*S(rI1)*m(1);
    //    kb =  Lr3*S(rI2)*m(3)  
    //         -Lr3*S(rI1)*m(2);
    //
    //    kc =  Lr2*S(rJ3)*m(3)
    //         -Lr2*S(rJ1)*m(4);
    //    kd = -Lr3*S(rJ2)*m(3)  
    //         +Lr3*S(rJ1)*m(5);

    static VectorND<6> m;
    m[0] =  0.5*pl[imx]/std::cos(ul(imx));
    m[2] = -0.5*pl[imy]/std::cos(ul(imy));
    m[1] =  0.5*pl[imz]/std::cos(ul(imz));

    m[3] =  0.5*pl[jmx]/std::cos(ul(jmx));
    m[5] = -0.5*pl[jmy]/std::cos(ul(jmy));
    m[4] =  0.5*pl[jmz]/std::cos(ul(jmz));


    static Matrix3D Sm;
    Sm.zero();
    Sm.addSpin(rI3,  m[3]);
    Sm.addSpin(rI1,  m[1]);
    static MatrixND<12,3> kbar;
    kbar.zero();
    kbar.addMatrixProduct(Lr2, Sm, -1.0);

    Sm.zero();
    Sm.addSpin(rI2,  m[3]);
    Sm.addSpin(rI1, -m[2]);
    kbar.addMatrixProduct(Lr3, Sm,  1.0);

    kg.assemble(kbar, 0, 3, 1.0);
    kg.assembleTranspose(kbar, 3, 0, 1.0);

    Sm.zero();
    Sm.addSpin(rJ3,  m[3]);
    Sm.addSpin(rJ1, -m[4]);
    kbar.zero();
    kbar.addMatrixProduct(Lr2, Sm, 1.0);

    Sm.zero();
    Sm.addSpin(rJ2, m[3]);
    Sm.addSpin(rJ1, m[5]);
    kbar.addMatrixProduct(Lr3, Sm,  -1.0);

    kg.assemble(kbar, 0, 9, 1.0);
    kg.assembleTranspose(kbar, 9, 0, 1.0);


    //
    // Ksigma4
    //
    {
      static Matrix3D ks33;
  
      ks33.zero();
      ks33.addSpinProduct(e2, rI3,  m[3]);
      ks33.addSpinProduct(e3, rI2, -m[3]);
      ks33.addSpinProduct(e2, rI1,  m[1]);
      ks33.addSpinProduct(e1, rI2, -m[1]);
      ks33.addSpinProduct(e3, rI1,  m[2]);
      ks33.addSpinProduct(e1, rI3, -m[2]);
      kg.assemble(ks33, 3, 3, 1.0);
    }

    //
    // Ksigma4
    //
    {
      static Matrix3D ks33;
      ks33.zero();
      ks33.addSpinProduct(e2, rJ3, -m[3]);
      ks33.addSpinProduct(e3, rJ2,  m[3]);
      ks33.addSpinProduct(e2, rJ1,  m[4]);
      ks33.addSpinProduct(e1, rJ2, -m[4]);
      ks33.addSpinProduct(e3, rJ1,  m[5]);
      ks33.addSpinProduct(e1, rJ3, -m[5]);
  
      kg.assemble(ks33, 9, 9, 1.0);
    }


    //
    // Ksigma5
    //
    //  Ks5 = [ Ks5_11   Ks5_12 | -Ks5_11   Ks5_14;
    //          Ks5_12'    O    | -Ks5_12'   O;
    //         -Ks5_11  -Ks5_12 |  Ks5_11  -Ks5_14;
    //          Ks5_14t     O   | -Ks5_14'   O];
    //
    //
    // v = (1/Ln)*(m(2)*rI2 + m(3)*rI3 + m(5)*rJ2 + m(6)*rJ3);
    //   = 1/Ln * (m[1]*rI2 + m[2]*rI3)
    //   + 1/Ln * (m[4]*rJ2 + m[5]*rJ3);
    //   = vi + vj
    //
    {
      OPS_STATIC Vector3D v;
      v.addVector(0.0, rI2, m[1]);
      v.addVector(1.0, rI3, m[2]);
      v.addVector(1.0, rJ2, m[4]);
      v.addVector(1.0, rJ3, m[5]);
      v /= Ln;

      // Ks5_11 = A*v*e1' + e1*v'*A + (e1'*v)*A;
      //        = A*vi*e1' + e1*vi'*A + (e1'*vi)*A
      //        + A*vj*e1' + e1*vj'*A + (e1'*vj)*A;
      //
      Matrix3D ks33;
      ks33.zero();
      ks33.addMatrix(A, e1.dot(v));

      static Matrix3D m33;
      m33.zero();
      m33.addTensorProduct(v, e1, 1.0);

      ks33.addMatrixProduct(A, m33, 1.0);

      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
          m33(i,j) = e1[i]*v[j];

      ks33.addMatrixProduct(m33, A, 1.0);

      kg.assemble(ks33, 0, 0,  1.0);
      kg.assemble(ks33, 0, 6, -1.0);
      kg.assemble(ks33, 6, 0, -1.0);
      kg.assemble(ks33, 6, 6,  1.0);
    }

    // Ks5_12 = -(m(2)*A*S(rI2) + m(3)*A*S(rI3));

    Matrix3D ks33;
    ks33.zero();
    ks33.addMatrixSpinProduct(A, rI2, -m[1]);
    ks33.addMatrixSpinProduct(A, rI3, -m[2]);

    kg.assemble(ks33, 0, 3,  1.0);
    kg.assemble(ks33, 6, 3, -1.0);
    kg.assembleTranspose(ks33, 3, 0,  1.0);
    kg.assembleTranspose(ks33, 3, 6, -1.0);

    //  Ks5_14 = -(m(5)*A*S(rJ2) + m(6)*A*S(rJ3));

    ks33.zero();
    ks33.addMatrixSpinProduct(A, rJ2, -m[4]);
    ks33.addMatrixSpinProduct(A, rJ3, -m[5]);

    kg.assemble(ks33, 0, 9,  1.0);
    kg.assemble(ks33, 6, 9, -1.0);

    kg.assembleTranspose(ks33, 9, 0,  1.0);
    kg.assembleTranspose(ks33, 9, 6, -1.0);

    // Ksigma -------------------------------
    OPS_STATIC Vector3D rm;

    rm = rI3;
    rm.addVector(1.0, rJ3, -1.0);
    kg.addMatrix(CrisfieldTransform::getKs2Matrix(A, e1, r1, Ln, r2, rm), m[3]);

//  rm = rJ2;
    rm.addVector(0.0, rJ2, -1.0);
    rm.addVector(1.0, rI2, -1.0);
    kg.addMatrix(CrisfieldTransform::getKs2Matrix(A, e1, r1, Ln, r3,  rm), m[3]);
    kg.addMatrix(CrisfieldTransform::getKs2Matrix(A, e1, r1, Ln, r2, rI1), m[1]);
    kg.addMatrix(CrisfieldTransform::getKs2Matrix(A, e1, r1, Ln, r3, rI1), m[2]);
    //
    kg.addMatrix(CrisfieldTransform::getKs2Matrix(A, e1, r1, Ln, r2, rJ1), m[4]);
    kg.addMatrix(CrisfieldTransform::getKs2Matrix(A, e1, r1, Ln, r3, rJ1), m[5]);

    //
    //  T' * diag (M .* tan(thetal))*T
    //

    for (int node=0; node<2; node++) {
      for (int k = 0; k < 3; k++) {
        const double factor =  pl[6*node+3+k] * std::tan(ul[(node ? jmx : imx) + k]);
        for (int i = 0; i < 12; i++) {
          const double Tki = T((node ? jmx : imx) + k,i);
          for (int j = 0; j < 12; j++)
            kg(i,j) += Tki * factor * T((node ? jmx : imx) + k, j);
        }
      }
    }

    return 0;
}


int
CorotFrameTransf3d03::getLocalAxes(Vector &XAxis, Vector &YAxis, Vector &ZAxis)
{
  ZAxis(0) = vz[0];
  ZAxis(1) = vz[1];
  ZAxis(2) = vz[2];
  if (nodes[0] == nullptr)
    return -1;

  static Vector dx(3);
  dx = (nodes[1]->getCrds() + nodeJOffset) - (nodes[0]->getCrds() + nodeIOffset);
  if (nodeIInitialDisp != 0) {
    dx(0) -= nodeIInitialDisp[0];
    dx(1) -= nodeIInitialDisp[1];
    dx(2) -= nodeIInitialDisp[2];
  }

  if (nodeJInitialDisp != 0) {
    dx(0) += nodeJInitialDisp[0];
    dx(1) += nodeJInitialDisp[1];
    dx(2) += nodeJInitialDisp[2];
  }

  // calculate the element length

  L = dx.Norm();

  if (L == 0.0) {
      opserr << "\nCorotFrameTransf3d03::computeElemtLengthAndOrien: 0 length\n";
      return -2;
  }

  // calculate the element local x axis components wrt to the global coordinates

  xAxis(0) = dx(0)/L;
  xAxis(1) = dx(1)/L;
  xAxis(2) = dx(2)/L;

  XAxis(0) = xAxis(0);
  XAxis(1) = xAxis(1);
  XAxis(2) = xAxis(2);

  //
  Vector3D yAxis = vz.cross(xAxis);

  const double ynorm = yAxis.norm();

  if (ynorm == 0.0) {
      opserr << "\nCorotFrameTransf3d03::getElementLengthAndOrientation";
      opserr << "\nvector v that defines plane xz is parallel to x axis\n";
      return -3;
  }

  yAxis /= ynorm;
  YAxis(0) = yAxis(0);
  YAxis(1) = yAxis(1);
  YAxis(2) = yAxis(2);

  Vector3D zAxis = xAxis.cross(yAxis);

  ZAxis(0) = zAxis(0);
  ZAxis(1) = zAxis(1);
  ZAxis(2) = zAxis(2);

  for (int i = 0; i < 3; i++) {
      R0(i,0) = xAxis[i];
      R0(i,1) = yAxis[i];
      R0(i,2) = zAxis[i];
  }

  return 0;
}


const Vector &
CorotFrameTransf3d03::getBasicTrialDisp()
{
  static Vector ub(6);
  LocalToBasic(ul, ub);
  return ub;
}


const Vector &
CorotFrameTransf3d03::getBasicIncrDeltaDisp()
{
  static VectorND<12> dul;

  // dul = ul - ulpr;
  dul = ul;
  dul.addVector(1.0, ulpr, -1.0);

  static Vector dub(6);
  LocalToBasic(dul, dub);
  return dub;
}


const Vector &
CorotFrameTransf3d03::getBasicIncrDisp()
{
  static VectorND<12> Dul;

  // Dul = ul - ulcommit;
  Dul = ul;
  Dul.addVector(1.0, ulcommit, -1.0);

  static Vector Dub(6);
  LocalToBasic(Dul, Dub);
  return Dub;
}


const Vector &
CorotFrameTransf3d03::getBasicTrialVel()
{
  opserr << "WARNING unimplemented method\n";
  static Vector dummy(6);
  return dummy;
}


const Vector &
CorotFrameTransf3d03::getBasicTrialAccel()
{
  opserr << "WARNING unimplemented method\n";
  static Vector dummy(6);
  return dummy;
}


const Matrix &
CorotFrameTransf3d03::getGlobalMatrixFromLocal(const Matrix &local)
{
  static Matrix Mg(12,12);
  blk3x12x3(R0, local, Mg);
  return Mg;
}

double
CorotFrameTransf3d03::getLengthGrad()
{
  const int di = nodes[0]->getCrdsSensitivity();
  const int dj = nodes[1]->getCrdsSensitivity();

  Vector3D dxi{0.0};
  Vector3D dxj{0.0};

  if (di != 0)
    dxi(di-1) = 1.0;
  if (dj != 0)
    dxj(dj-1) = 1.0;

  return 1/L * dX.dot(dxj - dxi);
}


void
CorotFrameTransf3d03::Print(OPS_Stream &s, int flag)
{

  if (flag == OPS_PRINT_CURRENTSTATE) {
      s << "\nFrameTransform: " << this->getTag() << " Type: CorotFrameTransf3d03";
      s << "\tvxz: " << Vector(vz);
      s << "\tnodeI Offset: " << nodeIOffset;
      s << "\tnodeJ Offset: " << nodeJOffset;
      return;
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
      s << OPS_PRINT_JSON_MATE_INDENT << "{";
      s << "\"name\": \"" << this->getTag() << "\", ";
      s << "\"type\": \"CorotFrameTransf3d03\"";
      s << ", \"vecxz\": [" 
        << vz(0) << ", " << vz(1) << ", " << vz(2) << "]";
      if (nodeIOffset != 0)
          s << ", \"iOffset\": [" << nodeIOffset[0] << ", " << nodeIOffset[1] << ", " << nodeIOffset[2] << "]";
      if (nodeJOffset != 0)
          s << ", \"jOffset\": [" << nodeJOffset[0] << ", " << nodeJOffset[1] << ", " << nodeJOffset[2] << "]";
      s << "}";
      return;
  }
}


const Vector &
CorotFrameTransf3d03::getPointGlobalCoordFromLocal(const Vector &xl)
{
  static Vector xg(3);
  opserr << " CorotFrameTransf3d03::getPointGlobalCoordFromLocal: not implemented yet" ;
  return xg;
}


const Vector &
CorotFrameTransf3d03::getPointGlobalDisplFromBasic(double xi, const Vector &uxb)
{
  static Vector uxg(3);
  opserr << " CorotFrameTransf3d03::getPointGlobalDisplFromBasic: not implemented yet" ;
  return uxg;
}


const Vector &
CorotFrameTransf3d03::getPointLocalDisplFromBasic(double xi, const Vector &uxb)
{
  static Vector uxg(3);
  opserr << " CorotFrameTransf3d03::getPointLocalDisplFromBasic: not implemented yet" ;

  return uxg;
}

int
CorotFrameTransf3d03::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static Vector data(48);
  if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
    opserr << " CorotFrameTransf3d03::recvSelf() - data could not be received\n" ;
    return -1;
  }

  for (int i = 0; i<7; i++)
    ulcommit[i] = data(i);

  for (int j = 0; j<3; j++) {
    Q_past[0].vector[j] = data(7+j);
    Q_past[1].vector[j] = data(11+j);
  }
  Q_past[0].scalar = data(10);
  Q_past[1].scalar = data(15);

  for (int k=0; k<3; k++) {
    xAxis(k) = data(15+k);
    vz[k] = data(18+k);
    nodeIOffset(k) = data(21+k);
    nodeJOffset(k) = data(24+k);
    alphaI(k) = data(27+k);
    alphaJ(k) = data(30+k);
  }

  int flag;
  flag = 0;
  for (int i=34; i<=39; i++)
    if (data(i) != 0.0)
      flag = 1;
  if (flag == 1) {
    if (nodeIInitialDisp == 0)
      nodeIInitialDisp = new double[6];
    for (int i=34, j = 0; i<=39; i++, j++)
      nodeIInitialDisp[j] = data(i);
  }

  flag = 0;
  for (int i=40; i<=45; i++)
    if (data(i) != 0.0)
      flag = 1;

  if (flag == 1) {
    if (nodeJInitialDisp == 0)
      nodeJInitialDisp = new double[6];
    for (int i=40, j = 0; i<=45; i++, j++)
      nodeJInitialDisp[j] = data(i);
  }
  L  = data(46);
  Ln = data(47);

  ul = ulcommit;
  Q_pres[0] = Q_past[0];
  Q_pres[1] = Q_past[1];

  initialDispChecked = true;
  return 0;
}

int
CorotFrameTransf3d03::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(48);
  for (int i = 0; i<7; i++)
    data[i] = ulcommit[i];

  // [7 - 15]: 4 times 2 quaternions
  for (int j = 0; j<3; j++) {
    data(7+j)  = Q_past[0].vector[j];
    data(11+j) = Q_past[1].vector[j];
  }
  data(10) = Q_past[0].scalar;
  data(15) = Q_past[1].scalar;

  for (int k=0; k<3; k++) {
    data(15+k) = xAxis(k);
    data(18+k) = vz[k];
    data(21+k) = nodeIOffset(k);
    data(24+k) = nodeJOffset(k);
    data(27+k) = alphaI(k);
    data(30+k) = alphaJ(k);
  }

  if (nodeIInitialDisp != 0) {
    data(34) = nodeIInitialDisp[0];
    data(35) = nodeIInitialDisp[1];
    data(36) = nodeIInitialDisp[2];
    data(37) = nodeIInitialDisp[3];
    data(38) = nodeIInitialDisp[4];
    data(39) = nodeIInitialDisp[5];
  } else {
    data(34)  = 0.0;
    data(35)  = 0.0;
    data(36) = 0.0;
    data(37) = 0.0;
    data(38) = 0.0;
    data(39) = 0.0;
  }

  if (nodeJInitialDisp != 0) {
    data(40) = nodeJInitialDisp[0];
    data(41) = nodeJInitialDisp[1];
    data(42) = nodeJInitialDisp[2];
    data(43) = nodeJInitialDisp[3];
    data(44) = nodeJInitialDisp[4];
    data(45) = nodeJInitialDisp[5];
  } else {
    data(40) = 0.0;
    data(41) = 0.0;
    data(42) = 0.0;
    data(43) = 0.0;
    data(44) = 0.0;
    data(45) = 0.0;
  }
  data(46) = L;
  data(47) = Ln;

  if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
    opserr << " CorotFrameTransf3d03::sendSelf() - data could not be sent\n" ;
    return -1;
  }

  return 0;
}

