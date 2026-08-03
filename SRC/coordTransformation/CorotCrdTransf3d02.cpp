//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//    See https://opensees.berkeley.edu/OpenSees/copyright.php for license.
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the implementation for the
// CorotCrdTransf3d02 class. CorotCrdTransf3d02 is a Corotational
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
#include <math.h>
#include <Triad.h>
#include <Vector.h>
#include <Matrix.h>
#include <VectorND.h>
#include <Vector3D.h>
#include <MatrixND.h>
#include <Matrix3D.h>
using namespace OpenSees;
#include <Node.h>
#include <Channel.h>
#include <Logging.h>
#include <CorotCrdTransf3d02.h>

#include "RotationLibrary.h"

#undef OPS_STATIC
#define OPS_STATIC // static 
// initialize static variables
Matrix CorotCrdTransf3d02::Tp(6,7);
Matrix CorotCrdTransf3d02::Tbl(6,12);
Matrix CorotCrdTransf3d02::kg(12,12);
Matrix CorotCrdTransf3d02::Lr2(12,3);
Matrix CorotCrdTransf3d02::Lr3(12,3);


// Permutation matrix (to renumber basic dof's)

// v = Tp * ul
//
//      |  thI  |   thJ  |
//      |0 1  2 | 3 4  5 | 6
// Tp=  [0 0  0   0 0  0   1;  0 // axial
//       0 1  0   0 0  0   0;  1 // rot z I
//       0 0  0   0 1  0   0;  2 // rot z J
//       0 0 -1   0 0  0   0;  3 // rot y I
//       0 0  0   0 0 -1   0;  4 // rot y J
//      -1 0  0   1 0  0   0]; 5 // torsion
// constexpr MatrixND<6,7> T = {{
//     {0,    0,    0,    0,    0,   -1 },
//     {0,    1,    0,    0,    0,    0 },
//     {0,    0,    0,   -1,    0,    0 },
//     {0,    0,    0,    0,    0,    1 },
//     {0,    0,    1,    0,    0,    0 },
//     {0,    0,    0,    0,   -1,    0 },
//     {1,    0,    0,    0,    0,    0 }}};


static inline void
getLMatrix(const Matrix3D& A, const Vector3D& e1, const Vector3D& r1, const Vector3D &ri, Matrix& L)
{
  OPS_STATIC Matrix3D L1, L2;
  OPS_STATIC Matrix3D rie1r1;
  OPS_STATIC Matrix3D e1e1r1;

  const double rie1 = ri.dot(e1);

  for (int k = 0; k < 3; k++) {
    const double e1r1k = (e1[k] + r1[k]);
    for (int j = 0; j < 3; j++) {
      rie1r1(j,k) = ri[j]*e1r1k;
      e1e1r1(j,k) = e1[j]*e1r1k;
    }
  }

  //L1  = ri'*e1 * A/2 + A*ri*(e1 + r1)'/2;
  L1.zero();
  L1.addMatrix(A, rie1*0.5);
  L1.addMatrixProduct(A, rie1r1, 0.5);

  // L2  = Sri/2 - ri'*e1*S(r1)/4 - Sri*e1*(e1 + r1)'/4;
  L2.zero();
  L2.addSpin(ri, 0.5);
  L2.addSpin(r1, -rie1/4.0);
  L2.addSpinMatrixProduct(ri, e1e1r1, -0.25);

  // L = [L1
  //      L2
  //     -L1
  //      L2];

  L.Zero();
  L.Assemble(L1, 0, 0,  1.0);
  L.Assemble(L2, 3, 0,  1.0);
  L.Assemble(L1, 6, 0, -1.0);
  L.Assemble(L2, 9, 0,  1.0);

}

static inline const MatrixND<12,12> &
// static inline const Matrix &
getKs2Matrix(Matrix3D& A, const Vector3D& e1, const Vector3D& r1, const double Ln, const Vector3D &ri, const Vector3D &z)
{
    static MatrixND<12,12> ks2;
//  static Matrix ks2(12,12);

    //  Ksigma2 = [ K11   K12 -K11   K12;
    //              K12t  K22 -K12t  K22;
    //             -K11  -K12  K11  -K12;
    //              K12t  K22 -K12t  K22];

    // U = (-1/2)*A*z*ri'*A + ri'*e1*A*z*e1'/(2*Ln)+...
    //      z'*(e1+r1)*A*ri*e1'/(2*Ln);

    double rite1 = 0;   // dot product ri . e1
    double zte1  = 0;   // dot product z  . e1
    double ztr1  = 0;   // dot product z  . r1

    for (int i = 0; i < 3; i++) {
      rite1 += ri(i)*e1[i];
      zte1  += z(i)*e1[i];
      ztr1  += z(i)*r1[i];
    }

    OPS_STATIC Matrix zrit(3,3), ze1t(3,3);
    OPS_STATIC Matrix rizt(3,3), rie1t(3,3);
    OPS_STATIC Matrix3D e1zt; // r1e1t(3,3),

//  const Matrix3D e1zt = e1.bun(z);

    // Chrystal's bun order
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 3; i++) {
        zrit(i,j)  = z(i)*ri(j);
        rizt(i,j)  = ri(i)*z(j);
        ze1t(i,j)  = z(i)*e1(j);
        e1zt(i,j)  = e1[i]*z(j);
        // r1e1t(i,j) = r1(i)*e1(j);
        rie1t(i,j) = ri(i)*e1(j);
      }
    }

    OPS_STATIC Matrix U(3,3);

    U.addMatrixTripleProduct(0.0, A, zrit, -0.5);
    U.addMatrixProduct(1.0, A, ze1t,   rite1/(2*Ln));
    U.addMatrixProduct(1.0, A, rie1t, (zte1 + ztr1)/(2*Ln));

    OPS_STATIC Matrix3D ks;
    OPS_STATIC Matrix3D m1;

    //K11 = U + U' + ri'*e1*(2*(e1'*z)+z'*r1)*A/(2*Ln);
    ks.zero();
    ks.addMatrix(U, 1.0);

    // add matrix U transpose
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        ks(i,j) += U(j,i);

    ks.addMatrix(A, rite1*(2*zte1 + ztr1)/(2*Ln));

    ks2.zero();
    ks2.assemble(ks, 0, 0,  1.0);
    ks2.assemble(ks, 0, 6, -1.0);
    ks2.assemble(ks, 6, 0, -1.0);
    ks2.assemble(ks, 6, 6,  1.0);

    //K12 = (1/4)*(-A*z*e1'*Sri - A*ri*z'*Sr1 - z'*(e1+r1)*A*Sri);
    m1.zero();
    m1.addMatrixProduct(A, ze1t, -1.0);
    ks.zero();
    ks.addMatrixSpinProduct(m1, ri, 0.25);

    m1.zero();
    m1.addMatrixProduct(A, rizt, -1.0);
    ks.addMatrixSpinProduct(m1, r1, 0.25);
    ks.addMatrixSpinProduct(A, ri, -0.25*(zte1+ztr1));

    ks2.assemble(ks, 0, 3,  1.0);
    ks2.assemble(ks, 0, 9,  1.0);
    ks2.assemble(ks, 6, 3, -1.0);
    ks2.assemble(ks, 6, 9, -1.0);

    ks2.assembleTranspose(ks, 3, 0,  1.0);
    ks2.assembleTranspose(ks, 3, 6, -1.0);
    ks2.assembleTranspose(ks, 9, 0,  1.0);
    ks2.assembleTranspose(ks, 9, 6, -1.0);

    // K22 = (1/8)*((-ri'*e1)*Sz*Sr1 + Sr1*z*e1'*Sri + ...
    //       Sri*e1*z'*Sr1 - (e1+r1)'*z*S(e1)*Sri + 2*Sz*Sri);

    ks.zero();
    ks.addSpinProduct(z, r1, -0.125*(rite1));

    m1.zero();
    m1.addSpinMatrixProduct( r1, ze1t, 1.0);
    ks.addMatrixSpinProduct( m1, ri, 0.125);

    m1.zero();
    m1.addSpinMatrixProduct(ri, e1zt, 1.0);
    ks.addMatrixSpinProduct(m1, r1, 0.125);

    ks.addSpinProduct(e1, ri, -0.125*(zte1 + ztr1));
    ks.addSpinProduct( z, ri, 0.25);

    // Ksigma2 = [ K11   K12 -K11   K12;
    //             K12t  K22 -K12t  K22;
    //            -K11  -K12  K11  -K12;
    //             K12t  K22 -K12t  K22];

    ks2.assemble(ks, 3, 3, 1.0);
    ks2.assemble(ks, 3, 9, 1.0);
    ks2.assemble(ks, 9, 3, 1.0);
    ks2.assemble(ks, 9, 9, 1.0);

    return ks2;
}


// constructor:
CorotCrdTransf3d02::CorotCrdTransf3d02(int tag, const Vector &vecInLocXZPlane,
                                   const Vector &rigJntOffsetI,
                                   const Vector &rigJntOffsetJ):
  CrdTransf(tag, CRDTR_TAG_CorotCrdTransf3d02),
  vAxis(3), nodeIOffset(3), nodeJOffset(3), xAxis(3),
  nodeIPtr(0), nodeJPtr(0), L(0), Ln(0),
  // alphaIq(4), alphaJq(4), alphaIqcommit(4), alphaJqcommit(4), 
  alphaI(3), alphaJ(3),
  ulcommit(7), ul(7),  ulpr(7), T(7,12),
  nodeIInitialDisp(0), nodeJInitialDisp(0), initialDispChecked(false)
{
    // check vector that defines local xz plane
    if (vecInLocXZPlane.Size() != 3 ) {
        opserr << "CorotCrdTransf3d02::CorotCrdTransf3d02:  Vector that defines local xz plane is invalid\n";
        opserr << "Size must be 3\n. Using (0,0,1)";
        vAxis(0) = 0;
        vAxis(1) = 0;
        vAxis(2) = 1;
    }
    else
        vAxis = vecInLocXZPlane;

    // check rigid joint offset for node I
    if (rigJntOffsetI.Size() != 3 ) {
        opserr << "CorotCrdTransf3d02::CorotCrdTransf3d02:  Invalid rigid joint offset vector for node I\n";
        opserr << "Size must be 3\n";
        nodeIOffset.Zero();
    }
    else
        nodeIOffset = rigJntOffsetI;

    // check rigid joint offset for node J
    if (rigJntOffsetJ.Size() != 3 ) {
        opserr << "CorotCrdTransf3d02::CorotCrdTransf3d02:  Invalid rigid joint offset vector for node J\n";
        opserr << "Size must be 3\n";
        nodeJOffset.Zero();
    }
    else
        nodeJOffset = rigJntOffsetJ;

    // temporary
    if (nodeIOffset.Norm() != 0 || nodeJOffset.Norm() != 0) {
        opserr << "CorotCrdTransf3d02::CorotCrdTransf3d02: rigid joint zones not implemented yet\n";
        opserr << "Using zero values\n";
        nodeIOffset.Zero();
        nodeJOffset.Zero();
    }

    // Permutation matrix (to renumber basic dof's)

    // v = Tp * ul
    //
    //       0 1  2 | 3 4  5 | 6
    // Tp=  [0 0  0   0 0  0   1;  0 // axial
    //       0 1  0   0 0  0   0;  1 // rot z I
    //       0 0  0   0 1  0   0;  2 // rot z J
    //       0 0 -1   0 0  0   0;  3 // rot y I
    //       0 0  0   0 0 -1   0;  4 // rot y J
    //      -1 0  0   1 0  0   0]; 5 // torsion
    //

    // using static matrix (one constant matrix for all objects)

    if (Tp(0, 6) == 0) {
      // initialize only once
        Tp(0, 6) =  1;
        Tp(1, 1) =  1;
        Tp(2, 4) =  1;
        Tp(3, 2) = -1;
        Tp(4, 5) = -1;
        Tp(5, 0) = -1;
        Tp(5, 3) =  1;
    }
}


// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
CorotCrdTransf3d02::CorotCrdTransf3d02():
  CrdTransf(0, CRDTR_TAG_CorotCrdTransf3d02),
  vAxis(3), nodeIOffset(3), nodeJOffset(3), xAxis(3),
  nodeIPtr(0), nodeJPtr(0), L(0), Ln(0),
//  alphaIq(4), alphaJq(4),  alphaIqcommit(4), alphaJqcommit(4), 
  alphaI(3), alphaJ(3),
  ulcommit(7), ul(7),  ulpr(7), T(7,12),
  nodeIInitialDisp(0), nodeJInitialDisp(0), initialDispChecked(false)
{
    // Permutation matrix (to renumber basic dof's)

    //       0 1  2 3 4  5 6
    //
    // Tp=  [0 0  0 0 0  0 1;  0
    //       0 1  0 0 0  0 0;  1
    //       0 0  0 0 1  0 0;  2
    //       0 0 -1 0 0  0 0;  3
    //       0 0  0 0 0 -1 0;  4
    //      -1 0  0 1 0  0 0]; 5

    // using static matrix (one constant matrix for all objects)

    if (Tp(0, 6) == 0) {
      // initialize only once
      Tp(0, 6) =  1;
      Tp(1, 1) =  1;
      Tp(2, 4) =  1;
      Tp(3, 2) = -1;
      Tp(4, 5) = -1;
      Tp(5, 0) = -1;
      Tp(5, 3) =  1;
    }
}


// destructor:
CorotCrdTransf3d02::~CorotCrdTransf3d02()
{
    if (nodeIInitialDisp != nullptr)
        delete [] nodeIInitialDisp;

    if (nodeJInitialDisp != nullptr)
        delete [] nodeJInitialDisp;
}


int
CorotCrdTransf3d02::commitState()
{
    ulcommit      = ul;
    alphaIqcommit = alphaIq;
    alphaJqcommit = alphaJq;

    return 0;
}


int
CorotCrdTransf3d02::revertToLastCommit()
{
    // determine global displacement increments from last iteration
    const Vector &dispI = nodeIPtr->getTrialDisp();
    const Vector &dispJ = nodeJPtr->getTrialDisp();

    for (int k = 0; k < 3; k++) {
      alphaI(k) =  dispI(k+3);
      alphaJ(k) =  dispJ(k+3);
    }

    if (nodeIInitialDisp != 0) {
      for (int j = 0; j<3; j++)
        alphaI(j) -= nodeIInitialDisp[j+3];
    }

    if (nodeJInitialDisp != 0) {
      for (int j = 0; j<3; j++)
        alphaJ(j) -= nodeJInitialDisp[j+3];
    }

    ul      = ulcommit;
    alphaIq = alphaIqcommit;
    alphaJq = alphaJqcommit;

    this->update();

    return 0;
}


int
CorotCrdTransf3d02::revertToStart()
{
  ul.Zero();
  alphaIq = VersorFromMatrix(R0);    // pseudo-vector for node 1
  alphaJq = VersorFromMatrix(R0);    // pseudo-vector for node J

  alphaI.Zero();
  alphaJ.Zero();

  this->update();
  return 0;
}


int
CorotCrdTransf3d02::initialize(Node *nodeIPointer, Node *nodeJPointer)
{
    int error;

    nodeIPtr = nodeIPointer;
    nodeJPtr = nodeJPointer;

    if ((!nodeIPtr) || (!nodeJPtr)) {
      opserr << "\nCorotCrdTransf3d02::initialize";
      opserr << "\ninvalid pointers to the element nodes\n";
      return -1;
    }

    // see if there is some initial displacements at nodes
    if (initialDispChecked == false) {
      const Vector &nodeIDisp = nodeIPtr->getDisp();
      const Vector &nodeJDisp = nodeJPtr->getDisp();
      for (int i = 0; i<6; i++)
      if (nodeIDisp(i) != 0.0) {
        nodeIInitialDisp = new double [6];
        for (int j = 0; j<6; j++)
          nodeIInitialDisp[j] = nodeIDisp(j);
        i = 6;
      }

      for (int j = 0; j<6; j++)
        if (nodeJDisp(j) != 0.0) {
          nodeJInitialDisp = new double [6];
          for (int i = 0; i<6; i++)
            nodeJInitialDisp[i] = nodeJDisp(i);
          j = 6;
        }
      initialDispChecked = true;
    }

    static Vector XAxis(3);
    static Vector YAxis(3);
    static Vector ZAxis(3);

    // get 3by3 rotation matrix
    if ((error = this->getLocalAxes(XAxis, YAxis, ZAxis)))
      return error;

    // compute initial pseudo-vectors for nodal triads
    alphaIq = VersorFromMatrix(R0); // pseudo-vector for node I
    alphaJq = VersorFromMatrix(R0); // pseudo-vector for node J

    this->commitState();

    return 0;
}

void inline
CorotCrdTransf3d02::compTransfMatrixBasicGlobal(const Triad& __restrict r, 
                                                const Triad& __restrict E, 
                                                const Triad& __restrict rI, 
                                                const Triad& __restrict rJ)
{

    // extract columns of rotation matrices
    const Vector3D &r1 = r[1], &r2 = r[2], &r3 = r[3],
                   &e1 = E[1], &e2 = E[2], &e3 = E[3],
                   &rI1=rI[1], &rI2=rI[2], &rI3=rI[3],
                   &rJ1=rJ[1], &rJ2=rJ[2], &rJ3=rJ[3];

    // compute the transformation matrix from the basic to the
    // global system
    //   A = (1/Ln)*(I - e1*e1');
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        A(i,j) = (double(i==j) - e1[i]*e1[j])/Ln;

    // This must be called up here
    getLMatrix(A, e1, r1, r2, Lr2);
    getLMatrix(A, e1, r1, r3, Lr3);

    //   T1 = [      O', (-S(rI3)*e2 + S(rI2)*e3)',        O', O']';
    //   T2 = [(A*rI2)', (-S(rI2)*e1 + S(rI1)*e2)', -(A*rI2)', O']';
    //   T3 = [(A*rI3)', (-S(rI3)*e1 + S(rI1)*e3)', -(A*rI3)', O']';
    //
    //   T4 = [      O', O',        O', (-S(rJ3)*e2 + S(rJ2)*e3)']';
    //   T5 = [(A*rJ2)', O', -(A*rJ2)', (-S(rJ2)*e1 + S(rJ1)*e2)']';
    //   T6 = [(A*rJ3)', O', -(A*rJ3)', (-S(rJ3)*e1 + S(rJ1)*e3)']';

    T.Zero();

    //   T1 = [      O', (-S(rI3)*e2 + S(rI2)*e3)',        O', O']';

    // (-S(rI3)*e2 + S(rI2)*e3)
    Vector3D Se  = rI2.cross(e3);
    Se -= rI3.cross(e2);

    for (int i = 0; i < 3; i++)
      T(0,i+3) =  Se(i);

    //   T2 = [(A*rI2)', (-S(rI2)*e1 + S(rI1)*e2)', -(A*rI2)', O']';

    Vector3D At = A*rI2;

    // (-S(rI2)*e1 + S(rI1)*e2)'
    Se  = rI1.cross(e2);
    Se -= rI2.cross(e1);

    for (int i = 0; i < 3; i++) {
        T(1,i  ) =  At[i];
        T(1,i+3) =  Se[i];
        T(1,i+6) = -At[i];
    }

    //   T3 = [(A*rI3)', (-S(rI3)*e1 + S(rI1)*e3)', -(A*rI3)', O']';

    At = A*rI3;
//  At.addMatrixVector(0.0, A, rI3, 1.0);
    
    // (-S(rI3)*e1 + S(rI1)*e3)
    Se  = rI1.cross(e3);
    Se -= rI3.cross(e1);

    for (int i = 0; i < 3; i++) {
        T(2,i  ) =  At(i);
        T(2,i+3) =  Se(i);
        T(2,i+6) = -At(i);
    }

    //   T4 = [      O', O',        O', (-S(rJ3)*e2 + S(rJ2)*e3)']';
    Se  = rJ2.cross(e3);  // -S(rJ3)*e2 + S(rJ2)*e3
    Se -= rJ3.cross(e2);
    for (int i = 0; i < 3; i++)
      T(3, i+9) =  Se(i);

    // T5 = [(A*rJ2)', O', -(A*rJ2)', (-S(rJ2)*e1 + S(rJ1)*e2)']';
    At = A*rJ2;
    Se  = rJ1.cross(e2); // (-S(rJ2)*e1 + S(rJ1)*e2)
    Se -= rJ2.cross(e1);
    for (int i = 0; i < 3; i++) {
        T(4, i  ) =  At[i];
        T(4, i+6) = -At[i];
        T(4, i+9) =  Se[i];
    }

    // T6 = [(A*rJ3)', O', -(A*rJ3)', (-S(rJ3)*e1 + S(rJ1)*e3)']';
    At  = A*rJ3;  
    Se  = rJ1.cross(e3);  // (-S(rJ3)*e1 + S(rJ1)*e3)
    Se -= rJ3.cross(e1);
    for (int i = 0; i < 3; i++) {
        T(5,i  ) =  At[i];
        T(5,i+6) = -At[i];
        T(5,i+9) =  Se[i];
    }

    // setup tranformation matrix

    // T(:,1) += Lr3*rI2 - Lr2*rI3;
    // T(:,2) +=           Lr2*rI1;
    // T(:,3) += Lr3*rI1          ;

    // T(:,4) += Lr3*rJ2 - Lr2*rJ3;
    // T(:,5) += Lr2*rJ1          ;      // ?????? check sign
    // T(:,6) += Lr3*rJ1          ;      // ?????? check sign

    for (int i = 0; i < 12; i++) {
      double T1i = 0;
      for (int k=0; k<3; k++)
        T1i += Lr2(i,k)*rI1[k];
      T(1,i) += T1i;
    }

    for (int i = 0; i < 12; i++) {
      double T4i = 0;
      for (int k=0; k<3; k++)
        T4i += Lr2(i,k)*rJ1[k]; // Lr(i);
      T(4,i) += T4i;
    }

    for (int i = 0; i < 12; i++) {
      double T0i = 0;
      for (int k=0; k<3; k++)
        T0i += Lr3(i,k)*rI2[k] - Lr2(i,k)*rI3[k];
      T(0,i) += T0i;
    }

    for (int i = 0; i < 12; i++) {
      double T3i = 0;
      for (int k=0; k<3; k++)
        T3i += Lr3(i,k)*rJ2[k] - Lr2(i,k)*rJ3[k];
      T(3,i) += T3i;
    }

    for (int i = 0; i < 12; i++) {
      double T2i = 0;
      for (int k=0; k<3; k++)
        T2i += Lr3(i,k)*rI1[k]; // Lr(i);
      T(2,i) += T2i;
    }

    for (int i = 0; i < 12; i++) {
      double T5i = 0;
      for (int k=0; k<3; k++)
        T5i += Lr3(i,k)*rJ1[k]; // Lr(i);
      T(5,i) += T5i;
    }

    for (int j = 0; j < 6; j++) {
//      const double c = 2 * cos(ul(j));
        const double c = 0.5 / cos(ul(j));
//      const double c = 0.5 / std::sqrt(1-sn[j]*sn[j]);
        for (int i = 0; i < 12; i++)
          T(j,i) *= c;
//        T(j,i) /= c;
    }

    // T(:,7) = [-e1' O' e1' O']';
    for (int i = 0; i < 3; i++) {
        T(6,i  ) = -e1[i];
        T(6,i+6) =  e1[i];
    }
}


// Set RI,RJ,Rbar, Ln, e and ul
int
CorotCrdTransf3d02::update()
{
    /********* OLD REMO - REPLACED BELOW TO FIX BUG ***************/
#if 0
    // determine global displacement increments from last iteration
    const Vector &dispIncrI = nodeIPtr->getIncrDeltaDisp();
    const Vector &dispIncrJ = nodeJPtr->getIncrDeltaDisp();

     // get the iterative spins dAlphaI and dAlphaJ
     // (rotational displacement increments at both nodes)

      static Vector dAlphaI(3);
      static Vector dAlphaJ(3);

      for (int k = 0; k < 3; k++) {
        dAlphaI(k) = dispIncrI(k+3);
        dAlphaJ(k) = dispIncrJ(k+3);
      }
#endif
    /**************************************************************/

    // determine global displacement increments from last iteration
    static Vector dispI(6);
    static Vector dispJ(6);
    dispI = nodeIPtr->getTrialDisp();
    dispJ = nodeJPtr->getTrialDisp();

    if (nodeIInitialDisp != 0) {
      for (int j = 0; j<6; j++)
        dispI(j) -= nodeIInitialDisp[j];
    }

    if (nodeJInitialDisp != 0) {
      for (int j = 0; j<6; j++)
        dispJ(j) -= nodeJInitialDisp[j];
    }

    // get the iterative spins dAlphaI and dAlphaJ
    // (rotational displacement increments at both nodes)
    {
      OPS_STATIC Vector3D dAlphaI, dAlphaJ;

      for (int k = 0; k < 3; k++) {
        dAlphaI[k] =  dispI(k+3) - alphaI[k];
        dAlphaJ[k] =  dispJ(k+3) - alphaJ[k];
        alphaI[k]  =  dispI(k+3);
        alphaJ[k]  =  dispJ(k+3);
      }

      // update the nodal triads TI and RJ using quaternions

      const VectorND<4> dAlphaIq = getQuaternionFromPseudoRotVector(dAlphaI);
      const VectorND<4> dAlphaJq = getQuaternionFromPseudoRotVector(dAlphaJ);

      alphaIq = quaternionProduct(alphaIq, dAlphaIq);
      alphaJq = quaternionProduct(alphaJq, dAlphaJq);

      getRotationMatrixFromQuaternion(alphaIq, this->RI);
      getRotationMatrixFromQuaternion(alphaJq, this->RJ);
    }

    //
    // compute the mean nodal triad
    //
    {
      Matrix3D dRgamma;
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
          dRgamma(i,j) = 0.0;
          for (int k = 0; k < 3; k++)
            dRgamma(i,j) += RJ(i,k) * RI(j,k);
        }

      // pseudo-vector for node J
      VectorND<4> gammaq = VersorFromMatrix(dRgamma);

      OPS_STATIC Vector3D gammaw;
      getTangScaledPseudoVectorFromQuaternion(gammaq, gammaw);

      gammaw *= 0.5;

      dRgamma = CaySO3(gammaw);

      Rbar.zero();
      Rbar.addMatrixProduct(dRgamma, RI, 1.0);
    }
    const Triad r(Rbar);

    // -----------------------------------------------
    // Node offsets
    // -----------------------------------------------

    // relative translation displacements
    OPS_STATIC Vector3D dJI;
    for (int k = 0; k < 3; k++)
      dJI[k] = dispJ(k) - dispI(k);

    // element projection
    OPS_STATIC Vector3D xJI;
    xJI  = nodeJPtr->getCrds() - nodeIPtr->getCrds();

    if (nodeIInitialDisp != nullptr) {
      xJI[0] -= nodeIInitialDisp[0];
      xJI[1] -= nodeIInitialDisp[1];
      xJI[2] -= nodeIInitialDisp[2];
    }

    if (nodeJInitialDisp != nullptr) {
      xJI[0] += nodeJInitialDisp[0];
      xJI[1] += nodeJInitialDisp[1];
      xJI[2] += nodeJInitialDisp[2];
    }

    OPS_STATIC Vector3D dx;
    // dx = xJI + dJI;
    dx  = xJI;
    dx += dJI;

    // calculate the deformed element length
    Ln = dx.norm();

    if (Ln == 0.0) {
      opserr << "\nCorotCrdTransf3d02::update: 0 deformed length\n";
      return -2;
    }

    // -----------------------------------------------
    // compute the base vectors e1, e2, e3
    // -----------------------------------------------
    const Triad rI{RI}, rJ{RJ};

    OPS_STATIC Vector3D e1, e2, e3;

    e1  = dx;
    e1 /= Ln;

    // 'rotate' the mean rotation matrix Rbar on to e1 to
    // obtain e2 and e3 (using the 'mid-point' procedure)

    // e2 = r2 - (e1 + r1)*((r2^e1)*0.5);
    // e3 = r3 - (e1 + r1)*((r3^e1)*0.5);

    OPS_STATIC Vector3D tmp;
    tmp  = e1;
    tmp += r[1];

    e2 = tmp;
    e3 = tmp;

    e2 *= 0.5*r[2].dot(e1);
    e2.addVector(-1.0,  r[2], 1.0);

    // e2 = r2 - (e1 + r1)*((r2^e1)*0.5);

    e3 *= r[3].dot(e1)*0.5;
    e3.addVector(-1.0,  r[3], 1.0);

    for (int k = 0; k < 3; k ++) {
      e(k,0) = e1[k];
      e(k,1) = e2[k];
      e(k,2) = e3[k];
    }

    // -----------------------------------------------
    // compute the basic deformations
    // -----------------------------------------------

    // save previous deformations
    ulpr = ul;

    // Rotations
    {
      Vector3D thetaI = LogC90(RI^e);
      for (int i=0; i<3; i++)
        ul(i) = thetaI[i];

      thetaI = LogC90(RJ^e);
      for (int i=3; i<6; i++)
        ul(i) = thetaI[i-3];
    }

    // Axial
    // ul = Ln - L;
    // ul(6) = 2 * ((xJI + dJI/2)^ dJI) / (Ln + L);  //mid-point formula
//  xJI.addVector(1.0, dJI, 0.5);
    ul(6) = Ln - L; // 2. * xJI.dot(dJI) / (Ln + L);  //mid-point formula

    // compute the transformation matrix
    this->compTransfMatrixBasicGlobal(r, Triad{e}, rI, rJ);

    return 0;
}


void
CorotCrdTransf3d02::compTransfMatrixLocalGlobal(Matrix &Tlg)
{
    // setup transformation matrix from local to global
    Tlg.Zero();

    Tlg(0,0) = Tlg(3,3) = Tlg(6,6) = Tlg(9,9)   = R0(0,0);
    Tlg(0,1) = Tlg(3,4) = Tlg(6,7) = Tlg(9,10)  = R0(1,0);
    Tlg(0,2) = Tlg(3,5) = Tlg(6,8) = Tlg(9,11)  = R0(2,0);
    Tlg(1,0) = Tlg(4,3) = Tlg(7,6) = Tlg(10,9)  = R0(0,1);
    Tlg(1,1) = Tlg(4,4) = Tlg(7,7) = Tlg(10,10) = R0(1,1);
    Tlg(1,2) = Tlg(4,5) = Tlg(7,8) = Tlg(10,11) = R0(2,1);
    Tlg(2,0) = Tlg(5,3) = Tlg(8,6) = Tlg(11,9)  = R0(0,2);
    Tlg(2,1) = Tlg(5,4) = Tlg(8,7) = Tlg(11,10) = R0(1,2);
    Tlg(2,2) = Tlg(5,5) = Tlg(8,8) = Tlg(11,11) = R0(2,2);
}


void
CorotCrdTransf3d02::compTransfMatrixBasicLocal(Matrix &Tbl)
{
  // setup transformation matrix from basic to local
  Tbl.Zero();

  // first get transformation matrix from basic to global
  static Matrix Tbg(6, 12);
  Tbg.addMatrixProduct(0.0, Tp, T, 1.0);

  // get inverse of transformation matrix from local to global
  static Matrix Tlg(12,12), TlgInv(12,12);
  this->compTransfMatrixLocalGlobal(Tlg);
  TlgInv.addMatrixTranspose(0.0, Tlg, 1.0);  // for square rot-matrix: Tlg^-1 = Tlg'

  // finally get transformation matrix from basic to local
  Tbl.addMatrixProduct(0.0, Tbg, TlgInv, 1.0);
}


const Vector &
CorotCrdTransf3d02::getBasicTrialDisp()
{
    static Vector ub(6);

    // use transformation matrix to renumber the degrees of freedom
    ub.addMatrixVector(0.0, Tp, ul, 1.0);

    return ub;
}


const Vector &
CorotCrdTransf3d02::getBasicIncrDeltaDisp()
{
    static Vector dub(6);
    static Vector dul(7);

    // dul = ul - ulpr;
    dul = ul;
    dul.addVector(1.0, ulpr, -1.0);

    // use transformation matrix to renumber the degrees of freedom
    dub.addMatrixVector(0.0, Tp, dul, 1.0);

    return dub;
}


const Vector &
CorotCrdTransf3d02::getBasicIncrDisp()
{
    static Vector Dub(6);
    static Vector Dul(7);

    // Dul = ul - ulcommit;
    Dul = ul;
    Dul.addVector (1.0, ulcommit, -1.0);

    // use transformation matrix to renumber the degrees of freedom
    Dub.addMatrixVector(0.0, Tp, Dul, 1.0);

    return Dub;
}


const Vector &
CorotCrdTransf3d02::getBasicTrialVel()
{
  opserr << "WARNING CorotCrdTransf3d02::getBasicTrialVel()"
      << " - has not been implemented yet. Returning zeros." << endln;

  static Vector dummy(6);
  return dummy;
}


const Vector &
CorotCrdTransf3d02::getBasicTrialAccel()
{
  opserr << "WARNING CorotCrdTransf3d02::getBasicTrialAccel()"
      << " - has not been implemented yet. Returning zeros." << endln;

  static Vector dummy(6);
  return dummy;
}


const Vector &
CorotCrdTransf3d02::getGlobalResistingForce(const Vector &pb, const Vector &p0)
{
    //cmp// this->update();

    static Vector pg(12);

    // if there are no element loads present
    if (p0 == 0.0) {
        // transform resisting forces from the basic system to local coordinates
        static Vector pl(7);
        pl.addMatrixTransposeVector(0.0, Tp, pb, 1.0);    // pl = Tp ^ pb;

        // transform resisting forces from local to global coordinates
        pg.addMatrixTransposeVector(0.0, T, pl, 1.0);   // pg = T ^ pl; residual
    }

    // if there are element loads present
    else {
#if 0
        // SLOWER!!!! THIS IS SAME APPROACH AS 2D CASE
        // ===========================================
        // transform resisting forces from the basic system to local coordinates
        this->compTransfMatrixBasicLocal(Tbl);
        static Vector pl(12);
        pl.addMatrixTransposeVector(0.0, Tbl, pb, 1.0);    // pl = Tbl ^ pb;

        // add end forces due to element p0 loads
        // assuming member loads are in local system
        pl(0) += p0(0);
        pl(1) += p0(1);
        pl(7) += p0(2);
        pl(2) += p0(3);
        pl(8) += p0(4);

        // transform resisting forces from local to global coordinates
        pg.addMatrixTransposeVector(0.0, Tlg, pl, 1.0);   // pg = Tlg ^ pl;
#endif
        // FASTER!!!! TRANSFORM REACTIONS AND ADD AT END
        // =============================================
        // transform resisting forces from the basic system to local coordinates
        static Vector pl(7);
        pl.addMatrixTransposeVector(0.0, Tp, pb, 1.0);    // pl = Tp ^ pb;

        // transform resisting forces from local to global coordinates
        pg.addMatrixTransposeVector(0.0, T, pl, 1.0);   // pg = T ^ pl; residual

        // add end forces due to element p0 loads
        // assuming member loads are in local system
        // TODO(cmp): Do this better!!! 
        static Vector pl0(12), pg0(12);
        pl0.Zero();
        pl0(0) = p0(0);
        pl0(1) = p0(1);
        pl0(7) = p0(2);
        pl0(2) = p0(3);
        pl0(8) = p0(4);
        static Matrix Tlg(12,12);
        this->compTransfMatrixLocalGlobal(Tlg);
        pg0.addMatrixTransposeVector(0.0, Tlg, pl0, 1.0);
        pg.addVector(1.0, pg0, 1.0);
    }

    return pg;
}

const Matrix &
CorotCrdTransf3d02::getGlobalStiffMatrix(const Matrix &kb, const Vector &pb)
{
    const Triad r{Rbar}, rI{RI}, rJ{RJ}, E{e};
    const Vector3D &r1 = r[1], &r2 = r[2], &r3 = r[3],
                   &e1 = E[1], &e2 = E[2], &e3 = E[3],
                   &rI1=rI[1], &rI2=rI[2], &rI3=rI[3],
                   &rJ1=rJ[1], &rJ2=rJ[2], &rJ3=rJ[3];

#if 0
// NOTE[cmp] CorotCrdTransf3d02::compTransfMatrixBasicGlobal must be called first
// to set Lr1, Lr2 and T
    this->update();
#else
    getLMatrix(A, e1, r1, r2, Lr2);
    getLMatrix(A, e1, r1, r3, Lr3);
#endif

    // transform tangent stiffness matrix from the basic system to local coordinates
    static Matrix kl(7,7);
    kl.addMatrixTripleProduct(0.0, Tp, kb, 1.0);      // kl = Tp ^ kb * Tp;

    // transform resisting forces from the basic system to local coordinates
    static Vector pl(7);
    pl.addMatrixTransposeVector(0.0, Tp, pb, 1.0);    // pl = Tp ^ pb;

    // transform tangent  stiffness matrix from local to global coordinates
    // static Matrix kg(12,12);

    // compute the tangent stiffness matrix in global coordinates
    kg.addMatrixTripleProduct(0.0, T, kl, 1.0);

    static Vector m(6);
    for (int i = 0; i < 6; i++)
      m(i) = 0.5*pl(i)/cos(ul(i));


    //   ks = t'*kl*t + ks1 + t * diag (m .* tan(thetal))*t' + ...
    //        m(4)*(ks2r2t3_u3 + ks2r3u2_t2) + ...
    //        m(2)*ks2r2t1 + m(3)*ks2r3t1 + ...
    //        m(5)*ks2r2u1 + m(6)*ks2r3u1 + ...
    //        ks3 + ks3' + ks4 + ks5;

    // ksigma1 -------------------------------
    //   ks1_11 =  a*pl(6);
    //   ks1 = [ ks1_11  o  -ks1_11  o;
    //             o     o      o    o;
    //          -ks1_11  o   ks1_11  o;
    //             o     o      o    o];


    kg.Assemble(A, 0, 0,  pl(6));
    kg.Assemble(A, 0, 6, -pl(6));
    kg.Assemble(A, 6, 0, -pl(6));
    kg.Assemble(A, 6, 6,  pl(6));


    // ksigma3 -------------------------------
    //  kbar2 = -Lr2*(m(3)*S(rI3) + m(1)*S(rI1)) + ...
    //           Lr3*(m(3)*S(rI2) - m(2)*S(rI1)) ;

    //  kbar4 =  Lr2*(m(3)*S(rJ3) - m(4)*S(rJ1)) - ...
    //           Lr3*(m(3)*S(rJ2) + m(5)*S(rJ1));

    //     ks3 = [o kbar2 o kbar4];

    static Matrix3D Sm;
    static Matrix kbar(12,3);

    Sm.zero();
    Sm.addSpin(rI3,  m(3));
    Sm.addSpin(rI1,  m(1));
    kbar.addMatrixProduct(0.0, Lr2, Sm, -1.0);

    Sm.zero();
    Sm.addSpin(rI2,  m[3]);
    Sm.addSpin(rI1, -m[2]);
    kbar.addMatrixProduct(1.0, Lr3, Sm,  1.0);

    kg.Assemble(kbar, 0, 3, 1.0);
    kg.AssembleTranspose(kbar, 3, 0, 1.0);

    Sm.zero();
    Sm.addSpin(rJ3,  m(3));
    Sm.addSpin(rJ1, -m(4));
    kbar.addMatrixProduct(0.0, Lr2, Sm, 1.0);

    Sm.zero();
    Sm.addSpin(rJ2, m(3));
    Sm.addSpin(rJ1, m(5));
    kbar.addMatrixProduct(1.0, Lr3, Sm,  -1.0);

    kg.Assemble(kbar, 0, 9, 1.0);
    kg.AssembleTranspose(kbar, 9, 0, 1.0);


    // Ksigma4 -------------------------------
    // Ks4_22 =  m(3)*( S(e2)*S(rI3) - S(e3)*S(rI2)) + ...
    //           m(1)*(-S(e1)*S(rI2) + S(e2)*S(rI1)) + ...
    //           m(2)*(-S(e1)*S(rI3) + S(e3)*S(rI1));

    // Ks4_44 = -m(3)*( S(e2)*S(rJ3) - S(e3)*S(rJ2)) + ...
    //           m(4)*(-S(e1)*S(rJ2) + S(e2)*S(rJ1)) + ...
    //           m(5)*(-S(e1)*S(rJ3) + S(e3)*S(rJ1));

    // Ks4 = [   O    O     O    O;
    //           O  Ks4_22  O    O;
    //           O    O     O    O;
    //           O    O     O  Ks4_44];

    static Matrix3D ks33;

    ks33.zero();
    ks33.addSpinProduct(e2, rI3,  m[3]);
    ks33.addSpinProduct(e3, rI2, -m[3]);
    ks33.addSpinProduct(e2, rI1,  m[1]);
    ks33.addSpinProduct(e1, rI2, -m[1]);
    ks33.addSpinProduct(e3, rI1,  m[2]);
    ks33.addSpinProduct(e1, rI3, -m[2]);

    kg.Assemble(ks33, 3, 3, 1.0);

    ks33.zero();
    ks33.addSpinProduct(e2, rJ3, -m[3]);
    ks33.addSpinProduct(e3, rJ2,  m[3]);
    ks33.addSpinProduct(e2, rJ1,  m[4]);
    ks33.addSpinProduct(e1, rJ2, -m[4]);
    ks33.addSpinProduct(e3, rJ1,  m[5]);
    ks33.addSpinProduct(e1, rJ3, -m[5]);

    kg.Assemble(ks33, 9, 9, 1.0);


    // Ksigma5 -------------------------------
    //
    //  Ks5 = [ Ks5_11   Ks5_12 -Ks5_11   Ks5_14;
    //          Ks5_12t    O    -Ks5_12t   O;
    //         -Ks5_11  -Ks5_12  Ks5_11  -Ks5_14;
    //          Ks5_14t     O   -Ks5_14t   O];

    // v = (1/Ln)*(m(2)*rI2 + m(3)*rI3 + m(5)*rJ2 + m(6)*rJ3);
    OPS_STATIC Vector3D v;
    v.addVector (0.0, rI2, m[1]);
    v.addVector (1.0, rI3, m[2]);
    v.addVector (1.0, rJ2, m[4]);
    v.addVector (1.0, rJ3, m[5]);
    v /= Ln;

    //Ks5_11 = A*v*e1' + e1*v'*A + (e1'*v)*A;
    static Matrix3D m33;
    double  e1tv = e1.dot(v);   // dot product e1. v

    ks33.zero();
    ks33.addMatrix(A, e1tv);

    m33.zero();
    m33.addTensorProduct(v, e1, 1.0);

    ks33.addMatrixProduct(A, m33, 1.0);

    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        m33(i,j) = e1[i]*v[j];

    ks33.addMatrixProduct(m33, A, 1.0);

    kg.Assemble(ks33, 0, 0,  1.0);
    kg.Assemble(ks33, 0, 6, -1.0);
    kg.Assemble(ks33, 6, 0, -1.0);
    kg.Assemble(ks33, 6, 6,  1.0);

    //Ks5_12 = -(m(2)*A*S(rI2) + m(3)*A*S(rI3));

    ks33.zero();
    ks33.addMatrixSpinProduct(A, rI2, -m[1]);
    ks33.addMatrixSpinProduct(A, rI3, -m[2]);

    kg.Assemble(ks33, 0, 3,  1.0);
    kg.Assemble(ks33, 6, 3, -1.0);
    kg.AssembleTranspose(ks33, 3, 0,  1.0);
    kg.AssembleTranspose(ks33, 3, 6, -1.0);

    //  Ks5_14 = -(m(5)*A*S(rJ2) + m(6)*A*S(rJ3));

    ks33.zero();
    ks33.addMatrixSpinProduct(A, rJ2, -m[4]);
    ks33.addMatrixSpinProduct(A, rJ3, -m[5]);

    kg.Assemble(ks33, 0, 9,  1.0);
    kg.Assemble(ks33, 6, 9, -1.0);

    kg.AssembleTranspose(ks33, 9, 0,  1.0);
    kg.AssembleTranspose(ks33, 9, 6, -1.0);

    // Ksigma -------------------------------
    OPS_STATIC Vector3D rm;

    rm = rI3;
    rm.addVector(1.0, rJ3, -1.0);
    kg.addMatrix(getKs2Matrix(A, e1, r1, Ln, r2, rm), m(3));

//  rm = rJ2;
    rm.addVector(0.0, rJ2, -1.0);
    rm.addVector(1.0, rI2, -1.0);
    kg.addMatrix(getKs2Matrix(A, e1, r1, Ln, r3,  rm), m(3));
    kg.addMatrix(getKs2Matrix(A, e1, r1, Ln, r2, rI1), m(1));
    kg.addMatrix(getKs2Matrix(A, e1, r1, Ln, r3, rI1), m(2));
    kg.addMatrix(getKs2Matrix(A, e1, r1, Ln, r2, rJ1), m(4));
    kg.addMatrix(getKs2Matrix(A, e1, r1, Ln, r3, rJ1), m(5));

    //  T * diag (M .* tan(thetal))*T'

    for (int k = 0; k < 6; k++) {
      const double factor = pl(k) * tan(ul(k));
      for (int i = 0; i < 12; i++) {
        const double Tki = T(k,i);
        for (int j = 0; j < 12; j++)
          kg(i,j) += Tki * factor * T(k,j);
      }
    }

    return kg;
}


const Matrix &
CorotCrdTransf3d02::getInitialGlobalStiffMatrix(const Matrix &kb)
{
    // transform tangent stiffness matrix from the basic system to local coordinates
    static Matrix kl(7,7);
    kl.addMatrixTripleProduct(0.0, Tp, kb, 1.0);      // kl = Tp ^ kb * Tp;

    // transform tangent  stiffness matrix from local to global coordinates
    //static Matrix kg(12,12);

    // compute the tangent stiffness matrix in global coordinates
    kg.addMatrixTripleProduct(0.0, T, kl, 1.0);

    return kg;
}


int
CorotCrdTransf3d02::getLocalAxes(Vector &XAxis, Vector &YAxis, Vector &ZAxis)
{
    // element projection

    static Vector dx(3);

    dx = (nodeJPtr->getCrds() + nodeJOffset) - (nodeIPtr->getCrds() + nodeIOffset);
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
        opserr << "\nCorotCrdTransf3d02::computeElemtLengthAndOrien: 0 length\n";
        return -2;
    }

    // calculate the element local x axis components (direction cossines)
    // wrt to the global coordinates
    xAxis = dx/L;

    XAxis(0) = xAxis(0);
    XAxis(1) = xAxis(1);
    XAxis(2) = xAxis(2);

    // calculate the cross-product y = v * x
    static Vector yAxis(3), zAxis(3);

    yAxis(0) = vAxis(1)*xAxis(2) - vAxis(2)*xAxis(1);
    yAxis(1) = vAxis(2)*xAxis(0) - vAxis(0)*xAxis(2);
    yAxis(2) = vAxis(0)*xAxis(1) - vAxis(1)*xAxis(0);

    double ynorm = yAxis.Norm();

    if (ynorm == 0) {
        opserr << "\nCorotCrdTransf3d02::getElementLengthAndOrientation";
        opserr << "\nvector v that defines plane xz is parallel to x axis\n";
        return -3;
    }

    yAxis /= ynorm;
    YAxis(0) = yAxis(0);
    YAxis(1) = yAxis(1);
    YAxis(2) = yAxis(2);

    // calculate the cross-product z = x * y

    zAxis(0) = xAxis(1)*yAxis(2) - xAxis(2)*yAxis(1);
    zAxis(1) = xAxis(2)*yAxis(0) - xAxis(0)*yAxis(2);
    zAxis(2) = xAxis(0)*yAxis(1) - xAxis(1)*yAxis(0);

    ZAxis(0) = zAxis(0);
    ZAxis(1) = zAxis(1);
    ZAxis(2) = zAxis(2);

    for (int i = 0; i < 3; i++) {
        R0(i,0) = xAxis(i);
        R0(i,1) = yAxis(i);
        R0(i,2) = zAxis(i);
    }

    return 0;
}


double
CorotCrdTransf3d02::getInitialLength()
{
  return L;
}


double
CorotCrdTransf3d02::getDeformedLength()
{
  return Ln;
}


CrdTransf *
CorotCrdTransf3d02::getCopy3d()
{
  // create a new instance of CorotCrdTransf3d02

  CorotCrdTransf3d02 *theCopy =
    new CorotCrdTransf3d02 (this->getTag(), vAxis, nodeIOffset, nodeJOffset);

  theCopy->nodeIPtr = nodeIPtr;
  theCopy->nodeJPtr = nodeJPtr;
  theCopy->xAxis = xAxis;
  theCopy->L  = L;
  theCopy->Ln = Ln;
  theCopy->R0 = R0;
  theCopy->alphaIq = alphaIq;
  theCopy->alphaJq = alphaJq;
  theCopy->alphaIqcommit = alphaIqcommit;
  theCopy->alphaJqcommit = alphaJqcommit;
  theCopy->ul = ul;
  theCopy->ulcommit = ulcommit;

  return theCopy;
}


int
CorotCrdTransf3d02::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(48);
  for (int i = 0; i<7; i++)
    data(i) = ulcommit(i);

  for (int j = 0; j<4; j++) {
    data(7+j) = alphaIqcommit[j];
    data(11+j) = alphaJqcommit[j];
  }

  for (int k=0; k<3; k++) {
    data(15+k) = xAxis(k);
    data(18+k) = vAxis(k);
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
    opserr << " CorotCrdTransf3d02::sendSelf() - data could not be sent\n" ;
    return -1;
  }

  return 0;
}


int
CorotCrdTransf3d02::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static Vector data(48);
  if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
    opserr << " CorotCrdTransf3d02::recvSelf() - data could not be received\n" ;
    return -1;
  }

  for (int i = 0; i<7; i++)
    ulcommit(i) = data(i);

  for (int j = 0; j<4; j++) {
    alphaIqcommit[j] = data(7+j);
    alphaJqcommit[j] = data(11+j);
  }

  for (int k=0; k<3; k++) {
    xAxis(k) = data(15+k);
    vAxis(k) = data(18+k);
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
  alphaIq = alphaIqcommit;
  alphaJq = alphaJqcommit;

  initialDispChecked = true;
  return 0;
}


const Matrix &
CorotCrdTransf3d02::getGlobalMatrixFromLocal(const Matrix &local)
{
    static Matrix Tlg(12,12);
    this->compTransfMatrixLocalGlobal(Tlg);  // OPTIMIZE LATER
    kg.addMatrixTripleProduct(0.0, Tlg, local, 1.0);  // OPTIMIZE LATER

    return kg;
}


const Vector &
CorotCrdTransf3d02::getPointGlobalCoordFromLocal(const Vector &xl)
{
    static Vector xg(3);
    opserr << " CorotCrdTransf3d02::getPointGlobalCoordFromLocal: not implemented yet" ;

    return xg;
}


const Vector &
CorotCrdTransf3d02::getPointGlobalDisplFromBasic(double xi, const Vector &uxb)
{
    static Vector uxg(3);
    opserr << " CorotCrdTransf3d02::getPointGlobalDisplFromBasic: not implemented yet" ;

    return uxg;
}


const Vector &
CorotCrdTransf3d02::getPointLocalDisplFromBasic(double xi, const Vector &uxb)
{
    static Vector uxg(3);
    opserr << " CorotCrdTransf3d02::getPointLocalDisplFromBasic: not implemented yet" ;

    return uxg;
}


void
CorotCrdTransf3d02::Print(OPS_Stream &s, int flag)
{

  if (flag == OPS_PRINT_CURRENTSTATE) {
          s << "\nCrdTransf: " << this->getTag() << " Type: CorotCrdTransf3d02";
          s << "\tvAxis: " << vAxis;
          s << "\tnodeI Offset: " << nodeIOffset;
          s << "\tnodeJ Offset: " << nodeJOffset;
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
      s << "\t\t\t{\"name\": \"" << this->getTag() << "\", \"type\": \"CorotCrdTransf3d02\"";
      s << ", \"vecInLocXZPlane\": [" << vAxis(0) << ", " << vAxis(1) << ", " << vAxis(2) << "]";
      if (nodeIOffset != 0)
          s << ", \"iOffset\": [" << nodeIOffset[0] << ", " << nodeIOffset[1] << ", " << nodeIOffset[2] << "]";
      if (nodeJOffset != 0)
          s << ", \"jOffset\": [" << nodeJOffset[0] << ", " << nodeJOffset[1] << ", " << nodeJOffset[2] << "]";
      s << "}";
  }
}

#if 0
// CMP: commented out, but keeping around because the comments
// are good
void
CorotCrdTransf3d02::compTransfMatrixBasicGlobalNew()
{
    // extract columns of rotation matrices

    OPS_STATIC Vector3D  r1,  r2,  r3,
                         e1,  e2,  e3,
                        rI1, rI2, rI3,
                        rJ1, rJ2, rJ3;

    for (int k = 0; k < 3; k ++) {
      r1[k] = Rbar(k,0);
      r2[k] = Rbar(k,1);
      r3[k] = Rbar(k,2);

      e1[k] = e(k,0);
      e2[k] = e(k,1);
      e3[k] = e(k,2);

      rI1[k] = RI(k,0);
      rI2[k] = RI(k,1);
      rI3[k] = RI(k,2);

      rJ1[k] = RJ(k,0);
      rJ2[k] = RJ(k,1);
      rJ3[k] = RJ(k,2);
    }

    // compute the transformation matrix from the basic to the
    // global system

    //   A = (1/Ln)*(I - e1*e1');
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        A(i,j) = (double(i==j) - e1[i]*e1[j])/Ln;


    Lr2 = this->getLMatrix (r2);
    Lr3 = this->getLMatrix (r3);

    static Matrix Sr1(3,3), Sr2(3,3), Sr3(3,3);
    static Vector Se(3), At(3);

    // O = zeros(3,1);
    // hI1 = [      O', (-S(rI3)*e2 + S(rI2)*e3)',        O', O']';
    // hI2 = [(A*rI3)', (-S(rI3)*e1 + S(rI1)*e3)', -(A*rI3)', O']';
    // hI3 = [(A*rI2)', (-S(rI2)*e1 + S(rI1)*e2)', -(A*rI2)', O']';
    // hJ1 = [      O', O',        O', (-S(rJ3)*e2 + S(rJ2)*e3)']';
    // hJ2 = [(A*rJ3)', O', -(A*rJ3)', (-S(rJ3)*e1 + S(rJ1)*e3)']';
    // hJ3 = [(A*rJ2)', O', -(A*rJ2)', (-S(rJ2)*e1 + S(rJ1)*e2)']';

    static Vector hI1(12), hI2(12), hI3(12),
                  hJ1(12), hJ2(12), hJ3(12);

    SO3::Spin(rI1, Sr1 );
    SO3::Spin(rI2, Sr2 );
    SO3::Spin(rI3, Sr3 );

    // hI1 = [      O', (-S(rI3)*e2 + S(rI2)*e3)',        O', O']';
    Se.addMatrixVector(0.0, Sr3, e2, -1.0);     // (-S(rI3)*e2 + S(rI2)*e3)
    Se.addMatrixVector(1.0, Sr2, e3,  1.0);
    for (int i = 0; i<3; i++)
      hI1(i+3) = Se(i);

    // hI2 = [(A*rI3)', (-S(rI3)*e1 + S(rI1)*e3)', -(A*rI3)', O']';
    At.addMatrixVector(0.0,  A, rI3,  1.0);
    Se.addMatrixVector(0.0, Sr3, e1, -1.0);     // (-S(rI3)*e1 + S(rI1)*e3)
    Se.addMatrixVector(1.0, Sr1, e3,  1.0);

    for (int i = 0; i < 3; i++) {
      hI2(i  ) =  At(i);
      hI2(i+3) =  Se(i);
      hI2(i+6) = -At(i);
    }

    // hI3 = [(A*rI2)', (-S(rI2)*e1 + S(rI1)*e2)', -(A*rI2)', O']';
    At.addMatrixVector(0.0, A, rI2, 1.0);
    Se.addMatrixVector(0.0, Sr2, e1, -1.0);     // (-S(rI2)*e1 + S(rI1)*e2)'
    Se.addMatrixVector(1.0, Sr1, e2,  1.0);

    for (int i = 0; i < 3; i++) {
      hI3(i  ) =  At(i);
      hI3(i+3) =  Se(i);
      hI3(i+6) = -At(i);
    }

    SO3::Spin(rJ1, Sr1);
    SO3::Spin(rJ2, Sr2);
    SO3::Spin(rJ3, Sr3);

    // hJ1 = [      O', O',        O', (-S(rJ3)*e2 + S(rJ2)*e3)']';
    Se.addMatrixVector(0.0, Sr3, e2, -1.0);    // -S(rJ3)*e2 + S(rJ2)*e3
    Se.addMatrixVector(1.0, Sr2, e3,  1.0);
    for (int i = 0; i < 3; i++)
        hJ1(i+9) =  Se(i);


    // hJ2 = [(A*rJ3)', O', -(A*rJ3)', (-S(rJ3)*e1 + S(rJ1)*e3)']';
    At.addMatrixVector(0.0, A, rJ3, 1.0);
    Se.addMatrixVector(0.0, Sr3, e1, -1.0);     // (-S(rJ3)*e1 + S(rJ1)*e3)
    Se.addMatrixVector(1.0, Sr1, e3,  1.0);

    for (int i = 0; i < 3; i++) {
      hJ2(i  ) =  At(i);
      hJ2(i+6) = -At(i);
      hJ2(i+9) =  Se(i);
    }

    // hJ3 = [(A*rJ2)', O', -(A*rJ2)', (-S(rJ2)*e1 + S(rJ1)*e2)']';
    At.addMatrixVector(0.0, A, rJ2, 1.0);
    Se.addMatrixVector(0.0, Sr2, e1, -1.0);     // (-S(rJ2)*e1 + S(rJ1)*e2)
    Se.addMatrixVector(1.0, Sr1, e2,  1.0);

    for (int i = 0; i < 3; i++) {
      hJ3(i  ) =  At(i);
      hJ3(i+6) = -At(i);
      hJ3(i+9) =  Se(i);
    }

    // f1 =  [-e1' O' e1' O'];
    // f2  = ( Lr2*rI1 + hI3)'./(2*(cos(thetalI(3))));
    // f3  = ( Lr2*rJ1 + hJ3)'./(2*(cos(thetalJ(3))));
    // f4  = (-Lr3*rI1 - hI2)'./(2*(cos(thetalI(2))));
    // f5  = (-Lr3*rJ1 - hJ2)'./(2*(cos(thetalJ(2))));
    // f6I = ( Lr3*rI2 - Lr2*rI3 + hI1)'./(2*(cos(thetalI(1))));
    // f6J = ( Lr3*rJ2 - Lr2*rJ3 + hJ1)'./(2*(cos(thetalJ(1))));

    // F = [f1;
    //     f2;
    //     f3;
    //     f4;
    //     f5;
    //     f6J-f6I];

    // T = F'
    T.Zero();
    static Vector Lr(12);

    // f1 =  [-e1' O' e1' O'];
    for (int i = 0; i<3; i++) {
        T(i  ,0) = -e1[i];
        T(i+3,0) = e1[i];
    }

    static Vector thetaI(3);
    static Vector thetaJ(3);


    thetaI(0) = ul(0);
    thetaI(1) = -ul(2);
    thetaI(2) = ul(1);

    thetaJ(0) = ul(3);
    thetaJ(1) = -ul(5);
    thetaJ(2) = ul(4);


    opserr << "thetaI: " << thetaI;
    opserr << "thetaJ: " << thetaJ;

    double c;

    // f2  = ( Lr2*rI1 + hI3)'./(2*(cos(thetalI(3))));
    Lr.addMatrixVector(0.0, Lr2, rI1,  1.0);
    Lr += hI3;
    c = 0.5*cos(thetaI(2));
    for (int i = 0; i<12; i++)
      T(1,i) = Lr(i)*c;

    // f3  = ( Lr2*rJ1 + hJ3)'./(2*(cos(thetalJ(3))));
    Lr.addMatrixVector(0.0, Lr2, rJ1,  1.0);
    Lr += hJ3;
    c = 0.5*cos(thetaJ(2));
    for (int i = 0; i<12; i++)
      T(2,i) = Lr(i)*c;

    // f4  = (-Lr3*rI1 - hI2)'./(2*(cos(thetalI(2))));
    Lr.addMatrixVector(0.0, Lr3, rI1,  -1.0);
    Lr -= hI2;
    c = 0.5*cos(thetaI(1));
    for (int i = 0; i<12; i++)
      T(3,i) = Lr(i)*c;

    // f5  = (-Lr3*rJ1 - hJ2)'./(2*(cos(thetalJ(2))));
    Lr.addMatrixVector(0.0, Lr3, rJ1,  -1.0);
    Lr -= hJ2;
    c = 0.5*cos(thetaJ(1));
    for (int i = 0; i<12; i++)
      T(4,i) = Lr(i)*c;

    // f6I = ( Lr3*rI2 - Lr2*rI3 + hI1)'./(2*(cos(thetalI(1))));
    Lr.addMatrixVector(0.0, Lr3, rI2,  1.0);
    Lr.addMatrixVector(1.0, Lr2, rI3, -1.0);
    Lr += hI1;
    c = 0.5*cos(thetaI(0));
    for (int i = 0; i<12; i++)
      T(5,i) = Lr(i)*c;

    // f6J = ( Lr3*rJ2 - Lr2*rJ3 + hJ1)'./(2*(cos(thetalJ(1))));
    Lr.addMatrixVector(0.0, Lr3, rJ2,  1.0);
    Lr.addMatrixVector(1.0, Lr2, rJ3, -1.0);
    Lr += hJ1;
    c = 0.5*cos(thetaI(0));
    for (int i = 0; i<12; i++)
      T(6,i) -= Lr(i)*c;
}
#endif

