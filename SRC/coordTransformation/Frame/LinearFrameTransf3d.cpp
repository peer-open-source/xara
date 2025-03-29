//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the implementation for the
// LinearFrameTransf3d class. LinearFrameTransf3d is a linear
// transformation for a planar frame between the global
// and basic coordinate systems
//
// Written: Remo Magalhaes de Souza
// Created: 04/2000
//
// Modified: 04/2005 Andreas Schellenberg (getBasicTrialVel, getBasicTrialAccel)
//
// Sensitivity: QGu UCSD 2009
//
#include <Vector.h>
#include <Matrix.h>
#include <Matrix3D.h>
#include <Node.h>
#include <Channel.h>
#include <Logging.h>
#include <string>
#include <LinearFrameTransf3d.h>
#include "blk3x12x3.h"

using namespace OpenSees;

static MatrixND<3,3>
FrameOrientationGradient(const Vector3D& xi, const Vector3D& xj, const Vector3D& vz, int di, int dj, int dv)
{
    Vector3D v1  = xj - xi;
    double L     = v1.norm();
    Vector3D e1  = v1/L;

    Vector3D v2  = vz.cross(e1);

    Vector3D e2 = v2 / v2.norm();
//  Vector3D v3 = e1.cross(e2);
//  Vector3D e3 = v3 / v3.norm();

    //
    Vector3D dvz{0.0};
    Vector3D dxi{0.0};
    Vector3D dxj{0.0};

    if (di != 0)
      dxi(di-1) = 1.0;
    if (dj != 0)
      dxj(dj-1) = 1.0;
    if (dv != 0)
      dvz(dv-1) = 1.0;

    //

    double   dL  = 1/L*(xj - xi).dot(dxj - dxi);
    Vector3D dv1 = dxj - dxi;
    Vector3D de1 = 1/(L*L)*(dv1*L - v1*dL);

    double L2    = v2.norm();
    Vector3D dv2 = dvz.cross(e1) + vz.cross(de1);
    double dL2   = 1/L2*v2.dot(dv2);
    Vector3D de2 = 1/(L2*L2)*(dv2*L2 - v2*dL2);

    Vector3D de3 = de1.cross(e2) + e1.cross(de2);

    MatrixND<3,3> dR;
    dR(0,0) = de1(0);
    dR(1,0) = de1(1);
    dR(2,0) = de1(2);

    dR(0,1) = de2(0);
    dR(1,1) = de2(1);
    dR(2,1) = de2(2);

    dR(0,2) = de3(0);
    dR(1,2) = de3(1);
    dR(2,2) = de3(2);

    return dR;

//  return np.stack([de1,de2,de3])
}


// initialize static variables
Matrix LinearFrameTransf3d::kg(12, 12);

// constructor:
LinearFrameTransf3d::LinearFrameTransf3d(int tag, const Vector &vecInLocXZPlane)
    : FrameTransform3d(tag, CRDTR_TAG_LinearFrameTransf3d),
      nodeIPtr(nullptr), nodeJPtr(nullptr),
      nodeIOffset(0), nodeJOffset(0), 
      L(0),
      nodeIInitialDisp(0),
      nodeJInitialDisp(0), initialDispChecked(false)
{
  R.zero();

  for (int i=0; i<3; i++)
    vz[i] = vecInLocXZPlane[i];

  R(0,2) = vz(0);
  R(1,2) = vz(1);
  R(2,2) = vz(2);
}

// constructor:
LinearFrameTransf3d::LinearFrameTransf3d(int tag, const Vector &vecInLocXZPlane,
                                         const Vector &rigJntOffset1,
                                         const Vector &rigJntOffset2)

  : FrameTransform3d(tag, CRDTR_TAG_LinearFrameTransf3d), 
    nodeIPtr(0), nodeJPtr(0),
    nodeIOffset(nullptr), nodeJOffset(nullptr), 
    L(0),
    nodeIInitialDisp(0),
    nodeJInitialDisp(0), 
    initialDispChecked(false)
{
  R.zero();

  for (int i=0; i<3; i++)
    vz[i] = vecInLocXZPlane[i];

  R(0,2) = vecInLocXZPlane(0);
  R(1,2) = vecInLocXZPlane(1);
  R(2,2) = vecInLocXZPlane(2);

  // check rigid joint offset for node I
  if (rigJntOffset1.Size() != 3) {
    opserr << "LinearFrameTransf3d::LinearFrameTransf3d:  Invalid rigid joint "
              "offset vector for node I\n";
    opserr << "Size must be 3\n";

  } else if (rigJntOffset1.Norm() > 0.0) {
    nodeIOffset    = new double[3];
    nodeIOffset[0] = rigJntOffset1(0);
    nodeIOffset[1] = rigJntOffset1(1);
    nodeIOffset[2] = rigJntOffset1(2);
  }

  // check rigid joint offset for node J
  if (rigJntOffset2.Size() != 3) {
    opserr << "LinearFrameTransf3d::LinearFrameTransf3d:  Invalid rigid joint "
              "offset vector for node J\n";
    opserr << "Size must be 3\n";

  } else if (rigJntOffset2.Norm() > 0.0) {
    nodeJOffset    = new double[3];
    nodeJOffset[0] = rigJntOffset2(0);
    nodeJOffset[1] = rigJntOffset2(1);
    nodeJOffset[2] = rigJntOffset2(2);
  }
}

// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
LinearFrameTransf3d::LinearFrameTransf3d()
    : FrameTransform3d(0, CRDTR_TAG_LinearFrameTransf3d), nodeIPtr(0), nodeJPtr(0),
      nodeIOffset(0), nodeJOffset(0), L(0), nodeIInitialDisp(0),
      nodeJInitialDisp(0), initialDispChecked(false)
{
  R.zero();
}

// destructor:
LinearFrameTransf3d::~LinearFrameTransf3d()
{
  if (nodeIOffset)
    delete[] nodeIOffset;
  if (nodeJOffset)
    delete[] nodeJOffset;
  if (nodeIInitialDisp != 0)
    delete[] nodeIInitialDisp;
  if (nodeJInitialDisp != 0)
    delete[] nodeJInitialDisp;
}

int
LinearFrameTransf3d::commitState()
{
  return 0;
}

int
LinearFrameTransf3d::revertToLastCommit()
{
  return 0;
}

int
LinearFrameTransf3d::revertToStart()
{
  return 0;
}

int
LinearFrameTransf3d::update()
{
  return 0;
}


int
LinearFrameTransf3d::initialize(Node *nodeIPointer, Node *nodeJPointer)
{
  int error;

  nodeIPtr = nodeIPointer;
  nodeJPtr = nodeJPointer;

  if ((!nodeIPtr) || (!nodeJPtr)) {
    opserr << "\nLinearFrameTransf3d::initialize";
    opserr << "\ninvalid pointers to the element nodes\n";
    return -1;
  }

  // see if there is some initial displacements at nodes
  if (initialDispChecked == false) {
    const Vector &nodeIDisp = nodeIPtr->getDisp();
    const Vector &nodeJDisp = nodeJPtr->getDisp();
    for (int i = 0; i < 6; i++)
      if (nodeIDisp(i) != 0.0) {
        nodeIInitialDisp = new double[6];
        for (int j = 0; j < 6; j++)
          nodeIInitialDisp[j] = nodeIDisp(j);
        i = 6;
      }

    for (int j = 0; j < 6; j++)
      if (nodeJDisp(j) != 0.0) {
        nodeJInitialDisp = new double[6];
        for (int i = 0; i < 6; i++)
          nodeJInitialDisp[i] = nodeJDisp(i);
        j = 6;
      }

    initialDispChecked = true;
  }

  // get element length and orientation
  if ((error = this->computeElemtLengthAndOrient()))
    return error;

  static Vector XAxis(3);
  static Vector YAxis(3);
  static Vector ZAxis(3);

  // fill 3by3 rotation matrix, R
  if ((error = this->getLocalAxes(XAxis, YAxis, ZAxis)))
    return error;

  return 0;
}


int
LinearFrameTransf3d::computeElemtLengthAndOrient()
{

  const Vector &XI = nodeIPtr->getCrds();
  const Vector &XJ = nodeJPtr->getCrds();

  for (int i=0; i<3; i++) {
    xi[i] = XI[i];
    xj[i] = XJ[i];
  }
  
  Vector3D dx = xj - xi;

  if (nodeJOffset != 0) {
    dx(0) += nodeJOffset[0];
    dx(1) += nodeJOffset[1];
    dx(2) += nodeJOffset[2];
  }

  if (nodeIOffset != 0) {
    dx(0) -= nodeIOffset[0];
    dx(1) -= nodeIOffset[1];
    dx(2) -= nodeIOffset[2];
  }

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
  L = dx.norm();

  if (L == 0.0)
    return -2;

  // Calculate the element local x axis components (direction cosines)
  // wrt to the global coordinates
  R(0,0) = dx(0) / L;
  R(1,0) = dx(1) / L;
  R(2,0) = dx(2) / L;

  return 0;
}

void
LinearFrameTransf3d::compTransfMatrixLocalGlobal(Matrix &Tlg)
{
  // setup transformation matrix from local to global
  Tlg.Zero();

  Tlg(0, 0) = Tlg(3, 3) = Tlg(6, 6) = Tlg( 9,  9) = R(0,0);
  Tlg(0, 1) = Tlg(3, 4) = Tlg(6, 7) = Tlg( 9, 10) = R(1,0);
  Tlg(0, 2) = Tlg(3, 5) = Tlg(6, 8) = Tlg( 9, 11) = R(2,0);

  Tlg(1, 0) = Tlg(4, 3) = Tlg(7, 6) = Tlg(10,  9) = R(0,1);
  Tlg(1, 1) = Tlg(4, 4) = Tlg(7, 7) = Tlg(10, 10) = R(1,1);
  Tlg(1, 2) = Tlg(4, 5) = Tlg(7, 8) = Tlg(10, 11) = R(2,1);

  Tlg(2, 0) = Tlg(5, 3) = Tlg(8, 6) = Tlg(11,  9) = R(0,2);
  Tlg(2, 1) = Tlg(5, 4) = Tlg(8, 7) = Tlg(11, 10) = R(1,2);
  Tlg(2, 2) = Tlg(5, 5) = Tlg(8, 8) = Tlg(11, 11) = R(2,2);
}

int
LinearFrameTransf3d::getLocalAxes(Vector &XAxis, Vector &YAxis, Vector &ZAxis)
{
  // Compute y = v cross x
  // Note: v(i) is stored in R(i,2)
  static Vector vAxis(3);
  vAxis(0) = R(0,2);
  vAxis(1) = R(1,2);
  vAxis(2) = R(2,2);
  if (nodeIPtr == nullptr) {
    ZAxis(0) = vAxis(0);
    ZAxis(1) = vAxis(1);
    ZAxis(2) = vAxis(2);
    return 0;
  }

  static Vector3D e1;
  e1[0] = R(0,0);
  e1[1] = R(1,0);
  e1[2] = R(2,0);
  XAxis(0) = e1[0];
  XAxis(1) = e1[1];
  XAxis(2) = e1[2];

  static Vector3D e2;
  e2(0) = vAxis(1) * e1(2) - vAxis(2) * e1(1);
  e2(1) = vAxis(2) * e1(0) - vAxis(0) * e1(2);
  e2(2) = vAxis(0) * e1(1) - vAxis(1) * e1(0);

  double ynorm = e2.norm();

  if (ynorm == 0) {
    opserr << "\nLinearFrameTransf3d::getLocalAxes";
    opserr << "\nvector v that defines plane xz is parallel to x axis\n";
    return -3;
  }

  e2 /= ynorm;

  YAxis(0) = e2[0];
  YAxis(1) = e2[1];
  YAxis(2) = e2[2];

  // Compute z = x cross y
  Vector3D e3 = e1.cross(e2);

  ZAxis(0) = e3[0];
  ZAxis(1) = e3[1];
  ZAxis(2) = e3[2];

  // Fill in transformation matrix
  R(0,1) = e2[0];
  R(1,1) = e2[1];
  R(2,1) = e2[2];

  R(0,2) = e3[0];
  R(1,2) = e3[1];
  R(2,2) = e3[2];

  formOffsets(R, nodeIOffset, nodeJOffset, this->RWI, this->RWJ);

  return 0;
}

double
LinearFrameTransf3d::getInitialLength()
{
  return L;
}

double
LinearFrameTransf3d::getDeformedLength()
{
  return L;
}


const Vector &
LinearFrameTransf3d::getBasicTrialDisp()
{
  // determine global displacements
  const Vector &disp1 = nodeIPtr->getTrialDisp();
  const Vector &disp2 = nodeJPtr->getTrialDisp();

  static VectorND<12> ug;
  for (int i = 0; i < 6; i++) {
    ug[i]     = disp1(i);
    ug[i + 6] = disp2(i);
  }

  if (nodeIInitialDisp != 0) {
    for (int j = 0; j < 6; j++)
      ug[j] -= nodeIInitialDisp[j];
  }

  if (nodeJInitialDisp != 0) {
    for (int j = 0; j < 6; j++)
      ug[j + 6] -= nodeJInitialDisp[j];
  }

  double oneOverL = 1.0 / L;

  static VectorND<6> ub;
  static Vector wrapper(ub);
  VectorND<12> ul = getLocal<2,6>(ug, R, nodeIOffset, nodeJOffset);
  ub = getBasic(ul, oneOverL);
  return wrapper;
}

const Vector &
LinearFrameTransf3d::getBasicIncrDisp()
{
  // determine global displacements
  const Vector &disp1 = nodeIPtr->getIncrDisp();
  const Vector &disp2 = nodeJPtr->getIncrDisp();

  static VectorND<12> ug;
  for (int i = 0; i < 6; i++) {
    ug[i]     = disp1(i);
    ug[i + 6] = disp2(i);
  }

  double oneOverL = 1.0 / L;

  static VectorND<6> ub;
  static Vector wrapper(ub);
  VectorND<12> ul = getLocal<2,6>(ug, R, nodeIOffset, nodeJOffset);
  ub = getBasic(ul, oneOverL);

  return wrapper;
}

const Vector &
LinearFrameTransf3d::getBasicIncrDeltaDisp()
{
  // determine global displacements
  const Vector &disp1 = nodeIPtr->getIncrDeltaDisp();
  const Vector &disp2 = nodeJPtr->getIncrDeltaDisp();

  static VectorND<12> ug;
  for (int i = 0; i < 6; i++) {
    ug[i]     = disp1(i);
    ug[i + 6] = disp2(i);
  }

  double oneOverL = 1.0 / L;

  static VectorND<6> ub;
  static Vector wrapper(ub);

  VectorND<12> ul = getLocal<2,6>(ug, R, nodeIOffset, nodeJOffset);
  ub = getBasic(ul, oneOverL);
  return wrapper;
}

const Vector &
LinearFrameTransf3d::getBasicTrialVel()
{
  // determine global velocities
  const Vector &vel1 = nodeIPtr->getTrialVel();
  const Vector &vel2 = nodeJPtr->getTrialVel();

  static VectorND<12> vg;
  for (int i = 0; i < 6; i++) {
    vg[i]     = vel1(i);
    vg[i + 6] = vel2(i);
  }

  double oneOverL = 1.0 / L;

  static VectorND<6> ub;
  static Vector wrapper(ub);

  VectorND<12> vl = getLocal<2,6>(vg, R, nodeIOffset, nodeJOffset);
  ub = getBasic(vl, oneOverL);
  return wrapper;
}

const Vector &
LinearFrameTransf3d::getBasicTrialAccel()
{
  // determine global accelerations
  const Vector &accel1 = nodeIPtr->getTrialAccel();
  const Vector &accel2 = nodeJPtr->getTrialAccel();

  static VectorND<12> ag;
  for (int i = 0; i < 6; i++) {
    ag[i]     = accel1(i);
    ag[i + 6] = accel2(i);
  }

  static VectorND<6> ub;
  static Vector wrapper(ub);

  VectorND<12> al = getLocal<2,6>(ag, R, nodeIOffset, nodeJOffset);
  ub = getBasic(al, 1/L);
  return wrapper;

//static Vector ub //(6);
//  = getBasic(ag, oneOverL);
//return ub;
}

VectorND<12>
LinearFrameTransf3d::pushResponse(VectorND<12>&pl)
{
  return pushConstant(pl);
}

MatrixND<12,12>
LinearFrameTransf3d::pushResponse(MatrixND<12,12>&kl, const VectorND<12>& pl)
{
  return pushConstant(kl);
}

VectorND<12>
LinearFrameTransf3d::pushConstant(const VectorND<12>&pl) const
{
  // transform vector from local to global coordinates
  VectorND<12> pg;
  for (int i=0; i<4; i++)
    for (int j=0; j<3; j++)
      pg[i*3+j] = R(j,0)*pl[3*i] + R(j,1)*pl[3*i+1] + R(j,2)*pl[3*i+2];


  if (nodeIOffset) {
    pg[3] += -nodeIOffset[2] * pg[1] + nodeIOffset[1] * pg[2];
    pg[4] +=  nodeIOffset[2] * pg[0] - nodeIOffset[0] * pg[2];
    pg[5] += -nodeIOffset[1] * pg[0] + nodeIOffset[0] * pg[1];
  }

  if (nodeJOffset) {
    pg[ 9] += -nodeJOffset[2] * pg[7] + nodeJOffset[1] * pg[8];
    pg[10] +=  nodeJOffset[2] * pg[6] - nodeJOffset[0] * pg[8];
    pg[11] += -nodeJOffset[1] * pg[6] + nodeJOffset[0] * pg[7];
  }

  return pg;
}

MatrixND<12,12>
LinearFrameTransf3d::pushConstant(const MatrixND<12,12>& kl)
{

  MatrixND<12,12> kg;

  // Transform local stiffness to global system
  // First compute kl*T_{lg}
  static double tmp[12][12];  // Temporary storage
  for (int m = 0; m < 12; m++) {
    tmp[m][0] = kl(m, 0) * R(0,0) + kl(m, 1) * R(0,1) + kl(m, 2) * R(0,2);
    tmp[m][1] = kl(m, 0) * R(1,0) + kl(m, 1) * R(1,1) + kl(m, 2) * R(1,2);
    tmp[m][2] = kl(m, 0) * R(2,0) + kl(m, 1) * R(2,1) + kl(m, 2) * R(2,2);

    tmp[m][3] = kl(m, 3) * R(0,0) + kl(m, 4) * R(0,1) + kl(m, 5) * R(0,2);
    tmp[m][4] = kl(m, 3) * R(1,0) + kl(m, 4) * R(1,1) + kl(m, 5) * R(1,2);
    tmp[m][5] = kl(m, 3) * R(2,0) + kl(m, 4) * R(2,1) + kl(m, 5) * R(2,2);

    if (nodeIOffset) {
      tmp[m][3] += kl(m, 0) * RWI[0][0] + kl(m, 1) * RWI[1][0] + kl(m, 2) * RWI[2][0];
      tmp[m][4] += kl(m, 0) * RWI[0][1] + kl(m, 1) * RWI[1][1] + kl(m, 2) * RWI[2][1];
      tmp[m][5] += kl(m, 0) * RWI[0][2] + kl(m, 1) * RWI[1][2] + kl(m, 2) * RWI[2][2];
    }

    tmp[m][6] = kl(m, 6) * R(0,0) + kl(m, 7) * R(0,1) + kl(m, 8) * R(0,2);
    tmp[m][7] = kl(m, 6) * R(1,0) + kl(m, 7) * R(1,1) + kl(m, 8) * R(1,2);
    tmp[m][8] = kl(m, 6) * R(2,0) + kl(m, 7) * R(2,1) + kl(m, 8) * R(2,2);

    tmp[m][9]  = kl(m, 9) * R(0,0) + kl(m, 10) * R(0,1) + kl(m, 11) * R(0,2);
    tmp[m][10] = kl(m, 9) * R(1,0) + kl(m, 10) * R(1,1) + kl(m, 11) * R(1,2);
    tmp[m][11] = kl(m, 9) * R(2,0) + kl(m, 10) * R(2,1) + kl(m, 11) * R(2,2);

    if (nodeJOffset) {
      tmp[m][ 9] += kl(m, 6) * RWJ[0][0] + kl(m, 7) * RWJ[1][0] + kl(m, 8) * RWJ[2][0];
      tmp[m][10] += kl(m, 6) * RWJ[0][1] + kl(m, 7) * RWJ[1][1] + kl(m, 8) * RWJ[2][1];
      tmp[m][11] += kl(m, 6) * RWJ[0][2] + kl(m, 7) * RWJ[1][2] + kl(m, 8) * RWJ[2][2];
    }
  }

  // Now compute T'_{lg}*(kl*T_{lg})
  for (int m = 0; m < 12; m++) {
    kg(0, m) = R(0,0) * tmp[0][m] + R(0,1) * tmp[1][m] + R(0,2) * tmp[2][m];
    kg(1, m) = R(1,0) * tmp[0][m] + R(1,1) * tmp[1][m] + R(1,2) * tmp[2][m];
    kg(2, m) = R(2,0) * tmp[0][m] + R(2,1) * tmp[1][m] + R(2,2) * tmp[2][m];

    kg(3, m) = R(0,0) * tmp[3][m] + R(0,1) * tmp[4][m] + R(0,2) * tmp[5][m];
    kg(4, m) = R(1,0) * tmp[3][m] + R(1,1) * tmp[4][m] + R(1,2) * tmp[5][m];
    kg(5, m) = R(2,0) * tmp[3][m] + R(2,1) * tmp[4][m] + R(2,2) * tmp[5][m];

    if (nodeIOffset) {
      kg(3, m) += RWI[0][0] * tmp[0][m] + RWI[1][0] * tmp[1][m] + RWI[2][0] * tmp[2][m];
      kg(4, m) += RWI[0][1] * tmp[0][m] + RWI[1][1] * tmp[1][m] + RWI[2][1] * tmp[2][m];
      kg(5, m) += RWI[0][2] * tmp[0][m] + RWI[1][2] * tmp[1][m] + RWI[2][2] * tmp[2][m];
    }

    kg( 6, m) = R(0,0) * tmp[6][m] + R(0,1) * tmp[7][m] + R(0,2) * tmp[8][m];
    kg( 7, m) = R(1,0) * tmp[6][m] + R(1,1) * tmp[7][m] + R(1,2) * tmp[8][m];
    kg( 8, m) = R(2,0) * tmp[6][m] + R(2,1) * tmp[7][m] + R(2,2) * tmp[8][m];

    kg( 9, m) = R(0,0) * tmp[9][m] + R(0,1) * tmp[10][m] + R(0,2) * tmp[11][m];
    kg(10, m) = R(1,0) * tmp[9][m] + R(1,1) * tmp[10][m] + R(1,2) * tmp[11][m];
    kg(11, m) = R(2,0) * tmp[9][m] + R(2,1) * tmp[10][m] + R(2,2) * tmp[11][m];

    if (nodeJOffset) {
      kg( 9, m) += RWJ[0][0] * tmp[6][m] + RWJ[1][0] * tmp[7][m] + RWJ[2][0] * tmp[8][m];
      kg(10, m) += RWJ[0][1] * tmp[6][m] + RWJ[1][1] * tmp[7][m] + RWJ[2][1] * tmp[8][m];
      kg(11, m) += RWJ[0][2] * tmp[6][m] + RWJ[1][2] * tmp[7][m] + RWJ[2][2] * tmp[8][m];
    }
  }
  return kg;
}



const Vector &
LinearFrameTransf3d::getGlobalResistingForce(const Vector &pb, const Vector &p0)
{
  // transform resisting forces from the basic system to local coordinates
  static VectorND<12> pl;

  double q0 = pb(0);
  double q1 = pb(1);
  double q2 = pb(2);
  double q3 = pb(3);
  double q4 = pb(4);
  double q5 = pb(5);

  double oneOverL = 1.0 / L;

  pl[0]  = -q0;                    // Ni
  pl[1]  =  oneOverL * (q1 + q2);  // Viy
  pl[2]  = -oneOverL * (q3 + q4);  // Viz
  pl[3]  = -q5;                    // Ti
  pl[4]  =  q3;
  pl[5]  =  q1;
  pl[6]  =  q0;                    // Nj
  pl[7]  = -pl[1];                 // Vjy
  pl[8]  = -pl[2];                 // Vjz
  pl[9]  = q5;                     // Tj
  pl[10] = q4;
  pl[11] = q2;


  pl[0] += p0[0];
  pl[1] += p0[1];
  pl[7] += p0[2];
  pl[2] += p0[3];
  pl[8] += p0[4];

  static VectorND<12> pg;
  static Vector wrapper(pg);

  pg  = pushResponse(pl);

  return wrapper;
}


const Matrix &
LinearFrameTransf3d::getGlobalStiffMatrix(const Matrix &KB, const Vector &pb)
{
  static double kb[6][6];     // Basic stiffness
  static MatrixND<12,12> kl;  // Local stiffness
  static double tmp[12][12];  // Temporary storage
  const  double oneOverL = 1.0 / L;

  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
      kb[i][j] = KB(i, j);

  // Transform basic stiffness to local system
  // First compute kb*T_{bl}
  for (int i = 0; i < 6; i++) {
    tmp[i][ 0] = -kb[i][0];
    tmp[i][ 1] =  oneOverL * (kb[i][1] + kb[i][2]);
    tmp[i][ 2] = -oneOverL * (kb[i][3] + kb[i][4]);
    tmp[i][ 3] = -kb[i][5];
    tmp[i][ 4] =  kb[i][3];
    tmp[i][ 5] =  kb[i][1];
    tmp[i][ 6] =  kb[i][0];
    tmp[i][ 7] = -tmp[i][1];
    tmp[i][ 8] = -tmp[i][2];
    tmp[i][ 9] =  kb[i][5];
    tmp[i][10] =  kb[i][4];
    tmp[i][11] =  kb[i][2];
  }

  // Now compute T'_{bl}*(kb*T_{bl})
  for (int i = 0; i < 12; i++) {
    kl( 0, i) = -tmp[0][i];
    kl( 1, i) =  oneOverL * (tmp[1][i] + tmp[2][i]);
    kl( 2, i) = -oneOverL * (tmp[3][i] + tmp[4][i]);
    kl( 3, i) = -tmp[5][i];
    kl( 4, i) = tmp[3][i];
    kl( 5, i) = tmp[1][i];
    kl( 6, i) = tmp[0][i];
    kl( 7, i) = -kl(1, i);
    kl( 8, i) = -kl(2, i);
    kl( 9, i) = tmp[5][i];
    kl(10, i) = tmp[4][i];
    kl(11, i) = tmp[2][i];
  }


  static MatrixND<12,12> Kg;
  static Matrix wrapper(Kg);
  Kg = pushConstant(kl);
  return wrapper;

}

const Matrix &
LinearFrameTransf3d::getInitialGlobalStiffMatrix(const Matrix &KB)
{
  static double kb[6][6];     // Basic stiffness
  static MatrixND<12,12> kl;  // Local stiffness
  static double tmp[12][12];  // Temporary storage
  double oneOverL = 1.0 / L;

  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
      kb[i][j] = KB(i, j);

  // Transform basic stiffness to local system
  // First compute kb*T_{bl}
  for (int i = 0; i < 6; i++) {
    tmp[i][0]  = -kb[i][0];
    tmp[i][1]  =  oneOverL * (kb[i][1] + kb[i][2]);
    tmp[i][2]  = -oneOverL * (kb[i][3] + kb[i][4]);
    tmp[i][3]  = -kb[i][5];
    tmp[i][4]  =  kb[i][3];
    tmp[i][5]  =  kb[i][1];
    tmp[i][6]  =  kb[i][0];
    tmp[i][7]  = -tmp[i][1];
    tmp[i][8]  = -tmp[i][2];
    tmp[i][9]  =  kb[i][5];
    tmp[i][10] =  kb[i][4];
    tmp[i][11] =  kb[i][2];
  }

  // Now compute T'_{bl}*(kb*T_{bl})
  for (int i = 0; i < 12; i++) {
    kl( 0, i) = -tmp[0][i];
    kl( 1, i) =  oneOverL * (tmp[1][i] + tmp[2][i]);
    kl( 2, i) = -oneOverL * (tmp[3][i] + tmp[4][i]);
    kl( 3, i) = -tmp[5][i];
    kl( 4, i) = tmp[3][i];
    kl( 5, i) = tmp[1][i];
    kl( 6, i) = tmp[0][i];
    kl( 7, i) = -kl(1, i);
    kl( 8, i) = -kl(2, i);
    kl( 9, i) = tmp[5][i];
    kl(10, i) = tmp[4][i];
    kl(11, i) = tmp[2][i];
  }

  static MatrixND<12,12> kg;
  static Matrix M(kg);

  kg = pushConstant(kl);

  return M;
}

FrameTransform3d *
LinearFrameTransf3d::getCopy()
{
  // create a new instance of LinearFrameTransf3d

  LinearFrameTransf3d *theCopy = nullptr;

  static Vector xz(3);
  xz(0) = R(0,2);
  xz(1) = R(1,2);
  xz(2) = R(2,2);

  Vector offsetI(3);
  Vector offsetJ(3);

  if (nodeIOffset) {
    offsetI(0) = nodeIOffset[0];
    offsetI(1) = nodeIOffset[1];
    offsetI(2) = nodeIOffset[2];
  }

  if (nodeJOffset) {
    offsetJ(0) = nodeJOffset[0];
    offsetJ(1) = nodeJOffset[1];
    offsetJ(2) = nodeJOffset[2];
  }

  theCopy = new LinearFrameTransf3d(this->getTag(), xz, offsetI, offsetJ);

  theCopy->nodeIPtr = nodeIPtr;
  theCopy->nodeJPtr = nodeJPtr;
  theCopy->L        = L;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      theCopy->R(j,i) = R(j,i);

  return theCopy;
}

const Matrix &
LinearFrameTransf3d::getGlobalMatrixFromLocal(const Matrix &ml)
{
  Matrix3D Rm;
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      Rm(i,j) = R(i,j);

  blk3x12x3(Rm, ml, kg);

//static MatrixND<12,12> Kg;
//static Matrix M(Kg);
//Kg = pushConstant(kl);

  return kg;
}

const Vector &
LinearFrameTransf3d::getPointGlobalCoordFromLocal(const Vector &xl)
{
  static Vector xg(3);
  return xg;
}

const Vector &
LinearFrameTransf3d::getPointGlobalDisplFromBasic(double xi, const Vector &uxb)
{
  static Vector uxg(3);
  return uxg;
}

const Vector &
LinearFrameTransf3d::getPointLocalDisplFromBasic(double xi, const Vector &uxb)
{
  // compute displacements at point xi, in local coordinates
  static Vector uxl(3);
  return uxl;
}


// Sensitivity

bool
LinearFrameTransf3d::isShapeSensitivity()
{
  int nodeParameterI = nodeIPtr->getCrdsSensitivity();
  int nodeParameterJ = nodeJPtr->getCrdsSensitivity();
  // TODO: implement dvz

  return (nodeParameterI != 0 || nodeParameterJ != 0);
}


double
LinearFrameTransf3d::getLengthGrad()
{
  const int di = nodeIPtr->getCrdsSensitivity();
  const int dj = nodeJPtr->getCrdsSensitivity();

  Vector3D dxi{0.0};
  Vector3D dxj{0.0};

  if (di != 0)
    dxi(di-1) = 1.0;
  if (dj != 0)
    dxj(dj-1) = 1.0;

  return 1/L*(xj - xi).dot(dxj - dxi);
}

double
LinearFrameTransf3d::getd1overLdh()
{
  return -getLengthGrad()/(L*L);
}

const Vector &
LinearFrameTransf3d::getGlobalResistingForceShapeSensitivity(const Vector &pb,
                                                           const Vector &p0,
                                                           int gradNumber)
{
  
  static Vector pg(12);
  pg.Zero();

  //
  // dp = T_{lg}' pl
  //
  int dv = 0; // TODO: Sensitivity
  int di = nodeIPtr->getCrdsSensitivity();
  int dj = nodeJPtr->getCrdsSensitivity();

  VectorND<12> pl = pushLocal(pb, L);
  
  pl[0] += p0[0];
  pl[1] += p0[1];
  pl[7] += p0[2];
  pl[2] += p0[3];
  pl[8] += p0[4];

  Matrix3D dR = FrameOrientationGradient(xi, xj, vz, di, dj, dv);
  for (int i=0; i<4; i++)
    for (int j=0; j<3; j++)
      pg(i*3+j) = dR(j,0) * pl(3*i) + dR(j,1) * pl(3*i+1) + dR(j,2) * pl(3*i+2);

  
  //
  // dp += T_{gl} dpl
  //
  double dL = this->getLengthGrad();
  double doneOverL = -dL/(L*L);
  VectorND<12> dpl{0.0};
  dpl[1]  =  doneOverL * (pb[1] + pb[2]);  // Viy
  dpl[2]  = -doneOverL * (pb[3] + pb[4]);  // Viz
  dpl[7]  = -dpl[1];                       // Vjy
  dpl[8]  = -dpl[2];                       // Vjz

  for (int i=0; i<4; i++)
    for (int j=0; j<3; j++)
      pg(i*3+j) += R(j,0) * dpl(3*i) + R(j,1) * dpl(3*i+1) + R(j,2) * dpl(3*i+2);

  return pg;
}


const Vector &
LinearFrameTransf3d::getBasicDisplFixedGrad()
{ 
  //
  // Form ug
  //
  const Vector &disp1 = nodeIPtr->getTrialDisp();
  const Vector &disp2 = nodeJPtr->getTrialDisp();

  VectorND<12> ug;
  for (int i = 0; i < 6; i++) {
    ug[i]     = disp1(i);
    ug[i + 6] = disp2(i);
  }

  if (nodeIInitialDisp != 0) {
    for (int j = 0; j < 6; j++)
      ug[j] -= nodeIInitialDisp[j];
  }

  if (nodeJInitialDisp != 0) {
    for (int j = 0; j < 6; j++)
      ug[j + 6] -= nodeJInitialDisp[j];
  }

  //
  // dub += (T_{bl}' T_{lg} + T_{bl} T_{lg}') * ug
  //
  int dv = 0; // TODO
  int di = nodeIPtr->getCrdsSensitivity();
  int dj = nodeJPtr->getCrdsSensitivity();

  static VectorND<6> dub;
  static Vector wrapper(dub);

//Matrix3D dR = FrameOrientationGradient(xi, xj, vz, di, dj, dv);

  VectorND<12> ul = getLocal<2,6>(ug, R, nodeIOffset, nodeJOffset);
  dub = getBasic(ul, 1/L);

  //
  //
  //
  dub[0] += 0;
  double dL = this->getLengthGrad();
  double doneOverL = -dL/(L*L);
  double tmp   = doneOverL * (ul[1] - ul[7]);
  dub[1] +=  tmp;
  dub[2] +=  tmp;
  tmp   = doneOverL * (ul[8] - ul[2]);
  dub[3] +=  tmp;
  dub[4] +=  tmp;

  return wrapper;

}

const Vector &
LinearFrameTransf3d::getBasicDisplTotalGrad(int gradNumber)
{

  VectorND<12> dug;
  for (int i = 0; i < 6; i++) {
    dug[i]     = nodeIPtr->getDispSensitivity((i + 1), gradNumber);
    dug[i + 6] = nodeJPtr->getDispSensitivity((i + 1), gradNumber);
  }

  static VectorND<6> dub;
  static Vector wrapper(dub);

  // dub = T_{bl} T_{lg} * ug'

  VectorND<12> dul = getLocal<2,6>(dug, R, nodeIOffset, nodeJOffset);
  dub = getBasic(dul, 1/L);

  wrapper += getBasicDisplFixedGrad();

  return wrapper;
}


void
LinearFrameTransf3d::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": \"" << this->getTag()
      << "\", \"type\": \"LinearFrameTransf3d\"";
    s << ", \"vecInLocXZPlane\": [" 
      << R(0,2) << ", " 
      << R(1,2) << ", "
      << R(2,2) << "]";
    if (nodeIOffset != 0)
      s << ", \"iOffset\": [" << nodeIOffset[0] << ", " << nodeIOffset[1]
        << ", " << nodeIOffset[2] << "]";
    if (nodeJOffset != 0)
      s << ", \"jOffset\": [" << nodeJOffset[0] << ", " << nodeJOffset[1]
        << ", " << nodeJOffset[2] << "]";
    s << "}";

    return;
  }

  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "\nFrameTransform: " << this->getTag() << " Type: LinearFrameTransf3d\n";
    s << "\tOrientation: " << Matrix(&R(0,0), 3,3) << "\n";
    if (nodeIOffset)
      s << "\tNode I offset: " << nodeIOffset[0] << " " << nodeIOffset[1] << " "
        << nodeIOffset[2] << "\n";
    if (nodeJOffset)
      s << "\tNode J offset: " << nodeJOffset[0] << " " << nodeJOffset[1] << " "
        << nodeJOffset[2] << "\n";
  }
}


int
LinearFrameTransf3d::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;

  static Vector data(23);
  data(0) = this->getTag();
  data(1) = L;

  if (nodeIOffset != 0) {
    data(2) = nodeIOffset[0];
    data(3) = nodeIOffset[1];
    data(4) = nodeIOffset[2];
  } else {
    data(2) = 0.0;
    data(3) = 0.0;
    data(4) = 0.0;
  }

  if (nodeJOffset != 0) {
    data(5) = nodeJOffset[0];
    data(6) = nodeJOffset[1];
    data(7) = nodeJOffset[2];
  } else {
    data(5) = 0.0;
    data(6) = 0.0;
    data(7) = 0.0;
  }

  if (nodeIInitialDisp != 0) {
    data(8)  = nodeIInitialDisp[0];
    data(9)  = nodeIInitialDisp[1];
    data(10) = nodeIInitialDisp[2];
    data(11) = nodeIInitialDisp[3];
    data(12) = nodeIInitialDisp[4];
    data(13) = nodeIInitialDisp[5];
  } else {
    data(8)  = 0.0;
    data(9)  = 0.0;
    data(10) = 0.0;
    data(11) = 0.0;
    data(12) = 0.0;
    data(13) = 0.0;
  }

  if (nodeJInitialDisp != 0) {
    data(14) = nodeJInitialDisp[0];
    data(15) = nodeJInitialDisp[1];
    data(16) = nodeJInitialDisp[2];
    data(17) = nodeJInitialDisp[3];
    data(18) = nodeJInitialDisp[4];
    data(19) = nodeJInitialDisp[5];
  } else {
    data(14) = 0.0;
    data(15) = 0.0;
    data(16) = 0.0;
    data(17) = 0.0;
    data(18) = 0.0;
    data(19) = 0.0;
  }

  data(20) = R(0,2);
  data(21) = R(1,2);
  data(22) = R(2,2);

  res += theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) {
    opserr << "LinearFrameTransf3d::sendSelf - failed to send Vector\n";

    return res;
  }

  return res;
}

int
LinearFrameTransf3d::recvSelf(int cTag, Channel &theChannel,
                            FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static Vector data(23);

  res += theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0) {
    opserr << "LinearFrameTransf3d::recvSelf - failed to receive Vector\n";

    return res;
  }

  this->setTag((int)data(0));
  L       = data(1);
  data(0) = this->getTag();
  data(1) = L;

  int flag;
  int i, j;

  flag = 0;
  for (int i = 2; i <= 4; i++)
    if (data(i) != 0.0)
      flag = 1;
  if (flag == 1) {
    if (nodeIOffset == 0)
      nodeIOffset = new double[3];
    for (i = 2, j = 0; i <= 4; i++, j++)
      nodeIOffset[j] = data(i);
  }

  flag = 0;
  for (i = 5; i <= 7; i++)
    if (data(i) != 0.0)
      flag = 1;
  if (flag == 1) {
    if (nodeJOffset == 0)
      nodeJOffset = new double[3];
    for (int i = 5, j = 0; i <= 7; i++, j++)
      nodeJOffset[j] = data(i);
  }

  flag = 0;
  for (int i = 8; i <= 13; i++)
    if (data(i) != 0.0)
      flag = 1;
  if (flag == 1) {
    if (nodeIInitialDisp == 0)
      nodeIInitialDisp = new double[6];
    for (i = 8, j = 0; i <= 13; i++, j++)
      nodeIInitialDisp[j] = data(i);
  }

  flag = 0;
  for (int i = 14; i <= 19; i++)
    if (data(i) != 0.0)
      flag = 1;
  if (flag == 1) {
    if (nodeJInitialDisp == 0)
      nodeJInitialDisp = new double[6];
    for (i = 14, j = 0; i <= 19; i++, j++)
      nodeJInitialDisp[j] = data(i);
  }

  R(0,2) = data(20);
  R(1,2) = data(21);
  R(2,2) = data(22);

  initialDispChecked = true;
  return res;
}
