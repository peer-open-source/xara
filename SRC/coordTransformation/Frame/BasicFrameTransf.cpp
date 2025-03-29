//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
//
#include <Vector.h>
#include <Matrix.h>
#include <Matrix3D.h>
#include <Node.h>
#include <Channel.h>
#include <Logging.h>
#include <BasicFrameTransf.h>
#include "blk3x12x3.h"

using namespace OpenSees;


BasicFrameTransf3d::BasicFrameTransf3d(FrameTransform<2,6> *t)

  : FrameTransform3d(t->getTag(), 0),
    t(*t)
{
  
}

BasicFrameTransf3d::~BasicFrameTransf3d()
{
  delete &t;
}

int
BasicFrameTransf3d::commitState()
{
  return t.commit();
}

int
BasicFrameTransf3d::revertToLastCommit()
{
  return t.revertToLastCommit();
}

int
BasicFrameTransf3d::revertToStart()
{
  return t.revertToStart();
}

int
BasicFrameTransf3d::update()
{
  return t.update();
}


int
BasicFrameTransf3d::initialize(Node *i, Node *j)
{
  std::array<Node*, 2> nodes = {i, j};
  return t.initialize(nodes);
}

int
BasicFrameTransf3d::getLocalAxes(Vector &XAxis, Vector &YAxis, Vector &ZAxis)
{
  Vector3D x, y, z;
  int s = t.getLocalAxes(x, y, z);
  for (int i=0; i<3; i++) {
    XAxis(i) = x[i];
    YAxis(i) = y[i];
    ZAxis(i) = z[i];
  }
  return s;
}

double
BasicFrameTransf3d::getInitialLength()
{
  return t.getInitialLength();
}

double
BasicFrameTransf3d::getDeformedLength()
{
  return t.getDeformedLength();
}


const Vector &
BasicFrameTransf3d::getBasicTrialDisp()
{
  static VectorND<6> ub;
  static Vector wrapper(ub);
  VectorND<12> ul;
  ul.insert<0, 3>(t.getNodePosition(0));
  ul.insert<3, 3>(t.getNodeRotationLogarithm(0));
  ul.insert<6, 3>(t.getNodePosition(1));
  ul.insert<9, 3>(t.getNodeRotationLogarithm(1));
  ub = getBasic(ul, 1/t.getInitialLength());
  return wrapper;
}

const Vector &
BasicFrameTransf3d::getBasicIncrDeltaDisp()
{
  static VectorND<6> ub;
  static Vector wrapper(ub);
  VectorND<12> ul = t.getStateVariation();
  ub = getBasic(ul, 1/t.getInitialLength());
  return wrapper;
}

const Vector &
BasicFrameTransf3d::getBasicIncrDisp()
{
  static VectorND<6> ub;
  static Vector wrapper(ub);
  opserr << "Unimplemented method\n";
  return wrapper;
}

const Vector &
BasicFrameTransf3d::getBasicTrialVel()
{
  static VectorND<6> ub;
  static Vector wrapper(ub);
  opserr << "Unimplemented method\n";
  return wrapper;
}

VectorND<12>
BasicFrameTransf3d::pushResponse(VectorND<12>&pl)
{
  return t.pushResponse(pl);
}

MatrixND<12,12>
BasicFrameTransf3d::pushResponse(MatrixND<12,12>&kl, const VectorND<12>& pl)
{
  return t.pushResponse(kl, pl);
}

VectorND<12>
BasicFrameTransf3d::pushConstant(const VectorND<12>&pl) const
{
  return t.pushConstant(pl);
}

MatrixND<12,12>
BasicFrameTransf3d::pushConstant(const MatrixND<12,12>& kl)
{
  return t.pushConstant(kl);
}



const Vector &
BasicFrameTransf3d::getGlobalResistingForce(const Vector &pb, const Vector &p0)
{
  // transform resisting forces from the basic system to local coordinates
  static VectorND<12> pl;

  double q0 = pb(0);
  double q1 = pb(1);
  double q2 = pb(2);
  double q3 = pb(3);
  double q4 = pb(4);
  double q5 = pb(5);

  double oneOverL = 1.0 / t.getInitialLength();

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
BasicFrameTransf3d::getGlobalStiffMatrix(const Matrix &KB, const Vector &pb)
{
  static double kb[6][6];     // Basic stiffness
  static MatrixND<12,12> kl;  // Local stiffness
  static double tmp[12][12];  // Temporary storage
  const  double oneOverL = 1.0 / t.getInitialLength();

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

  static VectorND<12> pl;

  double q0 = pb(0);
  double q1 = pb(1);
  double q2 = pb(2);
  double q3 = pb(3);
  double q4 = pb(4);
  double q5 = pb(5);

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

  static MatrixND<12,12> Kg;
  static Matrix wrapper(Kg);
  Kg = pushResponse(kl, pl);
  return wrapper;

}

const Matrix &
BasicFrameTransf3d::getInitialGlobalStiffMatrix(const Matrix &KB)
{
  static double kb[6][6];     // Basic stiffness
  static MatrixND<12,12> kl;  // Local stiffness
  static double tmp[6][12];  // Temporary storage
  double oneOverL = 1.0 / t.getInitialLength();

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
BasicFrameTransf3d::getCopy()
{

  BasicFrameTransf3d *theCopy = new BasicFrameTransf3d(t.getCopy());
  return theCopy;
}

const Matrix &
BasicFrameTransf3d::getGlobalMatrixFromLocal(const Matrix &M)
{
  //
  // Do diag(R)*M*diag(R)'
  //

  static MatrixND<12,12> Kout;
  static Matrix wrapper(Kout);
  wrapper = M;
  MatrixND<12,12> Kg = pushConstant(Kout);
  Kout = Kg;
  return wrapper;
}

const Vector &
BasicFrameTransf3d::getPointGlobalCoordFromLocal(const Vector &xl)
{
  static Vector xg(3);
  return xg;
}

const Vector &
BasicFrameTransf3d::getPointGlobalDisplFromBasic(double xi, const Vector &uxb)
{
  static Vector uxg(3);
  return uxg;
}

const Vector &
BasicFrameTransf3d::getPointLocalDisplFromBasic(double xi, const Vector &uxb)
{
  static Vector uxl(3);
  return uxl;
}

//
// Sensitivity
//
bool
BasicFrameTransf3d::isShapeSensitivity()
{
  return t.isShapeSensitivity();
}


double
BasicFrameTransf3d::getLengthGrad()
{
  return t.getLengthGrad();
}

double
BasicFrameTransf3d::getd1overLdh()
{
  double L = t.getInitialLength();
  return -getLengthGrad()/(L*L);
}

const Vector &
BasicFrameTransf3d::getGlobalResistingForceShapeSensitivity(const Vector &pb,
                                                           const Vector &p0,
                                                           int gradNumber)
{
  return t.getGlobalResistingForceShapeSensitivity(pb, p0, gradNumber);
}


const Vector &
BasicFrameTransf3d::getBasicDisplFixedGrad()
{
  static VectorND<6> dub;
  static Vector wrapper(dub);
  opserr << "WARNING unimplemented method\n";
  return wrapper;

}

const Vector &
BasicFrameTransf3d::getBasicDisplTotalGrad(int gradNumber)
{

  static VectorND<6> dub;
  static Vector wrapper(dub);
  opserr << "WARNING unimplemented method\n";
  return wrapper;
}


void
BasicFrameTransf3d::Print(OPS_Stream &s, int flag)
{
  t.Print(s, flag);
}


int
BasicFrameTransf3d::sendSelf(int cTag, Channel &theChannel)
{
  return -1;
}

int
BasicFrameTransf3d::recvSelf(int cTag, Channel &theChannel,
                            FEM_ObjectBroker &theBroker)
{
  return -1;
}
