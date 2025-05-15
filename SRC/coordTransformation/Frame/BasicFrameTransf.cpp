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
  Vector3D wi = t.getNodeRotationLogarithm(0),
           wj = t.getNodeRotationLogarithm(1);
  ub[0] = t.getNodePosition(1)[0];
  ub[1] = wi[2];
  ub[2] = wj[2];
  ub[3] = wi[1];
  ub[4] = wj[1];
  ub[5] = wj[0] - wi[0];
  return wrapper;
}

const Vector &
BasicFrameTransf3d::getBasicIncrDeltaDisp()
{
  constexpr static int ndf = 6;
  static VectorND<6> ub;
  static Vector wrapper(ub);
  VectorND<12> ul = t.getStateVariation();
  ub[0] =  ul[1*ndf+0];
  ub[1] =  ul[0*ndf+5];
  ub[2] =  ul[1*ndf+5];
  ub[3] =  ul[0*ndf+4];
  ub[4] =  ul[1*ndf+4];
  ub[5] =  ul[1*ndf+3] - ul[0*ndf+3];
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
BasicFrameTransf3d::getGlobalResistingForce(const Vector &q_pres, const Vector &p0)
{
  // transform resisting forces from the basic system to local coordinates
  
  VectorND<NDF*2> pl{};
  #ifdef DO_BASIC
  double L = theCoordTransf->getInitialLength();
  const double q1 = q_pres[imz],
                q2 = q_pres[jmz],
                q3 = q_pres[imy],
                q4 = q_pres[jmy];
  pl[0*NDF+1]  =  (q1 + q2)/L;      // Viy
  pl[0*NDF+2]  = -(q3 + q4)/L;      // Viz
  pl[1*NDF+1]  = -pl[1];            // Vjy
  pl[1*NDF+2]  = -pl[2];            // Vjz
  #endif
  pl[0*NDF+0]  = -q_pres[jnx];      // Ni
  pl[0*NDF+3]  = -q_pres[jmx];      // Ti
  pl[0*NDF+4]  =  q_pres[imy];
  pl[0*NDF+5]  =  q_pres[imz];
  pl[1*NDF+0]  =  q_pres[jnx];      // Nj
  pl[1*NDF+3]  =  q_pres[jmx];      // Tj
  pl[1*NDF+4]  =  q_pres[jmy];
  pl[1*NDF+5]  =  q_pres[jmz];

  VectorND<NDF*2> pf;
  pf.zero();
  pf[0*NDF + 0] = p0[0];
  pf[0*NDF + 1] = p0[1];
  pf[0*NDF + 2] = p0[3];
  pf[1*NDF + 1] = p0[2];
  pf[1*NDF + 2] = p0[4];

  static VectorND<NDF*2> pg;
  static Vector wrapper(pg);

  pg  = t.pushResponse(pl);
  pg += t.pushConstant(pf);
  return wrapper;
}


const Matrix &
BasicFrameTransf3d::getGlobalStiffMatrix(const Matrix &kb, const Vector &q_pres)
{

  VectorND<NDF*2> pl{};
#ifdef DO_BASIC
  double L = t.getInitialLength();
  const double q1 = q_pres[imz],
               q2 = q_pres[jmz],
               q3 = q_pres[imy],
               q4 = q_pres[jmy];
  pl[0*NDF+1]  =  (q1 + q2)/L;      // Viy
  pl[0*NDF+2]  = -(q3 + q4)/L;      // Viz
  pl[1*NDF+1]  = -pl[1];            // Vjy
  pl[1*NDF+2]  = -pl[2];            // Vjz
#endif
  pl[0*NDF+0]  = -q_pres[jnx];      // Ni
  pl[0*NDF+3]  = -q_pres[jmx];      // Ti
  pl[0*NDF+4]  =  q_pres[imy];
  pl[0*NDF+5]  =  q_pres[imz];
  pl[1*NDF+0]  =  q_pres[jnx];      // Nj
  pl[1*NDF+3]  =  q_pres[jmx];      // Tj
  pl[1*NDF+4]  =  q_pres[jmy];
  pl[1*NDF+5]  =  q_pres[jmz];
  

  MatrixND<2*NDF,2*NDF> kl;
  kl.zero();
#ifndef DO_BASIC
  for (int i=0; i<NDF*2; i++) {
    int ii = std::abs(iq[i]);
    double c = 1.0;
    if (ii >= NBV)
      continue;

    for (int j=0; j<NDF*2; j++) {
      int jj = std::abs(iq[j]);
      if (jj >= NBV)
        continue;

      kl(i,j) = kb(ii, jj)*c;
    }
  }

  for (int i = 0; i < 2*NDF; i++) {
    kl(0*NDF+0, i) = kl(i, 0*NDF+0) = i==0? kl(NDF+0, NDF+0): (i==3? kl(NDF+0, NDF+3) : -kl( NDF+0, i));
    kl(0*NDF+3, i) = kl(i, 0*NDF+3) = i==0? kl(NDF+3, NDF+0): (i==3? kl(NDF+3, NDF+3) : -kl( NDF+3, i));
  }

#else
  // Transform basic stiffness to local system
  static double tmp[NBV][12];
  // First compute kb*T_{bl}
  for (int i = 0; i < 6; i++) {
    tmp[i][ 0] = -kb(i, 0);
    tmp[i][0*NDF+1] =  (kb(i, 1) + kb(i, 2))/L;
    tmp[i][0*NDF+2] = -(kb(i, 3) + kb(i, 4))/L;
    tmp[i][1*NDF+1] = -tmp[i][1];
    tmp[i][1*NDF+2] = -tmp[i][2];
    tmp[i][ 3] = -kb(i, 5);
    tmp[i][ 4] =  kb(i, 3);
    tmp[i][ 5] =  kb(i, 1);
    tmp[i][ 6] =  kb(i, 0);
    tmp[i][ 9] =  kb(i, 5);
    tmp[i][10] =  kb(i, 4);
    tmp[i][11] =  kb(i, 2);
  }

  kl.zero();
  // Now compute T'_{bl}*(kb*T_{bl})
  for (int i = 0; i < 12; i++) {
    kl(0*NDF+0, i) =  -tmp[0][i];
    kl(0*NDF+1, i) =  (tmp[1][i] + tmp[2][i])/L;
    kl(0*NDF+2, i) = -(tmp[3][i] + tmp[4][i])/L;
    kl(1*NDF+1, i) = -kl(1, i);
    kl(1*NDF+2, i) = -kl(2, i);
    kl(0*NDF+3, i) =  -tmp[5][i]; // 
    kl(0*NDF+4, i) =   tmp[3][i];
    kl(0*NDF+5, i) =   tmp[1][i];
    kl(1*NDF+0, i) =   tmp[0][i];
    kl(1*NDF+3, i) =   tmp[5][i];
    kl(1*NDF+4, i) =   tmp[4][i];
    kl(1*NDF+5, i) =   tmp[2][i];
  }
#endif

  static MatrixND<12,12> Kg;
  static Matrix Wrapper(Kg);

  Kg = t.pushResponse(kl, pl);

  return Wrapper;
}

const Matrix &
BasicFrameTransf3d::getInitialGlobalStiffMatrix(const Matrix &KB)
{
  static double kb[6][6];     // Basic stiffness
  static MatrixND<12,12> kl;  // Local stiffness
  double tmp[6][12]{};   // Temporary storage

  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
      kb[i][j] = KB(i, j);

  // Transform basic stiffness to local system
  // First compute kb*T_{bl}
  for (int i = 0; i < 6; i++) {
    tmp[i][0]  = -kb[i][0];
#ifdef DO_BASIC
double oneOverL = 1.0 / t.getInitialLength();
    tmp[i][1]  =  oneOverL * (kb[i][1] + kb[i][2]);
    tmp[i][2]  = -oneOverL * (kb[i][3] + kb[i][4]);
    tmp[i][7]  = -tmp[i][1];
    tmp[i][8]  = -tmp[i][2];
#endif
    tmp[i][3]  = -kb[i][5];
    tmp[i][4]  =  kb[i][3];
    tmp[i][5]  =  kb[i][1];
    tmp[i][6]  =  kb[i][0];
    tmp[i][9]  =  kb[i][5];
    tmp[i][10] =  kb[i][4];
    tmp[i][11] =  kb[i][2];
  }

  // Now compute T'_{bl}*(kb*T_{bl})
  for (int i = 0; i < 12; i++) {
    kl( 0, i) = -tmp[0][i];
#ifdef DO_BASIC
    kl( 1, i) =  oneOverL * (tmp[1][i] + tmp[2][i]);
    kl( 2, i) = -oneOverL * (tmp[3][i] + tmp[4][i]);
    kl( 7, i) = -kl(1, i);
    kl( 8, i) = -kl(2, i);
#endif
    kl( 3, i) = -tmp[5][i];
    kl( 4, i) = tmp[3][i];
    kl( 5, i) = tmp[1][i];
    kl( 6, i) = tmp[0][i];
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
  // return t.getGlobalResistingForceShapeSensitivity(pb, p0, gradNumber);

  static VectorND<6> dub;
  static Vector wrapper(dub);
  opserr << "WARNING unimplemented method\n";
  return wrapper;
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
