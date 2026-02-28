//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
//
// TODO:
//  - Move formation of km from setDomain; it needs to be re-formed
//    after updateParameter
//
// Written: cc,cmp 05/2024
//
#include <Frame/Prism.h>
#include <Vector.h>
#include <VectorND.h>
#include <Matrix.h>
#include <MatrixND.h>
#include <PrismFrame2d.h>
#include <ElementalLoad.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <CrdTransf.h>
#include <FrameSection.h>
#include <Information.h>
#include <Parameter.h>
#include <ElementResponse.h>
#include <math.h>

using namespace OpenSees;

Matrix PrismFrame2d::K(6,6);
Vector PrismFrame2d::P(6);

PrismFrame2d::PrismFrame2d()
  : Element(0,ELE_TAG_ElasticBeam2d), 
    A(0.0), E(0.0), Iz(0.0), G(0.0), 
    alpha(0.0), depth(0.0), rho(0.0), 
    mass_flag(0), 
    release(0),
    Q(6),
    connectedExternalNodes(2), 
    theCoordTransf(nullptr)
{
  q.zero();
  q0.zero();
  p0.zero();
  kg.zero();

  for (int i=0; i<2; i++)
    theNodes[i] = nullptr;      
}

PrismFrame2d::PrismFrame2d(int tag, double a, double e, double i, 
                           int Nd1, int Nd2, CrdTransf &coordTransf,
                           double Alpha, double depth_, double r, int cm,
                           int rel, int geom_flag_)
  :Element(tag,ELE_TAG_ElasticBeam2d), 
   A(a), E(e), G(0.0), Iz(i), Ay(0.0), 
   alpha(Alpha), depth(depth_), rho(r), 
   mass_flag(cm), release(rel), 
   geom_flag(geom_flag_),
   shear_flag(0),
   Q(6), connectedExternalNodes(2), theCoordTransf(nullptr)
{
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;
    
  theCoordTransf = coordTransf.getCopy2d();

  // Make no release if input not 0, 1, 2, or 3
  if (release < 0 || release > 3)
    release = 0;
  
  q.zero();
  q0.zero();
  p0.zero();
  kg.zero();

  // set node pointers to NULL
  theNodes[0] = nullptr;
  theNodes[1] = nullptr;
}

PrismFrame2d::PrismFrame2d(int tag, int Nd1, int Nd2, 
                           FrameSection &section,  
                           CrdTransf &coordTransf, 
                           double Alpha, double depth_, 
                           double r, int cm, bool use_mass, int rel,
                           int geom_flag_, int shear_flag)
  : Element(tag,ELE_TAG_ElasticBeam2d), 
    G(0.0), Ay(0.0),
    alpha(Alpha), depth(depth_), 
    rho(r), mass_flag(cm), 
    release(rel),
    geom_flag(geom_flag_),
    shear_flag(shear_flag),
    Q(6), connectedExternalNodes(2), theCoordTransf(nullptr)
{
#if 0
  section.getIntegral(Field::Unit,   State::Init, A);
  section.getIntegral(Field::UnitY,  State::Init, Ay);
  section.getIntegral(Field::UnitYY, State::Init, Iz);


  const Matrix &sectTangent = section.getInitialTangent();
  const ID &sectCode = section.getType();
  for (int i=0; i<sectCode.Size(); i++) {
    int code = sectCode(i);
    switch(code) {
    case SECTION_RESPONSE_P:
      E = sectTangent(i,i)/A;
      break;
    case SECTION_RESPONSE_VY:
      if (Ay != 0.0)
        G  = sectTangent(i,i)/Ay;
      break;
    default:
      break;
    }
  }
  
#else
  Frame::Prism shape_data(2,3);
  section.getShape(shape_data);
  A  =  shape_data.A;
  Iz = *shape_data.Iz;
  E  = *shape_data.E;
  G  = *shape_data.G;
  if (!shear_flag) {
    Ay = 0.0;
  } 
  else {
    Ay = *shape_data.Ay;
  }
#endif
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;
  
  theCoordTransf = coordTransf.getCopy2d();

  // Make no release if input not 0, 1, 2, or 3
  if (release < 0 || release > 3)
    release = 0;
  
  q0.zero();
  p0.zero();

  // set node pointers to NULL
  theNodes[0] = nullptr;
  theNodes[1] = nullptr;
}

PrismFrame2d::~PrismFrame2d()
{
  if (theCoordTransf)
    delete theCoordTransf;
}

int
PrismFrame2d::getNumExternalNodes() const
{
  return 2;
}

const ID &
PrismFrame2d::getExternalNodes() 
{
  return connectedExternalNodes;
}

Node **
PrismFrame2d::getNodePtrs() 
{
  return theNodes;
}

int
PrismFrame2d::getNumDOF()
{
  return 6;
}

void
PrismFrame2d::setDomain(Domain *theDomain)
{
  if (theDomain == nullptr) {
    opserr << "PrismFrame2d::setDomain -- Domain is null\n";
    return;
  }
    
  theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
  theNodes[1] = theDomain->getNode(connectedExternalNodes(1));    
  
  if (theNodes[0] == nullptr) {
    opserr << "PrismFrame2d::setDomain -- Node 1: " << connectedExternalNodes(0) << " does not exist\n";
    return;
  }
                            
  if (theNodes[1] == nullptr) {
    opserr << "PrismFrame2d::setDomain -- Node 2: " << connectedExternalNodes(1) << " does not exist\n";
    return;
  }

  int dofNd1 = theNodes[0]->getNumberDOF();
  int dofNd2 = theNodes[1]->getNumberDOF();    

      
  if (theDomain != nullptr)
    this->Element::link(*theDomain);
  
  if (theCoordTransf->initialize(theNodes[0], theNodes[1]) != 0) {
    opserr << "PrismFrame2d::setDomain -- Error initializing coordinate transformation\n";
    return;
  }
  
  L = theCoordTransf->getInitialLength();

  if (L == 0.0)
    opserr << "PrismFrame2d::setDomain -- Element has zero length\n";

  if (G != 0 && Ay != 0)
    phi = 12.0 * E * Iz / (L * L * G * Ay);
  else
    phi = 0.0;

  formBasicStiffness(km);
  ke = km; // 
}

inline void
PrismFrame2d::formBasicStiffness(OpenSees::MatrixND<3,3>& kb) const
{
  //
  // Compute basic material stiffness matrix
  //
  kb.zero();

  kb(0,0) = E*A/L;

  double EI = E*Iz;
  if (release == 0) {
    kb(1,1) = kb(2,2) = EI*(4+phi)/(L*(1+phi));
    kb(2,1) = kb(1,2) = EI*(2-phi)/(L*(1+phi));    
  }
  if (release == 1) { // release I
    kb(2,2) = 3.0*EI/L;
  }
  if (release == 2) { // release J
    kb(1,1) = 3.0*EI/L;
  }
  
  double r[3] = {(double)(release&0b100), 
                 (double)(release&0b001),
                 (double)(release&0b010)};

  MatrixND<3,3> ah {{
     1-r[0],        0,                  0           ,
       0   ,      1-r[1],       -0.5*(1-r[2])*r[1]  ,
       0   ,-0.5*(1-r[1])*r[2],       1-r[2]        }};
}

int
PrismFrame2d::commitState()
{
  int retVal = 0;
  // call base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "PrismFrame2d::commitState () - failed in base class";
  }
  retVal += theCoordTransf->commitState();
  return retVal;
}

int
PrismFrame2d::revertToLastCommit()
{
  return theCoordTransf->revertToLastCommit();
}

int
PrismFrame2d::revertToStart()
{
  return theCoordTransf->revertToStart();
}

int
PrismFrame2d::update()
{
  int ok = theCoordTransf->update();

  const Vector &v = theCoordTransf->getBasicTrialDisp();
  L = theCoordTransf->getInitialLength();

  double N = E*A/L*v[0];

  switch (geom_flag) {
    case 0:
    break;

    case 1:
      kg.zero();
      kg(1,1) = kg(2,2) =  4.0*N*L/30.0;
      kg(1,2) = kg(2,1) = -1.0*N*L/30.0;
      ke = km + kg;

    case 2:
    {
      // initialize with linear coefficients
      double A = 4*E*Iz/L,
             B = 2*E*Iz/L;
      double psi = L*std::sqrt(std::fabs(N)/(E*Iz));
      if (N < -1.e-4) {
        // Axial force is compressive
        double cs = std::cos(psi);
        double sn = std::sin(psi);
        double C  = psi*E*Iz/L/(2.0*(1 - cs) - psi*sn);
        A   = C*(sn - psi*cs);
        B   = C*(psi - sn);

      } else if (N > 1.e-4) {
        // Axial force is tensile
        double snh = std::sinh(psi);
        double csh = std::cosh(psi);
        double C   = psi*E*Iz/L/(2*(csh - 1) - psi*snh);
        A   = C*(snh - psi*csh);
        B   = C*(psi - snh);

      } else {
      // if axial force is near zero keep linear coefficients
      
      }

      if (std::fabs(N) > 1e-8) {
        ke(1,1) = ke(2,2) =  A;
        ke(1,2) = ke(2,1) =  B;
      } else {
        ke = km;
      }
      break;
    }
  }

#if 0
  // Account for releases
  if (release != 0) {
    switch (release) {
      case 0:
        k(1,1) = 0;
      case 1:
        k(3,3)   = k(2,2) - std::pow(k(2,3),2)/k(2,2);
        k(2,2:3) = zeros(1,2);
        k(3,2)   = 0;
      case 2:
        k(2,2)   = k(3,3)-k(2,3)^2/k(3,3);
        k(3,2:3) = zeros(1,2);
        k(2,3)   = 0;
    }
  }
#endif

  q = ke*v;
  q += q0;

  return ok;
}

const Matrix &
PrismFrame2d::getTangentStiff()
{
  return theCoordTransf->getGlobalStiffMatrix(Matrix(ke), q);
}

const Matrix &
PrismFrame2d::getInitialStiff()
{
  return theCoordTransf->getInitialGlobalStiffMatrix(Matrix(km));
}

const Matrix &
PrismFrame2d::getMass()
{
  bool mass_initialized;
  if (rho > 0.0)  {
    // get initial element length
    if (mass_flag == 0)  {
        static MatrixND<6,6> M;
        static Matrix Wrapper(M);
        M.zero();
        // lumped mass matrix
        double m = 0.5*rho*L;
        M(0,0) = M(1,1) = M(3,3) = M(4,4) = m;
        return Wrapper;

    } else {

      // consistent mass matrix
      Matrix mlTrn(6, 6), mlRot(6, 6), ml(6, 6);
      mlTrn.Zero();
      mlRot.Zero();
      ml.Zero();
      double c1x  = rho * L / 210.0;
      mlTrn(0, 0) = mlTrn(3, 3) = c1x * 70.0;
      mlTrn(0, 3) = mlTrn(3, 0) = c1x * 35.0;
      double c1z                = c1x / pow(1.0 + phi, 2);
      mlTrn(1, 1) = mlTrn(4, 4) = c1z * (70.0 * phi * phi + 147.0 * phi + 78.0);
      mlTrn(1, 4) = mlTrn(4, 1) = c1z * (35.0 * phi * phi + 63.0 * phi + 27.0);
      mlTrn(2, 2) = mlTrn(5, 5) = c1z * L * L / 4.0 * (7.0 * phi * phi + 14.0 * phi + 8.0);
      mlTrn(2, 5) = mlTrn(5, 2) = -c1z * L * L / 4.0 * (7.0 * phi * phi + 14.0 * phi + 6.0);
      mlTrn(1, 2) = mlTrn(2, 1) = c1z * L / 4.0 * (35.0 * phi * phi + 77.0 * phi + 44.0);
      mlTrn(4, 5) = mlTrn(5, 4) = -mlTrn(1, 2);
      mlTrn(1, 5) = mlTrn(5, 1) = -c1z * L / 4.0 * (35.0 * phi * phi + 63.0 * phi + 26.0);
      mlTrn(2, 4) = mlTrn(4, 2) = -mlTrn(1, 5);
      double c2z                = rho / A * Iz / (30.0 * L * pow(1.0 + phi, 2));
      mlRot(1, 1) = mlRot(4, 4) = c2z * 36.0;
      mlRot(1, 4) = mlRot(4, 1) = -mlRot(1, 1);
      mlRot(2, 2) = mlRot(5, 5) = c2z * L * L * (10.0 * phi * phi + 5.0 * phi + 4.0);
      mlRot(2, 5) = mlRot(5, 2) = c2z * L * L * (5.0 * phi * phi - 5.0 * phi - 1.0);
      mlRot(1, 2) = mlRot(2, 1) = mlRot(1, 5) = mlRot(5, 1) = -c2z * L * (15.0 * phi - 3.0);
      mlRot(2, 4) = mlRot(4, 2) = mlRot(4, 5) = mlRot(5, 4) = -mlRot(1, 2);

      // add translational and rotational parts
      ml = mlTrn + mlRot;
      // transform from local to global system
      return theCoordTransf->getGlobalMatrixFromLocal(Matrix(M));
    }

    mass_initialized = true;
  }

  return K;
}

void 
PrismFrame2d::zeroLoad()
{
  Q.Zero();
  q0.zero();
  p0.zero();
  return;
}

int 
PrismFrame2d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);

  if (type == LOAD_TAG_Beam2dUniformLoad) {
    double wt = data(0)*loadFactor;  // Transverse (+ve upward)
    double wa = data(1)*loadFactor;  // Axial (+ve from node I to J)

    double V = 0.5*wt*L;
    double P = wa*L;

    // Reactions in basic system
    p0[0] -= P;
    p0[1] -= V;
    p0[2] -= V;

    // Fixed end forces in basic system
    q0[0] -= 0.5*P;
    if (release == 0) {
      double M = V*L/6.0; // wt*L*L/12
      q0[1] -= M;
      q0[2] += M;
    }
    if (release == 1) {
      q0[2] += wt*L*L/8;
    }
    if (release == 2) {
      q0[1] -= wt*L*L/8;
    }
    if (release == 3) {
      // Nothing to do
    }
  }

  else if (type == LOAD_TAG_Beam2dPartialUniformLoad) {
    // These equations should works for partial trapezoidal load
    double waa = data(2)*loadFactor;  // Axial
    double wab = data(3)*loadFactor;  // Axial
    double wya = data(0)*loadFactor;  // Transverse
    double wyb = data(1)*loadFactor;  // Transverse
    double a   = data(4)*L;
    double b   = data(5)*L;

        // auxiliary values
    double ba = b-a;
    double ba2 = pow(b, 2.0) - pow(a, 2.0);
    double ba3 = pow(b, 3.0) - pow(a, 3.0);
    double ba4 = pow(b, 4.0) - pow(a, 4.0);
    double ba5 = pow(b, 5.0) - pow(a, 5.0);
    double z1 = wya + (wya*a)/ba - (wyb*a)/ba;
    double wybpa = wya+wyb;
    double wybma = wyb-wya;
    double L2 = pow(L, 2.0);

    // equivalent nodal forces
    double Fyt = 0.5*wybpa*ba;
    double V2 = (1.0/L)*(wya*ba*(a+0.5*ba)+0.5*wybma*ba*(a+(2.0/3.0)*ba));
    double V1 = Fyt-V2;
    double M1 = (0.5*z1*ba2) + (wybma*ba3/(3.0*ba)) 
              - (z1*ba3*2.0/(3.0*L))  - (wybma*ba4/(2.0*L*ba)) 
              + (z1*ba4/(4.0*L2))     + (wybma*ba5/(5.0*L2*ba));
    double M2 = (-1.0*z1*ba3/(3.0*L)) - (wybma*ba4/(4.0*L*ba)) 
              + (z1*ba4/(4.0*L2))     + (wybma*ba5/(5.0*L2*ba));
    double P = waa*ba + 0.5*(wab-waa)*ba;
    double PJ = (1.0/L)*(waa*ba*(a + 0.5*ba) + 0.5*(wab-waa)*ba*(a + (2.0/3.0)*ba));

    // Reactions in basic system
    p0[0] -=  P;
    p0[1] -= V1;
    p0[2] -= V2;

    // Fixed end forces in basic system
    q0[0] -= PJ;
    q0[1] -= M1;
    q0[2] -= M2;
  }

  else if (type == LOAD_TAG_Beam2dPointLoad) {
    double P = data(0)*loadFactor;
    double N = data(1)*loadFactor;
    double aOverL = data(2);

    if (aOverL < 0.0 || aOverL > 1.0)
      return 0;

    double a = aOverL*L;
    double b = L-a;

    // Reactions in basic system
    p0[0] -= N;
    double V1 = P*(1.0-aOverL);
    double V2 = P*aOverL;
    p0[1] -= V1;
    p0[2] -= V2;

    double L2 = 1.0/(L*L);
    double a2 = a*a;
    double b2 = b*b;

    // Fixed end forces in basic system
    q0[0] -= N*aOverL;
    if (release == 0) {
      double M1 = -a * b2 * P * L2;
      double M2 = a2 * b * P * L2;
      q0[1] += M1;
      q0[2] += M2;
    }
    if (release == 1) {
      q0[2] += 0.5*P*a*b*L2*(a+L);
    }
    if (release == 2) {
      q0[1] -= 0.5*P*a*b*L2*(b+L);
    }
    if (release == 3) {
      // Nothing to do
    }    
  }
  
  else if (type == LOAD_TAG_Beam2dTempLoad) {
    double Ttop1 = data(0)* loadFactor;
    double Tbot1 = data(1)* loadFactor;
    double Ttop2 = data(2)* loadFactor;
    double Tbot2 = data(3)* loadFactor;
        
    // fixed end forces due to a linear thermal load
    double dT1 = Ttop1-Tbot1;
    double dT = (Ttop2-Tbot2)-(Ttop1-Tbot1);
    double a = alpha/depth;  // constant based on temp difference at top and bottom, 
    // coefficient of thermal expansion and beam depth
    double M1 = a*E*Iz*(-dT1+(4.0/3.0)*dT); //Fixed End Moment end 1
    double M2 = a*E*Iz*(dT1+(5.0/3.0)*dT); //Fixed End Moment end 2
    double F = alpha*(((Ttop2+Ttop1)/2+(Tbot2+Tbot1)/2)/2)*E*A; // Fixed End Axial Force
    double M1M2divL =(M1+M2)/L; // Fixed End Shear
    
    // Reactions in basic system
    p0[0] += 0.0;
    p0[1] += M1M2divL;
    p0[2] -= M1M2divL;

    // Fixed end forces in basic system
    q0[0] -= F;
    q0[1] += M1;
    q0[2] += M2;
  }

  else {
    opserr << "PrismFrame2d::addLoad()  -- load type unknown for element with tag: " << this->getTag() << "\n";
    return -1;
  }

  return 0;
}

int
PrismFrame2d::addInertiaLoadToUnbalance(const Vector &accel)
{
  if (rho == 0.0)
    return 0;

  // get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
        
  if (3 != Raccel1.Size() || 3 != Raccel2.Size()) {
    opserr << "PrismFrame2d::addInertiaLoadToUnbalance matrix and vector sizes are incompatable\n";
    return -1;
  }
    
  // want to add ( - fact * M R * accel ) to unbalance
  if (mass_flag == 0)  {
    // take advantage of lumped mass matrix
    double m = 0.5*rho*L;

    Q(0) -= m * Raccel1(0);
    Q(1) -= m * Raccel1(1);

    Q(3) -= m * Raccel2(0);
    Q(4) -= m * Raccel2(1);
  } else  {
    // use matrix vector multip. for consistent mass matrix
    static Vector Raccel(6);
    for (int i=0; i<3; i++)  {
      Raccel(i)   = Raccel1(i);
      Raccel(i+3) = Raccel2(i);
    }
    Q.addMatrixVector(1.0, this->getMass(), Raccel, -1.0);
  }
  
  return 0;
}

const Vector &
PrismFrame2d::getResistingForceIncInertia()
{        
  P = this->getResistingForce();
  
  // subtract external load P = P - Q
  P.addVector(1.0, Q, -1.0);
  
  // add the damping forces if rayleigh damping
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    P.addVector(1.0, this->getRayleighDampingForces(), 1.0);
    
  if (rho == 0.0)
    return P;

  // add inertia forces from element mass
  const Vector &accel1 = theNodes[0]->getTrialAccel();
  const Vector &accel2 = theNodes[1]->getTrialAccel();    
  
  if (mass_flag == 0)  {
    // take advantage of lumped mass matrix
    double m = 0.5*rho*L;

    P(0) += m * accel1(0);
    P(1) += m * accel1(1);

    P(3) += m * accel2(0);
    P(4) += m * accel2(1);
  } else  {
    // use matrix vector multip. for consistent mass matrix
    static Vector accel(6);
    for (int i=0; i<3; i++)  {
      accel(i)   = accel1(i);
      accel(i+3) = accel2(i);
    }
    P.addMatrixVector(1.0, this->getMass(), accel, 1.0);
  }
  
  return P;
}


const Vector &
PrismFrame2d::getResistingForce()
{
  theCoordTransf->update();
  
  const Vector &v = theCoordTransf->getBasicTrialDisp();
  q  = ke*v;
  q += q0;

  // Vector for reactions in basic system
  Vector p0Vec(p0);
  
  P = theCoordTransf->getGlobalResistingForce(q, p0Vec);

  return P;
}

int
PrismFrame2d::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;

    static Vector data(17);
    
    data(0) = A;
    data(1) = E; 
    data(2) = Iz; 
    data(3) = rho;
    data(4) = mass_flag;
    data(5) = this->getTag();
    data(6) = connectedExternalNodes(0);
    data(7) = connectedExternalNodes(1);
    data(8) = theCoordTransf->getClassTag();
            
    int dbTag = theCoordTransf->getDbTag();
    
    if (dbTag == 0) {
      dbTag = theChannel.getDbTag();
      if (dbTag != 0)
        theCoordTransf->setDbTag(dbTag);
    }

    data(9) = dbTag;
    data(10) = alpha;
    data(11) = depth;

    data(12) = alphaM;
    data(13) = betaK;
    data(14) = betaK0;
    data(15) = betaKc;
    data(16) = release;
    
    // Send the data vector
    res += theChannel.sendVector(this->getDbTag(), cTag, data);
    if (res < 0) {
      opserr << "PrismFrame2d::sendSelf -- could not send data Vector\n";
      return res;
    }

    // Ask the CoordTransf to send itself
    res += theCoordTransf->sendSelf(cTag, theChannel);
    if (res < 0) {
      opserr << "PrismFrame2d::sendSelf -- could not send CoordTransf\n";
      return res;
    }
    
    return res;
}

int
PrismFrame2d::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
      
  static Vector data(17);

  res += theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0) {
    opserr << "PrismFrame2d::recvSelf -- could not receive data Vector\n";
    return res;
  }

  A = data(0);
  E = data(1);
  Iz = data(2); 
  alpha = data(10);
  depth = data(11);

  alphaM = data(12);
  betaK  = data(13);
  betaK0 = data(14);
  betaKc = data(15);
  release = (int)data(16);
  
  rho = data(3);
  mass_flag = (int)data(4);
  this->setTag((int)data(5));
  connectedExternalNodes(0) = (int)data(6);
  connectedExternalNodes(1) = (int)data(7);

  // Check if the CoordTransf is null; if so, get a new one
  int crdTag = (int)data(8);
  if (theCoordTransf == 0) {
    theCoordTransf = theBroker.getNewCrdTransf(crdTag);
    if (theCoordTransf == 0) {
      opserr << "PrismFrame2d::recvSelf -- could not get a CrdTransf2d\n";
      exit(-1);
    }
  }
  
  // Check that the CoordTransf is of the right type; if not, delete
  // the current one and get a new one of the right type
  if (theCoordTransf->getClassTag() != crdTag) {
    delete theCoordTransf;
    theCoordTransf = theBroker.getNewCrdTransf(crdTag);
    if (theCoordTransf == 0) {
      opserr << "PrismFrame2d::recvSelf -- could not get a CrdTransf2d\n";
      exit(-1);
    }
  }
      
  // Now, receive the CoordTransf
  theCoordTransf->setDbTag((int)data(9));
  res += theCoordTransf->recvSelf(cTag, theChannel, theBroker);
  if (res < 0) {
    opserr << "PrismFrame2d::recvSelf -- could not receive CoordTransf\n";
    return res;
  }
  
  return res;
}

void
PrismFrame2d::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_ELEM_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"PrismFrame2d\", ";
    s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
    s << "\"E\": " << E << ", ";
    s << "\"G\": " << G << ", ";
    s << "\"A\": " << A << ", ";
    s << "\"Ay\": " << Ay << ", ";
    s << "\"Iz\": "<< Iz << ", ";
    s << "\"massperlength\": "<< rho << ", ";
    s << "\"release\": "<< release << ", ";
    s << "\"kinematics\": "<< geom_flag << ", ";
    s << "\"mass_flag\": "<< mass_flag << ", ";
    s << "\"crdTransformation\": \"" << theCoordTransf->getTag() << "\"}";
    return;
  }

  if (flag == -1) {
    int eleTag = this->getTag();
    s << "EL_BEAM\t" << eleTag << "\t";
    s << 0 << "\t" << 0 << "\t" << connectedExternalNodes(0) << "\t" << connectedExternalNodes(1) ;
    s << "0\t0.0000000\n";
  }

  if (flag == OPS_PRINT_CURRENTSTATE) {
    this->getResistingForce();
    s << "\nPrismFrame2d: " << this->getTag() << "\n";
    s << "\tConnected Nodes: " << connectedExternalNodes ;
    s << "\tCoordTransf: " << theCoordTransf->getTag() << "\n";
    s << "\tmass density:  " << rho << ", mass_flag: " << mass_flag << "\n";
    s << "\trelease code:  " << release << "\n";
    double P  = q[0];
    double M1 = q[1];
    double M2 = q[2];
    double V = (M1+M2)/L;
    s << "\tEnd 1 Forces (P V M): " << -P+p0[0]
      << " " << V+p0[1] << " " << M1 << "\n";
    s << "\tEnd 2 Forces (P V M): " << P
      << " " << -V+p0[2] << " " << M2 << "\n";
  }

}


Response*
PrismFrame2d::setResponse(const char **argv, int argc, OPS_Stream &output)
{

  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType","PrismFrame2d");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);

    // global forces
  if (strcmp(argv[0],"force") == 0 || 
      strcmp(argv[0],"forces") == 0 ||
      strcmp(argv[0],"globalForce") == 0 || 
      strcmp(argv[0],"globalForces") == 0) {

    output.tag("ResponseType","Px_1");
    output.tag("ResponseType","Py_1");
    output.tag("ResponseType","Mz_1");
    output.tag("ResponseType","Px_2");
    output.tag("ResponseType","Py_2");
    output.tag("ResponseType","Mz_2");

    theResponse =  new ElementResponse(this, 2, Vector(6));
  
  // local forces
  } else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0) {

    output.tag("ResponseType","N_1");
    output.tag("ResponseType","V_1");
    output.tag("ResponseType","M_1");
    output.tag("ResponseType","N_2");
    output.tag("ResponseType","V_2");
    output.tag("ResponseType","M_2");
    
    theResponse = new ElementResponse(this, 3, Vector(6));

  // basic forces
  } else if (strcmp(argv[0],"basicForce") == 0 || strcmp(argv[0],"basicForces") == 0) {

    output.tag("ResponseType","N");
    output.tag("ResponseType","M_1");
    output.tag("ResponseType","M_2");
    
    theResponse = new ElementResponse(this, 4, Vector(3));

  // deformations
  }  else if (strcmp(argv[0],"deformatons") == 0 || 
              strcmp(argv[0],"basicDeformations") == 0 ||
              strcmp(argv[0],"basicDeformation") == 0) {
    
    output.tag("ResponseType","eps");
    output.tag("ResponseType","theta1");
    output.tag("ResponseType","theta2");
    theResponse = new ElementResponse(this, 5, Vector(3));
  
  // chord rotation -
  } else if (strcmp(argv[0],"chordRotation") == 0 || strcmp(argv[0],"chordDeformation") == 0) {

    output.tag("ResponseType","eps");
    output.tag("ResponseType","theta1");
    output.tag("ResponseType","theta2");

    theResponse =  new ElementResponse(this, 5, Vector(3));
  }

  if (theResponse == 0)
    theResponse = theCoordTransf->setResponse(argv, argc, output);

  output.endTag(); // ElementOutput
  
  return theResponse;
}

int
PrismFrame2d::getResponse (int responseID, Information &info)
{
  this->getResistingForce();

  switch (responseID) {
  case 1: // stiffness
    return info.setMatrix(this->getTangentStiff());

  case 2: // global forces
    return info.setVector(this->getResistingForce());
    
  case 3: // local forces
  {
    // Axial
    double N = q[0];
    P(3) =  N;
    P(0) = -N+p0[0];
    // Moment
    double M1 = q[1];
    double M2 = q[2];
    P(2) = M1;
    P(5) = M2;
    // Shear
    double V = (M1+M2)/L;
    P(1) =  V + p0[1];
    P(4) = -V + p0[2];
    return info.setVector(P);
  }
    
  case 4: // basic forces
    return info.setVector(q);

  case 5:
    return info.setVector(theCoordTransf->getBasicTrialDisp());

  default:
    return -1;
  }
}

int
PrismFrame2d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  // E of the beam interior
  if (strcmp(argv[0],"E") == 0) {
    param.setValue(E);
    return param.addObject(1, this);
  }
  // A of the beam interior
  if (strcmp(argv[0],"A") == 0) {
    param.setValue(A);
    return param.addObject(2, this);
  }
  // I of the beam interior
  if (strcmp(argv[0],"I") == 0) {
    param.setValue(Iz);
    return param.addObject(3, this);
  }
  // mass per length
  if (strcmp(argv[0],"rho") == 0) {
    param.setValue(rho);
    return param.addObject(4, this);
  }
  // moment release
  if (strcmp(argv[0],"release") == 0) {
    param.setValue(release);
    return param.addObject(5, this);
  }  
  
  return -1;
}

int
PrismFrame2d::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
  case -1:
    return -1;
  case 1:
    E  = info.theDouble;
    return 0;
  case 2:
    A  = info.theDouble;
    return 0;
  case 3:
    Iz = info.theDouble;
    return 0;
  case 4:
    rho = info.theDouble;
    return 0;
  case 5:
    release = (int)info.theDouble;
    if (release < 0 || release > 3)
      release = 0;
    return 0;
  default:
    return -1;
  }

  formBasicStiffness(km);
  ke = km; // 
}

