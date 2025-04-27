/* ****************************************************************** **
**    Openers - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
//
// Purpose: This file contains the class definition for ElasticBeam2d.
// ElasticBeam2d is a 3d beam element. As such it can only
// connect to a node with 6-dof. 
//
// Written: fmk 11/95
// Revised:
//
#include <ElasticBeam2d.h>
#include <ElementalLoad.h>

#include <Domain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <CrdTransf.h>
#include <SectionForceDeformation.h>
#include <Information.h>
#include <Parameter.h>
#include <ElementResponse.h>

#include <math.h>
#include <stdlib.h>
#include <elementAPI.h>
#include <string>
#include <ElementIter.h>

#include <map>

Matrix ElasticBeam2d::K(6,6);
Vector ElasticBeam2d::P(6);
// Matrix ElasticBeam2d::kb(3,3);

ElasticBeam2d::ElasticBeam2d()
  : Element(0,ELE_TAG_ElasticBeam2d), 
    A(0.0), E(0.0), I(0.0), alpha(0.0), depth(0.0), rho(0.0), cMass(0), release(0),
    Q(6), connectedExternalNodes(2), theCoordTransf(0)
{
  // does nothing
  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;

  // set node pointers to NULL
  for (int i=0; i<2; i++)
    theNodes[i] = 0;      
}

ElasticBeam2d::ElasticBeam2d(int tag, double a, double e, double i, 
                             int Nd1, int Nd2, CrdTransf &coordTransf,
                             double Alpha, double depth_, double r, int cm,
                             int rel)
  :Element(tag,ELE_TAG_ElasticBeam2d), 
   A(a), E(e), I(i), alpha(Alpha), depth(depth_), rho(r), cMass(cm), release(rel),
  Q(6), connectedExternalNodes(2), theCoordTransf(0)
{
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;
    
  theCoordTransf = coordTransf.getCopy2d();
    
  if (!theCoordTransf) {
    opserr << "ElasticBeam2d::ElasticBeam2d -- failed to get copy of coordinate transformation\n";
    exit(01);
  }

  // Make no release if input not 0, 1, 2, or 3
  if (release < 0 || release > 3)
    release = 0;
  
  q0[0] = 0.0;
  q0[1] = 0.0;
  q0[2] = 0.0;

  p0[0] = 0.0;
  p0[1] = 0.0;
  p0[2] = 0.0;

  // set node pointers to NULL
  theNodes[0] = 0;
  theNodes[1] = 0;
}

ElasticBeam2d::ElasticBeam2d(int tag, int Nd1, int Nd2, SectionForceDeformation &section,  
                             CrdTransf &coordTransf, double Alpha, double depth_, double r, int cm, int rel)
  :Element(tag,ELE_TAG_ElasticBeam2d), alpha(Alpha), depth(depth_), rho(r), cMass(cm), release(rel),
  Q(6), connectedExternalNodes(2), theCoordTransf(nullptr)
{
  E = 1.0;
  rho = r;
  cMass = cm;

  const Matrix &sectTangent = section.getInitialTangent();
  const ID &sectCode = section.getType();
  for (int i=0; i<sectCode.Size(); i++) {
    int code = sectCode(i);
    switch(code) {
    case SECTION_RESPONSE_P:
      A = sectTangent(i,i);
      break;
    case SECTION_RESPONSE_MZ:
      I = sectTangent(i,i);
      break;
    default:
      break;
    }
  }
  
  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;
  
  theCoordTransf = coordTransf.getCopy2d();
  
  if (!theCoordTransf) {
    opserr << "ElasticBeam2d::ElasticBeam2d -- failed to get copy of coordinate transformation\n";
    exit(-1);
  }

  // Make no release if input not 0, 1, 2, or 3
  if (release < 0 || release > 3)
    release = 0;
  
  q0.zero();
  p0.zero();

  // set node pointers to NULL
  theNodes[0] = nullptr;
  theNodes[1] = nullptr;
}

ElasticBeam2d::~ElasticBeam2d()
{
    if (theCoordTransf)
        delete theCoordTransf;
}

int
ElasticBeam2d::getNumExternalNodes() const
{
    return 2;
}

const ID &
ElasticBeam2d::getExternalNodes() 
{
    return connectedExternalNodes;
}

Node **
ElasticBeam2d::getNodePtrs() 
{
  return theNodes;
}

int
ElasticBeam2d::getNumDOF()
{
    return 6;
}

void
ElasticBeam2d::setDomain(Domain *theDomain)
{
    if (theDomain == nullptr) {
      opserr << "ElasticBeam2d::setDomain -- Domain is null\n";
      return;
    }
    
    theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
    theNodes[1] = theDomain->getNode(connectedExternalNodes(1));    
    
    if (theNodes[0] == nullptr) {
      opserr << "ElasticBeam2d::setDomain -- Node 1: " << connectedExternalNodes(0) << " does not exist\n";
      return;
    }
                              
    if (theNodes[1] == nullptr) {
      opserr << "ElasticBeam2d::setDomain -- Node 2: " << connectedExternalNodes(1) << " does not exist\n";
      return;
    }

    int dofNd1 = theNodes[0]->getNumberDOF();
    int dofNd2 = theNodes[1]->getNumberDOF();    
    
    if (dofNd1 != 3) {
      opserr << "ElasticBeam2d::setDomain -- Node 1: " << connectedExternalNodes(0) 
             << " has incorrect number of DOF\n";
      exit(-1);
    }
    
    if (dofNd2 != 3) {
      opserr << "ElasticBeam2d::setDomain -- Node 2: " << connectedExternalNodes(1) 
             << " has incorrect number of DOF\n";
      exit(-1);
    }
        
    this->DomainComponent::setDomain(theDomain);
    
    if (theCoordTransf->initialize(theNodes[0], theNodes[1]) != 0) {
        opserr << "ElasticBeam2d::setDomain -- Error initializing coordinate transformation\n";
        exit(-1);
    }
    
    double L = theCoordTransf->getInitialLength();

    if (L == 0.0) {
      opserr << "ElasticBeam2d::setDomain -- Element has zero length\n";
      exit(-1);
    }

    //
    // Compute basic stiffness matrix
    //
    kb.zero();

    kb(0,0) = E*A/L;  

    double EI     = E*I;
    if (release == 0) {
      kb(1,1) = kb(2,2) = 4.0*EI/L;
      kb(2,1) = kb(1,2) = 2.0*EI/L;    
    }
    if (release == 1) { // release I
      kb(2,2) = 3.0*EI/L;
    }
    if (release == 2) { // release J
      kb(1,1) = 3.0*EI/L;
    }
}

int
ElasticBeam2d::commitState()
{
  int retVal = 0;
  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "ElasticBeam2d::commitState () - failed in base class";
  }    
  retVal += theCoordTransf->commitState();
  return retVal;
}

int
ElasticBeam2d::revertToLastCommit()
{
  return theCoordTransf->revertToLastCommit();
}

int
ElasticBeam2d::revertToStart()
{
  return theCoordTransf->revertToStart();
}

int
ElasticBeam2d::update()
{
  return theCoordTransf->update();
}

const Matrix &
ElasticBeam2d::getTangentStiff()
{
  const Vector &v = theCoordTransf->getBasicTrialDisp();

  q = kb*v;
  q += q0;

#if 0
  double L = theCoordTransf->getInitialLength();
  double EoverL   = E/L;
  double EAoverL  = A*EoverL;                        // EA/L
  // determine q = kv + q0
  q[0] = EAoverL*v(0);
  kb.zero();
  kb(0,0) = EAoverL;  
  if (release == 0) {
    double EIoverL2 = 2.0*I*EoverL;                // 2EI/L
    double EIoverL4 = 2.0*EIoverL2;                // 4EI/L
    q[1] = EIoverL4*v(1) + EIoverL2*v(2);
    q[2] = EIoverL2*v(1) + EIoverL4*v(2);
    kb(1,1) = kb(2,2) = EIoverL4;
    kb(2,1) = kb(1,2) = EIoverL2;    
  }
  if (release == 1) { // release I
    q[1] = 0.0;
    double EIoverL3 = 3.0*I*EoverL;
    q[2] = EIoverL3*v(2);
    kb(2,2) = EIoverL3;
  }
  if (release == 2) { // release J
    q[2] = 0.0;
    double EIoverL3 = 3.0*I*EoverL;
    q[1] = EIoverL3*v(1);
    kb(1,1) = EIoverL3;
  }
  if (release == 3) { // both I and J
    q[1] = 0.0;
    q[2] = 0.0;
  }
  q[0] += q0[0];
  q[1] += q0[1];
  q[2] += q0[2];
#endif  
  return theCoordTransf->getGlobalStiffMatrix(kb, q);
}

const Matrix &
ElasticBeam2d::getInitialStiff()
{
#if 0
  double L = theCoordTransf->getInitialLength();

  double EoverL   = E/L;
  double EAoverL  = A*EoverL;                        // EA/L

  kb.zero();
  kb(0,0) = EAoverL;
  if (release == 0) {
    double EIoverL2 = 2.0*I*EoverL;                // 2EI/L
    double EIoverL4 = 2.0*EIoverL2;                // 4EI/L
    kb(1,1) = kb(2,2) = EIoverL4;
    kb(2,1) = kb(1,2) = EIoverL2;
  }
  if (release == 1) { // release I
    kb(2,2) = 3.0*I*EoverL;
  }
  if (release == 2) { // release J
    kb(1,1) = 3.0*I*EoverL;
  }
#endif
  return theCoordTransf->getInitialGlobalStiffMatrix(kb);
}

const Matrix &
ElasticBeam2d::getMass()
{ 
    K.Zero();
    
    if (rho > 0.0)  {
        // get initial element length
        double L = theCoordTransf->getInitialLength();
        if (cMass == 0)  {
            // lumped mass matrix
            double m = 0.5*rho*L;
            K(0,0) = K(1,1) = K(3,3) = K(4,4) = m;

        } else  {
            // consistent mass matrix
            static Matrix ml(6,6);
            double m = rho*L/420.0;
            ml(0,0) = ml(3,3) = m*140.0;
            ml(0,3) = ml(3,0) = m*70.0;

            ml(1,1) = ml(4,4) = m*156.0;
            ml(1,4) = ml(4,1) = m*54.0;
            ml(2,2) = ml(5,5) = m*4.0*L*L;
            ml(2,5) = ml(5,2) = -m*3.0*L*L;
            ml(1,2) = ml(2,1) = m*22.0*L;
            ml(4,5) = ml(5,4) = -ml(1,2);
            ml(1,5) = ml(5,1) = -m*13.0*L;
            ml(2,4) = ml(4,2) = -ml(1,5);
            
            // transform local mass matrix to global system
            K = theCoordTransf->getGlobalMatrixFromLocal(ml);
        }
    }
    
    return K;
}

void 
ElasticBeam2d::zeroLoad()
{
  Q.Zero();
  q0.zero();
  p0.zero();
  return;
}

int 
ElasticBeam2d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  double L = theCoordTransf->getInitialLength();

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
    double M1 = (0.5*z1*ba2) + (wybma*ba3/(3.0*ba)) - (z1*ba3*2.0/(3.0*L)) - (wybma*ba4/(2.0*L*ba)) + (z1*ba4/(4.0*L2)) + (wybma*ba5/(5.0*L2*ba));
    double M2 = (-1.0*z1*ba3/(3.0*L)) - (wybma*ba4/(4.0*L*ba)) + (z1*ba4/(4.0*L2)) + (wybma*ba5/(5.0*L2*ba));
    double P = waa*ba + 0.5*(wab-waa)*ba;
    double PJ = (1.0/L)*(waa*ba*(a+0.5*ba)+0.5*(wab-waa)*ba*(a+(2.0/3.0)*ba));

    // Reactions in basic system
    p0[0] -= P;
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
    double M1 = a*E*I*(-dT1+(4.0/3.0)*dT); //Fixed End Moment end 1
    double M2 = a*E*I*(dT1+(5.0/3.0)*dT); //Fixed End Moment end 2
    double F = alpha*(((Ttop2+Ttop1)/2+(Tbot2+Tbot1)/2)/2)*E*A; // Fixed End Axial Force
    double M1M2divL =(M1+M2)/L; // Fixed End Shear
    
    // Reactions in basic system
    p0[0] += 0;
    p0[1] += M1M2divL;
    p0[2] -= M1M2divL;

    // Fixed end forces in basic system
    q0[0] -= F;
    q0[1] += M1;
    q0[2] += M2;
  }

  else {
    opserr << "ElasticBeam2d::addLoad()  -- load type unknown for element with tag: " << this->getTag() << endln;
    return -1;
  }

  return 0;
}

int
ElasticBeam2d::addInertiaLoadToUnbalance(const Vector &accel)
{
  if (rho == 0.0)
    return 0;

  // get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
        
  if (3 != Raccel1.Size() || 3 != Raccel2.Size()) {
    opserr << "ElasticBeam2d::addInertiaLoadToUnbalance matrix and vector sizes are incompatable\n";
    return -1;
  }
    
  // want to add ( - fact * M R * accel ) to unbalance
  if (cMass == 0)  {
    // take advantage of lumped mass matrix
    double L = theCoordTransf->getInitialLength();
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
ElasticBeam2d::getResistingForceIncInertia()
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
  
  if (cMass == 0)  {
    // take advantage of lumped mass matrix
    double L = theCoordTransf->getInitialLength();
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
ElasticBeam2d::getResistingForce()
{
  theCoordTransf->update();
  
  const Vector &v = theCoordTransf->getBasicTrialDisp();
  q  = kb*v;
  q += q0;
#if 0
  double L = theCoordTransf->getInitialLength();

  double EoverL   = E/L;
  double EAoverL  = A*EoverL;                        // EA/L
  
  // determine q = kv + q0
  q[0] = EAoverL*v(0);
  if (release == 0) {
    double EIoverL2 = 2.0*I*EoverL;                // 2EI/L
    double EIoverL4 = 2.0*EIoverL2;                // 4EI/L
    q[1] = EIoverL4*v(1) + EIoverL2*v(2);
    q[2] = EIoverL2*v(1) + EIoverL4*v(2);
  }
  if (release == 1) {
    q[1] = 0.0;
    q[2] = 3.0*I*EoverL*v(2);
  }
  if (release == 2) {
    q[1] = 3.0*I*EoverL*v(1);
    q[2] = 0.0;
  }
  if (release == 3) {
    q[1] = 0.0;
    q[2] = 0.0;
  }
  
  q[0] += q0[0];
  q[1] += q0[1];
  q[2] += q0[2];
#endif 

  // Vector for reactions in basic system
  Vector p0Vec(p0);
  
  P = theCoordTransf->getGlobalResistingForce(q, p0Vec);
  if (rho != 0)
    P.addVector(1.0, Q, -1.0);

  return P;
}

int
ElasticBeam2d::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;

    static Vector data(17);
    
    data(0) = A;
    data(1) = E; 
    data(2) = I; 
    data(3) = rho;
    data(4) = cMass;
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
      opserr << "ElasticBeam2d::sendSelf -- could not send data Vector\n";
      return res;
    }

    // Ask the CoordTransf to send itself
    res += theCoordTransf->sendSelf(cTag, theChannel);
    if (res < 0) {
      opserr << "ElasticBeam2d::sendSelf -- could not send CoordTransf\n";
      return res;
    }
    
    return res;
}

int
ElasticBeam2d::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res = 0;
        
    static Vector data(17);

    res += theChannel.recvVector(this->getDbTag(), cTag, data);
    if (res < 0) {
      opserr << "ElasticBeam2d::recvSelf -- could not receive data Vector\n";
      return res;
    }

    A = data(0);
    E = data(1);
    I = data(2); 
    alpha = data(10);
    depth = data(11);

    alphaM = data(12);
    betaK  = data(13);
    betaK0 = data(14);
    betaKc = data(15);
    release = (int)data(16);
    
    rho = data(3);
    cMass = (int)data(4);
    this->setTag((int)data(5));
    connectedExternalNodes(0) = (int)data(6);
    connectedExternalNodes(1) = (int)data(7);

    // Check if the CoordTransf is null; if so, get a new one
    int crdTag = (int)data(8);
    if (theCoordTransf == 0) {
      theCoordTransf = theBroker.getNewCrdTransf(crdTag);
      if (theCoordTransf == 0) {
        opserr << "ElasticBeam2d::recvSelf -- could not get a CrdTransf2d\n";
        exit(-1);
      }
    }
    
    // Check that the CoordTransf is of the right type; if not, delete
    // the current one and get a new one of the right type
    if (theCoordTransf->getClassTag() != crdTag) {
      delete theCoordTransf;
      theCoordTransf = theBroker.getNewCrdTransf(crdTag);
      if (theCoordTransf == 0) {
        opserr << "ElasticBeam2d::recvSelf -- could not get a CrdTransf2d\n";
        exit(-1);
      }
    }
        
    // Now, receive the CoordTransf
    theCoordTransf->setDbTag((int)data(9));
    res += theCoordTransf->recvSelf(cTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "ElasticBeam2d::recvSelf -- could not receive CoordTransf\n";
      return res;
    }
    
    return res;
}

void
ElasticBeam2d::Print(OPS_Stream &s, int flag)
{
  // to update forces!
  this->getResistingForce();

  if (flag == -1) {
    int eleTag = this->getTag();
    s << "EL_BEAM\t" << eleTag << "\t";
    s << 0 << "\t" << 0 << "\t" << connectedExternalNodes(0) << "\t" << connectedExternalNodes(1) ;
    s << "0\t0.0000000\n";
  }

  if (flag == OPS_PRINT_CURRENTSTATE) {
    this->getResistingForce();
    s << "\nElasticBeam2d: " << this->getTag() << endln;
    s << "\tConnected Nodes: " << connectedExternalNodes ;
    s << "\tCoordTransf: " << theCoordTransf->getTag() << endln;
    s << "\tmass density:  " << rho << ", cMass: " << cMass << endln;
    s << "\trelease code:  " << release << endln;
    double P  = q[0];
    double M1 = q[1];
    double M2 = q[2];
    double L = theCoordTransf->getInitialLength();
    double V = (M1+M2)/L;
    s << "\tEnd 1 Forces (P V M): " << -P+p0[0]
      << " " << V+p0[1] << " " << M1 << endln;
    s << "\tEnd 2 Forces (P V M): " << P
      << " " << -V+p0[2] << " " << M2 << endln;
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_ELEM_INDENT << "{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"ElasticBeam2d\", ";
    s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1) << "], ";
        s << "\"E\": " << E << ", ";
        s << "\"A\": "<< A << ", ";
    s << "\"Iz\": "<< I << ", ";
    s << "\"massperlength\": "<< rho << ", ";
    s << "\"mass_flag\": "<< cMass << ", ";
    s << "\"release\": "<< release << ", ";
    s << "\"crdTransformation\": \"" << theCoordTransf->getTag() << "\"}";
  }
}


Response*
ElasticBeam2d::setResponse(const char **argv, int argc, OPS_Stream &output)
{

  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType","ElasticBeam2d");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);

    // global forces
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0 ||
      strcmp(argv[0],"globalForce") == 0 || strcmp(argv[0],"globalForces") == 0) {

    output.tag("ResponseType","Px_1");
    output.tag("ResponseType","Py_1");
    output.tag("ResponseType","Mz_1");
    output.tag("ResponseType","Px_2");
    output.tag("ResponseType","Py_2");
    output.tag("ResponseType","Mz_2");

    theResponse =  new ElementResponse(this, 2, P);
  
  // local forces
  }    else if (strcmp(argv[0],"localForce") == 0 || strcmp(argv[0],"localForces") == 0) {

    output.tag("ResponseType","N_1");
    output.tag("ResponseType","V_1");
    output.tag("ResponseType","M_1");
    output.tag("ResponseType","N_2");
    output.tag("ResponseType","V_2");
    output.tag("ResponseType","M_2");
    
    theResponse = new ElementResponse(this, 3, P);

  // basic forces
  }    else if (strcmp(argv[0],"basicForce") == 0 || strcmp(argv[0],"basicForces") == 0) {

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
ElasticBeam2d::getResponse (int responseID, Information &eleInfo)
{
  double N, M1, M2, V;
  double L = theCoordTransf->getInitialLength();
  this->getResistingForce();

  switch (responseID) {
  case 1: // stiffness
    return eleInfo.setMatrix(this->getTangentStiff());
    
  case 2: // global forces
    return eleInfo.setVector(this->getResistingForce());
    
  case 3: // local forces
    // Axial
    N = q[0];
    P(3) =  N;
    P(0) = -N+p0[0];
    // Moment
    M1 = q[1];
    M2 = q[2];
    P(2) = M1;
    P(5) = M2;
    // Shear
    V = (M1+M2)/L;
    P(1) =  V+p0[1];
    P(4) = -V+p0[2];
    return eleInfo.setVector(P);
    
  case 4: // basic forces
    return eleInfo.setVector(q);

  case 5:
    return eleInfo.setVector(theCoordTransf->getBasicTrialDisp());

  default:
    return -1;
  }
}

int
ElasticBeam2d::setParameter(const char **argv, int argc, Parameter &param)
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
    param.setValue(I);
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
ElasticBeam2d::updateParameter (int parameterID, Information &info)
{
  switch (parameterID) {
  case -1:
    return -1;
  case 1:
    E = info.theDouble;
    return 0;
  case 2:
    A = info.theDouble;
    return 0;
  case 3:
    I = info.theDouble;
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
}

