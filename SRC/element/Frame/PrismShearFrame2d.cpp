//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Prismatic Timoshenko beam
//
// Adapted from ElasticTimoshenkoBeamColumn
//
// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 03/13
// Revision: B
//
#include <PrismShearFrame2d.h>

#include <Domain.h>
#include <Node.h>
#include <CrdTransf.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Renderer.h>
#include <Information.h>
#include <Parameter.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>

#include <math.h>

// initialize the class wide variables
Matrix PrismShearFrame2d::theMatrix(6, 6);
Vector PrismShearFrame2d::theVector(6);


PrismShearFrame2d::PrismShearFrame2d(int tag, int Nd1, int Nd2, double e, double g, double a,
                                     double iz, double avy, CrdTransf& coordTransf, double r,
                                     int cm)
 : Element(tag, ELE_TAG_PrismShearFrame2d),
   connectedExternalNodes(2),
   theCoordTransf(nullptr),
   E(e),
   G(g),
   A(a),
   Iz(iz),
   Avy(avy),
   rho(r),
   cMass(cm),
   geom_flag(0),
   phi(0.0),
   L(0.0),
   ul(6),
   ql(6),
   kl(6, 6),
   klgeo(6, 6),
   Tgl(6, 6),
   Ki(6, 6),
   M(6, 6),
   theLoad(6)
{

  connectedExternalNodes(0) = Nd1;
  connectedExternalNodes(1) = Nd2;

  // set node pointers to NULL
  for (int i = 0; i < 2; i++)
    theNodes[i] = nullptr;

  // get a copy of the coordinate transformation
  theCoordTransf = coordTransf.getCopy2d();
  if (!theCoordTransf) {
    opserr << "PrismShearFrame2d::PrismShearFrame2d() - "
           << "failed to get copy of coordinate transformation.\n";
    exit(-1);
  }

  // get coordinate transformation type and save flag
  if (strncmp(theCoordTransf->getClassType(), "Linear", 6) == 0) {
    geom_flag = 0;
  } else if (strncmp(theCoordTransf->getClassType(), "PDelta", 6) == 0) {
    geom_flag = 1;
  } else if (strncmp(theCoordTransf->getClassType(), "Corot", 5) == 0) {
    geom_flag = 1;
    opserr << "\nWARNING PrismShearFrame2d::PrismShearFrame2d()"
           << " - Element: " << this->getTag() << endln
           << "Unsupported Corotational transformation assigned.\n"
           << "Using PDelta transformation instead.\n";
  }

  // zero fixed end forces vector
  ql0.zero();
}


PrismShearFrame2d::PrismShearFrame2d()
 : Element(0, ELE_TAG_PrismShearFrame2d),
   connectedExternalNodes(2),
   theCoordTransf(0),
   E(0.0),
   G(0.0),
   A(0.0),
   Iz(0.0),
   Avy(0.0),
   rho(0.0),
   cMass(0),
   geom_flag(0),
   phi(0.0),
   L(0.0),
   ul(6),
   kl(6, 6),
   klgeo(6, 6),
   Tgl(6, 6),
   Ki(6, 6),
   M(6, 6),
   theLoad(6)
{

  // set node pointers to NULL
  for (int i = 0; i < 2; i++)
    theNodes[i] = 0;

  // zero fixed end forces vector
  ql0.Zero();
}


PrismShearFrame2d::~PrismShearFrame2d()
{
  if (theCoordTransf)
    delete theCoordTransf;
}


int
PrismShearFrame2d::getNumExternalNodes() const
{
  return 2;
}


const ID&
PrismShearFrame2d::getExternalNodes()
{
  return connectedExternalNodes;
}


Node**
PrismShearFrame2d::getNodePtrs()
{
  return theNodes;
}


int
PrismShearFrame2d::getNumDOF()
{
  return 6;
}


void
PrismShearFrame2d::setDomain(Domain* theDomain)
{
  // check Domain is not null - invoked when object removed from a domain
  if (!theDomain) {
    theNodes[0] = 0;
    theNodes[1] = 0;

    return;
  }

  // first set the node pointers
  theNodes[0] = theDomain->getNode(connectedExternalNodes(0));
  theNodes[1] = theDomain->getNode(connectedExternalNodes(1));

  // if can't find both - send a warning message
  if (!theNodes[0] || !theNodes[1]) {
    if (!theNodes[0]) {
      opserr << "PrismShearFrame2d::setDomain() - Nd1: " << connectedExternalNodes(0)
             << " does not exist in the model for";
    } else {
      opserr << "PrismShearFrame2d::setDomain() - Nd2: " << connectedExternalNodes(1)
             << " does not exist in the model for";
    }
    opserr << " element: " << this->getTag() << ".\n";

    return;
  }

  // now determine the number of dof and the dimension
  int dofNd1 = theNodes[0]->getNumberDOF();
  int dofNd2 = theNodes[1]->getNumberDOF();

  // if differing dof at the ends - print a warning message
  if (dofNd1 != 3) {
    opserr << "PrismShearFrame2d::setDomain() - node 1: " << connectedExternalNodes(0)
           << " has incorrect number of DOF (not 3).\n";
    return;
  }
  if (dofNd2 != 3) {
    opserr << "PrismShearFrame2d::setDomain() - node 2: " << connectedExternalNodes(1)
           << " has incorrect number of DOF (not 3).\n";
    return;
  }

  // call the base class method
  if (theDomain != nullptr)
    this->Element::link(*theDomain);

  // initialize the coordinate transformation
  if (theCoordTransf->initialize(theNodes[0], theNodes[1]) != 0) {
    opserr << "PrismShearFrame2d::setDomain() - "
           << "error initializing coordinate transformation.\n";
    return;
  }

  // set up the transformation matrix for orientation
  this->setUp();
}


int
PrismShearFrame2d::commitState()
{
  int errCode = 0;

  // commit the base class
  errCode += this->Element::commitState();

  // no need to commit coordinate transformation
  // since it is only used to get type of transf
  // errCode += theCoordTransf->commitState();

  return errCode;
}


int
PrismShearFrame2d::revertToLastCommit()
{
  return 0;
}


int
PrismShearFrame2d::revertToStart()
{
  return 0;
}


int
PrismShearFrame2d::update()
{
  return 0;
}


const Matrix&
PrismShearFrame2d::getTangentStiff()
{
  // zero the matrix
  theMatrix.Zero();

  if (geom_flag == 0) {
    // transform from local to global system
    theMatrix.addMatrixTripleProduct(0.0, Tgl, kl, 1.0);

  } else {
    // initialize local stiffness matrix
    static Matrix klTot(6, 6);
    klTot.addMatrix(0.0, kl, 1.0);

    // get global trial displacements
    const Vector& ui = theNodes[0]->getTrialDisp();
    const Vector& uj = theNodes[1]->getTrialDisp();
    static Vector ug(6);
    for (int i = 0; i < 3; i++) {
      ug(i)     = ui(i);
      ug(i + 3) = uj(i);
    }

    // transform response from the global to the local system
    ul.addMatrixVector(0.0, Tgl, ug, 1.0);

    // get the resisting forces in local system
    ql.addMatrixVector(0.0, kl, ul, 1.0);

    // add geometric stiffness to local stiffness
    if (ql(3) != 0.0)
      klTot.addMatrix(1.0, klgeo, ql(3));

    // transform from local to global system
    theMatrix.addMatrixTripleProduct(0.0, Tgl, klTot, 1.0);
  }

  return theMatrix;
}


const Matrix&
PrismShearFrame2d::getInitialStiff()
{
  return Ki;
}


const Matrix&
PrismShearFrame2d::getMass()
{
  return M;
}


void
PrismShearFrame2d::zeroLoad()
{
  theLoad.Zero();
  ql0.Zero();
  return;
}


int
PrismShearFrame2d::addLoad(ElementalLoad* load, double loadFactor)
{
  int type;
  const Vector& data = load->getData(type, loadFactor);

  if (type == LOAD_TAG_Beam2dUniformLoad) {
    double wt = data(0) * loadFactor; // Transverse (+ve upward)
    double wa = data(1) * loadFactor; // Axial (+ve from node I to J)

    double V = 0.5 * wt * L;
    double M = V * L / 6.0; // wt*L*L/12
    double P = 0.5 * wa * L;

    // fixed end forces in local system
    ql0(0) -= P;
    ql0(1) -= V;
    ql0(2) -= M;
    ql0(3) -= P;
    ql0(4) -= V;
    ql0(5) += M;
  }

  else {
    opserr << "PrismShearFrame2d::addLoad() - "
           << "load type unknown for element: " << this->getTag() << ".\n";
    return -1;
  }

  return 0;
}


int
PrismShearFrame2d::addInertiaLoadToUnbalance(const Vector& accel)
{
  // check for quick return
  if (rho == 0.0)
    return 0;

  // assemble Raccel vector
  const Vector& Raccel1 = theNodes[0]->getRV(accel);
  const Vector& Raccel2 = theNodes[1]->getRV(accel);
  static Vector Raccel(6);
  for (int i = 0; i < 3; i++) {
    Raccel(i)     = Raccel1(i);
    Raccel(i + 3) = Raccel2(i);
  }

  // want to add ( - fact * M R * accel ) to unbalance
  theLoad.addMatrixVector(1.0, M, Raccel, -1.0);

  return 0;
}


const Vector&
PrismShearFrame2d::getResistingForce()
{
  // zero the residual
  theVector.Zero();

  // get global trial displacements
  const Vector& ui = theNodes[0]->getTrialDisp();
  const Vector& uj = theNodes[1]->getTrialDisp();
  static Vector ug(6);
  for (int i = 0; i < 3; i++) {
    ug(i)     = ui(i);
    ug(i + 3) = uj(i);
  }

  // transform response from the global to the local system
  ul.addMatrixVector(0.0, Tgl, ug, 1.0);

  // get the resisting forces in local system
  ql.addMatrixVector(0.0, kl, ul, 1.0);

  // add P-Delta moments to local forces
  if ((ql(3) != 0.0) && (geom_flag == 1))
    ql.addMatrixVector(1.0, klgeo, ul, ql(3));

  // add effects of element loads, ql = ql(ul) + ql0
  ql.addVector(1.0, ql0, 1.0);

  // determine resisting forces in global system
  theVector.addMatrixTransposeVector(0.0, Tgl, ql, 1.0);

  return theVector;
}


const Vector&
PrismShearFrame2d::getResistingForceIncInertia()
{
  // first get the resisting forces
  theVector = this->getResistingForce();

  // subtract external load
  theVector.addVector(1.0, theLoad, -1.0);

  // add the damping forces from rayleigh damping
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    theVector.addVector(1.0, this->getRayleighDampingForces(), 1.0);

  // check for quick return
  if (rho == 0.0)
    return theVector;

  // add inertia forces from element mass
  const Vector& accel1 = theNodes[0]->getTrialAccel();
  const Vector& accel2 = theNodes[1]->getTrialAccel();
  static Vector accel(6);
  for (int i = 0; i < 3; i++) {
    accel(i)     = accel1(i);
    accel(i + 3) = accel2(i);
  }
  theVector.addMatrixVector(1.0, M, accel, 1.0);

  return theVector;
}


int
PrismShearFrame2d::sendSelf(int commitTag, Channel& sChannel)
{
  int res = 0;

  static Vector data(16);
  data(0)  = this->getTag();
  data(1)  = connectedExternalNodes(0);
  data(2)  = connectedExternalNodes(1);
  data(3)  = E;
  data(4)  = G;
  data(5)  = A;
  data(6)  = Iz;
  data(7)  = Avy;
  data(8)  = rho;
  data(9)  = cMass;
  data(10) = alphaM;
  data(11) = betaK;
  data(12) = betaK0;
  data(13) = betaKc;
  data(14) = theCoordTransf->getClassTag();

  int dbTag = theCoordTransf->getDbTag();
  if (dbTag == 0) {
    dbTag = sChannel.getDbTag();
    if (dbTag != 0)
      theCoordTransf->setDbTag(dbTag);
  }
  data(15) = dbTag;

  // send the data vector
  res += sChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "PrismShearFrame2d::sendSelf() - could not send data Vector.\n";
    return res;
  }

  // ask the CoordTransf to send itself
  res += theCoordTransf->sendSelf(commitTag, sChannel);
  if (res < 0) {
    opserr << "PrismShearFrame2d::sendSelf() - could not send CoordTransf.\n";
    return res;
  }

  return res;
}


int
PrismShearFrame2d::recvSelf(int commitTag, Channel& rChannel, FEM_ObjectBroker& theBroker)
{
  int res = 0;

  static Vector data(16);
  res += rChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "PrismShearFrame2d::recvSelf() - could not receive data Vector.\n";
    return res;
  }

  this->setTag((int)data(0));
  connectedExternalNodes(0) = (int)data(1);
  connectedExternalNodes(1) = (int)data(2);
  E                         = data(3);
  G                         = data(4);
  A                         = data(5);
  Iz                        = data(6);
  Avy                       = data(7);
  rho                       = data(8);
  cMass                     = (int)data(9);
  alphaM                    = data(10);
  betaK                     = data(11);
  betaK0                    = data(12);
  betaKc                    = data(13);

  // check if the CoordTransf is null; if so, get a new one
  int crdTag = (int)data(14);
  if (theCoordTransf == 0) {
    theCoordTransf = theBroker.getNewCrdTransf(crdTag);
    if (theCoordTransf == 0) {
      opserr << "PrismShearFrame2d::recvSelf() - could not get a CrdTransf2d.\n";
      return -1;
    }
  }

  // check that the CoordTransf is of the right type; if not, delete
  // the current one and get a new one of the right type
  if (theCoordTransf->getClassTag() != crdTag) {
    delete theCoordTransf;
    theCoordTransf = theBroker.getNewCrdTransf(crdTag);
    if (theCoordTransf == 0) {
      opserr << "PrismShearFrame2d::recvSelf() - could not get a CrdTransf2d.\n";
      return -1;
    }
  }

  // receive the CoordTransf
  theCoordTransf->setDbTag((int)data(15));
  res += theCoordTransf->recvSelf(commitTag, rChannel, theBroker);
  if (res < 0) {
    opserr << "PrismShearFrame2d::recvSelf() - could not receive CoordTransf.\n";
    return res;
  }

  // get coordinate transformation type and save flag
  if (strncmp(theCoordTransf->getClassType(), "Linear", 6) == 0) {
    geom_flag = 0;
  } else if (strncmp(theCoordTransf->getClassType(), "PDelta", 6) == 0) {
    geom_flag = 1;
  } else if (strncmp(theCoordTransf->getClassType(), "Corot", 5) == 0) {
    geom_flag = 1;
    opserr << "\nWARNING PrismShearFrame2d::recvSelf()"
           << " - Element: " << this->getTag() << endln
           << "Unsupported Corotational transformation assigned.\n"
           << "Using PDelta transformation instead.\n";
  }

  // revert the CoordTransf to its last committed state
  theCoordTransf->revertToLastCommit();

  return res;
}


void
PrismShearFrame2d::Print(OPS_Stream& s, int flag)
{
  if (flag == OPS_PRINT_CURRENTSTATE) {
    // print everything
    s << "Element: " << this->getTag();
    s << "  type: PrismShearFrame2d";
    s << "  iNode: " << connectedExternalNodes(0);
    s << "  jNode: " << connectedExternalNodes(1) << endln;
    s << "  E: " << E << "  G: " << G << endln;
    s << "  A: " << A << "  Iz: " << Iz << "  Avy: " << Avy << endln;
    s << "  coordTransf: " << theCoordTransf->getClassType() << endln;
    s << "  rho: " << rho << "  cMass: " << cMass << endln;
    // determine resisting forces in global system
    s << "  resisting force: " << this->getResistingForce() << endln;
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"PrismShearFrame2d\", ";
    s << "\"nodes\": [" 
      << connectedExternalNodes(0) << ", " 
      << connectedExternalNodes(1) << "], ";
    s << "\"E\": " << E << ", ";
    s << "\"G\": " << G << ", ";
    s << "\"A\": " << A << ", ";
    s << "\"Avy\": " << Avy << ", ";
    s << "\"Iz\": " << Iz << ", ";
    s << "\"massperlength\": " << rho << ", ";
    s << "\"crdTransformation\": \"" << theCoordTransf->getTag() << "\"}";
  }
}


Response*
PrismShearFrame2d::setResponse(const char** argv, int argc, OPS_Stream& output)
{
  Response* theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType", "PrismShearFrame2d");
  output.attr("eleTag", this->getTag());
  output.attr("node1", connectedExternalNodes[0]);
  output.attr("node2", connectedExternalNodes[1]);

  // global forces
  if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0 ||
      strcmp(argv[0], "globalForce") == 0 || strcmp(argv[0], "globalForces") == 0) {
    output.tag("ResponseType", "Px_1");
    output.tag("ResponseType", "Py_1");
    output.tag("ResponseType", "Mz_1");
    output.tag("ResponseType", "Px_2");
    output.tag("ResponseType", "Py_2");
    output.tag("ResponseType", "Mz_2");

    theResponse = new ElementResponse(this, 1, theVector);
  }
  // local forces
  else if (strcmp(argv[0], "localForce") == 0 || strcmp(argv[0], "localForces") == 0) {
    output.tag("ResponseType", "N_1");
    output.tag("ResponseType", "V_1");
    output.tag("ResponseType", "M_1");
    output.tag("ResponseType", "N_2");
    output.tag("ResponseType", "V_2");
    output.tag("ResponseType", "M_2");

    theResponse = new ElementResponse(this, 2, theVector);
  }

  output.endTag(); // ElementOutput

  return theResponse;
}


int
PrismShearFrame2d::getResponse(int responseID, Information& info)
{
  switch (responseID) {
  case 1: // global forces
    return info.setVector(this->getResistingForce());

  case 2: // local forces
    theVector.Zero();
    // determine resisting forces in local system
    theVector = ql;

    return info.setVector(theVector);

  default: return -1;
  }
}

int
PrismShearFrame2d::setParameter(const char** argv, int argc, Parameter& param)
{
  if (argc < 1)
    return -1;

  // E of the beam
  if (strcmp(argv[0], "E") == 0)
    return param.addObject(1, this);

  // G of the beam
  if (strcmp(argv[0], "G") == 0)
    return param.addObject(2, this);

  // A of the beam
  if (strcmp(argv[0], "A") == 0)
    return param.addObject(3, this);

  // Iz of the beam
  if (strcmp(argv[0], "Iz") == 0)
    return param.addObject(4, this);

  // Avy of the beam
  if (strcmp(argv[0], "Avy") == 0)
    return param.addObject(5, this);

  return -1;
}


int
PrismShearFrame2d::updateParameter(int parameterID, Information& info)
{
  switch (parameterID) {
  case -1: return -1;
  case 1:  E = info.theDouble; return 0;
  case 2:  G = info.theDouble; return 0;
  case 3:  A = info.theDouble; return 0;
  case 4:  Iz = info.theDouble; return 0;
  case 5:  Avy = info.theDouble; return 0;
  default: return -1;
  }
}


void
PrismShearFrame2d::setUp()
{
  // element projection
  static Vector dx(2);

  const Vector& ndICoords = theNodes[0]->getCrds();
  const Vector& ndJCoords = theNodes[1]->getCrds();

  dx = ndJCoords - ndICoords;

  // determine the element length
  L = theCoordTransf->getInitialLength();

  if (L == 0.0) {
    opserr << "PrismShearFrame2d::setUp()  - "
           << "element: " << this->getTag() << " has zero length.\n";
    return;
  }

  // create transformation matrix from global to local system
  Tgl.Zero();
  Tgl(0, 0) = Tgl(1, 1) = Tgl(3, 3) = Tgl(4, 4) = dx(0) / L;
  Tgl(0, 1) = Tgl(3, 4) = dx(1) / L;
  Tgl(1, 0) = Tgl(4, 3) = -dx(1) / L;
  Tgl(2, 2) = Tgl(5, 5) = 1.0;

  // determine ratio of bending to shear stiffness
  phi = 12.0 * E * Iz / (L * L * G * Avy);

  // compute initial stiffness matrix in local system
  kl.Zero();
  kl(0, 0) = kl(3, 3) = E * A / L;
  kl(0, 3) = kl(3, 0) = -kl(0, 0);
  double a1z          = E * Iz / (L * L * L * (1.0 + phi));
  kl(1, 1) = kl(4, 4) = a1z * 12.0;
  kl(1, 4) = kl(4, 1) = -kl(1, 1);
  kl(2, 2) = kl(5, 5) = a1z * L * L * (4.0 + phi);
  kl(2, 5) = kl(5, 2) = a1z * L * L * (2.0 - phi);
  kl(1, 2) = kl(2, 1) = kl(1, 5) = kl(5, 1) = a1z * L * 6.0;
  kl(2, 4) = kl(4, 2) = kl(4, 5) = kl(5, 4) = -kl(1, 2);

  // compute geometric stiffness matrix in local system
  klgeo.Zero();
  if (geom_flag == 1) {
    double b1z  = 1.0 / (30.0 * L * pow(1.0 + phi, 2));
    klgeo(1, 1) = klgeo(4, 4) = b1z * (30.0 * phi * phi + 60.0 * phi + 36.0);
    klgeo(1, 4) = klgeo(4, 1) = -klgeo(1, 1);
    klgeo(2, 2) = klgeo(5, 5) = b1z * L * L * (2.5 * phi * phi + 5.0 * phi + 4.0);
    klgeo(2, 5) = klgeo(5, 2) = -b1z * L * L * (2.5 * phi * phi + 5.0 * phi + 1.0);
    klgeo(1, 2) = klgeo(2, 1) = klgeo(1, 5) = klgeo(5, 1) = b1z * L * 3.0;
    klgeo(2, 4) = klgeo(4, 2) = klgeo(4, 5) = klgeo(5, 4) = -klgeo(1, 2);
  }

  // compute initial stiffness matrix in global system
  Ki.addMatrixTripleProduct(0.0, Tgl, kl, 1.0);

  // compute mass matrix in global system
  M.Zero();
  if (rho > 0.0) {
    if (cMass == 0) {
      // lumped mass matrix
      double m = 0.5 * rho * L;
      for (int i = 0; i < 2; i++) {
        M(i, i)         = m;
        M(i + 3, i + 3) = m;
      }
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
      M.addMatrixTripleProduct(0.0, Tgl, ml, 1.0);
    }
  }
}
