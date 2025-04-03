//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// This element is adapted from TimoshenkoBeamColumn3d
//
// Written: MHS
// Created: Feb 2001
//
#include <CubicFrame3d.h>
#include <Node.h>
#include <FrameSection.h>
#include <FrameTransform.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <CompositeResponse.h>
#include <ElementalLoad.h>
#include <BeamIntegration.h>
#include <Parameter.h>
#include <math.h>

Matrix CubicFrame3d::K(12, 12);
Vector CubicFrame3d::P(12);
using namespace OpenSees;

#define ELE_TAG_CubicFrame3d 0

CubicFrame3d::CubicFrame3d(int tag, 
                           std::array<int, 2>& nodes, 
                           std::vector<FrameSection*>& sections,
                           BeamIntegration& bi,
                           FrameTransform3d& coordTransf, 
                           double r)
 : Element(tag, ELE_TAG_CubicFrame3d),
   numSections(sections.size()),
   theSections(nullptr),
   theCoordTransf(nullptr),
   beamInt(nullptr),
   Q(12),
   q(6),
   density(r),
   parameterID(0),
   connectedExternalNodes(2)
{
  // Allocate arrays of pointers to SectionForceDeformations
  theSections = new FrameSection*[numSections];

  for (int i = 0; i < numSections; i++) {
    // Get copies of the material for each integration point
    theSections[i] = sections[i]->getFrameCopy(scheme);
  }

  beamInt = bi.getCopy();


  theCoordTransf = coordTransf.getCopy();

  // Set connected external node IDs
  connectedExternalNodes(0) = nodes[0];
  connectedExternalNodes(1) = nodes[1];


  theNodes[0] = nullptr;
  theNodes[1] = nullptr;

  q0.zero();
  p0.zero();
}

CubicFrame3d::CubicFrame3d()
 : Element(0, ELE_TAG_CubicFrame3d),
   numSections(0),
   theSections(nullptr),
   theCoordTransf(nullptr),
   beamInt(nullptr),
   connectedExternalNodes(2),
   Q(12),
   q(6),
   density(0.0),
   parameterID(0)
{
  q0.zero();
  p0.zero();

  theNodes[0] = nullptr;
  theNodes[1] = nullptr;
}

CubicFrame3d::~CubicFrame3d()
{
  for (int i = 0; i < numSections; i++) {
    if (theSections[i])
      delete theSections[i];
  }

  // Delete the array of pointers to SectionForceDeformation pointer arrays
  if (theSections)
    delete[] theSections;

  if (theCoordTransf)
    delete theCoordTransf;

  if (beamInt != 0)
    delete beamInt;
}

int
CubicFrame3d::getNumExternalNodes() const
{
  return 2;
}

const ID&
CubicFrame3d::getExternalNodes()
{
  return connectedExternalNodes;
}

Node**
CubicFrame3d::getNodePtrs()
{
  return theNodes;
}

int
CubicFrame3d::getNumDOF()
{
  return 12;
}

void
CubicFrame3d::setDomain(Domain* theDomain)
{
  // Check Domain is not null. This happens when element is removed from a domain.
  // In this case just set null pointers to null and return.
  if (theDomain == nullptr) {
    for (int i=0; i<NEN; i++)
      theNodes[i] = nullptr;
    return;
  }

  for (int i=0; i<NEN; i++) {
    // Retrieve the node from the domain using its tag.
    // If no node is found, then return
    theNodes[i] = theDomain->getNode(connectedExternalNodes(i));
    if (theNodes[i] == nullptr)
      return;

    // If node is found, ensure node has the proper number of DOFs
    int dofs = theNodes[i]->getNumberDOF();
    if (dofs != NDF) {
      opserr << "WARNING " << this->getClassType() << " element " << this->getTag() 
             << " does not have " << NDF << " DOFs at node " 
             << theNodes[i]->getTag() << "\n";
      return;
    }
  }

  if (theCoordTransf->initialize(theNodes[0], theNodes[NEN-1])) {
    // Add some error check
  }

  double L = theCoordTransf->getInitialLength();

  if (L == 0.0) {
    // Add some error check
  }

  for (int i = 0; i < numSections; i++) {
    const Matrix& ks0 = theSections[i]->getInitialTangent();
    int order         = theSections[i]->getOrder();
    const ID& code    = theSections[i]->getType();

    double EI = 0.0;
    double GA = 0.0;
    for (int k = 0; k < nsr; k++) {
      if (code(k) == SECTION_RESPONSE_MZ)
        EI += ks0(k, k);
      if (code(k) == SECTION_RESPONSE_VY)
        GA += ks0(k, k);
    }
    phizs[i] = 0.0;
    if (GA != 0.0)
      phizs[i] = 12 * EI / (GA * L * L);

    EI = 0.0;
    GA = 0.0;
    for (int k = 0; k < nsr; k++) {
      if (code(k) == SECTION_RESPONSE_MY)
        EI += ks0(k, k);
      if (code(k) == SECTION_RESPONSE_VZ)
        GA += ks0(k, k);
    }
    phiys[i] = 0.0;
    if (GA != 0.0)
      phiys[i] = 12 * EI / (GA * L * L);
  }

  beamInt->getSectionLocations(numSections, L, xi);
  beamInt->getSectionWeights(numSections, L, wt);

  this->DomainComponent::setDomain(theDomain);

  this->update();
}

int
CubicFrame3d::commitState()
{
  int status = 0;

  // call element commitState to do any base class stuff
  if ((status = this->Element::commitState()) != 0) {
    opserr << "CubicFrame3d::commitState () - failed in base class";
  }

  // Loop over the integration points and commit the material states
  for (int i = 0; i < numSections; i++)
    status += theSections[i]->commitState();

  status += theCoordTransf->commitState();

  return status;
}

int
CubicFrame3d::revertToLastCommit()
{
  int status = 0;

  // Loop over the integration points and revert to last committed state
  for (int i = 0; i < numSections; i++)
    status += theSections[i]->revertToLastCommit();

  status += theCoordTransf->revertToLastCommit();

  return status;
}

int
CubicFrame3d::revertToStart()
{
  int status = 0;

  // Loop over the integration points and revert states to start
  for (int i = 0; i < numSections; i++)
    status += theSections[i]->revertToStart();

  status += theCoordTransf->revertToStart();

  return status;
}

int
CubicFrame3d::update()
{
  int err = 0;

  // Update the transformation
  theCoordTransf->update();

  // Get basic deformations
  const Vector& v = theCoordTransf->getBasicTrialDisp();

  double L        = theCoordTransf->getInitialLength();
  double jsx = 1.0 / L;

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    int order      = theSections[i]->getOrder();
    const ID& code = theSections[i]->getType();


    double xi6  = 6.0 * xi[i];
    double phiz = phizs[i];
    double phiy = phiys[i];

    VectorND<nsr> e;
    for (int j = 0; j < nsr; j++) {
      switch (scheme[j]) {
      case SECTION_RESPONSE_P: 
        e(j) = jsx * v(0);
        break;
      case SECTION_RESPONSE_VY:
        e(j) = 0.5 * phiz / (1 + phiz) * v(1) + 0.5 * phiz / (1 + phiz) * v(2);
        break;
      case SECTION_RESPONSE_VZ:
        e(j) = 0.5 * phiy / (1 + phiy) * v(3) + 0.5 * phiy / (1 + phiy) * v(4);
        break;
      case SECTION_RESPONSE_T: 
        e(j) = jsx * v(5); 
        break;
      case SECTION_RESPONSE_MY:
        e(j) = jsx / (1 + phiy) * ((xi6 - 4.0 - phiy) * v(3) + (xi6 - 2.0 + phiy) * v(4));
        break;
      case SECTION_RESPONSE_MZ:
        e(j) = jsx / (1 + phiz) * ((xi6 - 4.0 - phiz) * v(1) + (xi6 - 2.0 + phiz) * v(2));
        break;
      default:
        e(j) = 0.0; 
        break;
      }
    }

    // Set the section deformations
    err += theSections[i]->setTrialState<nsr,scheme>(e);
  }

  return err;
}

const Matrix&
CubicFrame3d::getTangentStiff()
{
  static MatrixND<6,6> kb;
  static Matrix wrapper(kb);


  double L   = theCoordTransf->getInitialLength();
  double jsx = 1.0 / L;


  kb.zero();
  q.Zero();
  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    MatrixND<nsr,6> ka;
    ka.zero();

    double xi6  = 6.0 * xi[i];
    double phiz = phizs[i];
    double phiy = phiys[i];

    // Get the section tangent stiffness and stress resultant
     MatrixND<nsr,nsr> ks = theSections[i]->getTangent<nsr,scheme>(State::Pres);
    const VectorND<nsr> s = theSections[i]->getResultant<nsr,scheme>();
    // Perform numerical integration
    //kb.addMatrixTripleProduct(1.0, *B, ks, wts(i)/L);
    double wti = wt[i] * jsx;
    for (int j = 0; j < nsr; j++) {
      double tmp;
      switch (scheme[j]) {
      case SECTION_RESPONSE_P:
        for (int k = 0; k < nsr; k++)
          ka(k, 0) += ks(k, j) * wti;
        break;
      case SECTION_RESPONSE_MZ:
        for (int k = 0; k < nsr; k++) {
          tmp = ks(k, j) * wti;
          ka(k, 1) += 1.0 / (1 + phiz) * (xi6 - 4.0 - phiz) * tmp;
          ka(k, 2) += 1.0 / (1 + phiz) * (xi6 - 2.0 + phiz) * tmp;
        }
        break;
      case SECTION_RESPONSE_MY:
        for (int k = 0; k < nsr; k++) {
          tmp = ks(k, j) * wti;
          ka(k, 3) += 1.0 / (1 + phiy) * (xi6 - 4.0 - phiy) * tmp;
          ka(k, 4) += 1.0 / (1 + phiy) * (xi6 - 2.0 + phiy) * tmp;
        }
        break;
      case SECTION_RESPONSE_VY:
        for (int k = 0; k < nsr; k++) {
          tmp = ks(k, j) * wti;
          ka(k, 1) += 0.5 * phiz * L / (1 + phiz) * tmp;
          ka(k, 2) += 0.5 * phiz * L / (1 + phiz) * tmp;
        }
        break;
      case SECTION_RESPONSE_VZ:
        for (int k = 0; k < nsr; k++) {
          tmp = ks(k, j) * wti;
          ka(k, 3) += 0.5 * phiy * L / (1 + phiy) * tmp;
          ka(k, 4) += 0.5 * phiy * L / (1 + phiy) * tmp;
        }
        break;
      case SECTION_RESPONSE_T:
        for (int k = 0; k < nsr; k++)
          ka(k, 5) += ks(k, j) * wti;
        break;
      default: break;
      }
    }
    for (int j = 0; j < nsr; j++) {
      double tmp;
      switch (scheme[j]) {
      case SECTION_RESPONSE_P:
        for (int k = 0; k < 6; k++)
          kb(0, k) += ka(j, k);
        break;
      case SECTION_RESPONSE_MZ:
        for (int k = 0; k < 6; k++) {
          tmp = ka(j, k);
          kb(1, k) += 1.0 / (1 + phiz) * (xi6 - 4.0 - phiz) * tmp;
          kb(2, k) += 1.0 / (1 + phiz) * (xi6 - 2.0 + phiz) * tmp;
        }
        break;
      case SECTION_RESPONSE_MY:
        for (int k = 0; k < 6; k++) {
          tmp = ka(j, k);
          kb(3, k) += 1.0 / (1 + phiy) * (xi6 - 4.0 - phiy) * tmp;
          kb(4, k) += 1.0 / (1 + phiy) * (xi6 - 2.0 + phiy) * tmp;
        }
        break;
      case SECTION_RESPONSE_VY:
        for (int k = 0; k < 6; k++) {
          tmp = ka(j, k);
          kb(1, k) += 0.5 * phiz * L / (1 + phiz) * tmp;
          kb(2, k) += 0.5 * phiz * L / (1 + phiz) * tmp;
        }
        break;
      case SECTION_RESPONSE_VZ:
        for (int k = 0; k < 6; k++) {
          tmp = ka(j, k);
          kb(3, k) += 0.5 * phiy * L / (1 + phiy) * tmp;
          kb(4, k) += 0.5 * phiy * L / (1 + phiy) * tmp;
        }
        break;
      case SECTION_RESPONSE_T:
        for (int k = 0; k < 6; k++)
          kb(5, k) += ka(j, k);
        break;
      default: break;
      }
    }

    //q.addMatrixTransposeVector(1.0, *B, s, wts(i));
    for (int j = 0; j < nsr; j++) {
      double si = s[j] * wt[i];
      switch (scheme[j]) {
      case SECTION_RESPONSE_P: 
        q(0) += si; break;
      case SECTION_RESPONSE_MZ:
        q(1) += 1.0 / (1 + phiz) * (xi6 - 4.0 - phiz) * si;
        q(2) += 1.0 / (1 + phiz) * (xi6 - 2.0 + phiz) * si;
        break;
      case SECTION_RESPONSE_MY:
        q(3) += 1.0 / (1 + phiy) * (xi6 - 4.0 - phiy) * si;
        q(4) += 1.0 / (1 + phiy) * (xi6 - 2.0 + phiy) * si;
        break;
      case SECTION_RESPONSE_VY:
        q(1) += 0.5 * phiz * L / (1 + phiz) * si;
        q(2) += 0.5 * phiz * L / (1 + phiz) * si;
        break;
      case SECTION_RESPONSE_VZ:
        q(3) += 0.5 * phiy * L / (1 + phiy) * si;
        q(4) += 0.5 * phiy * L / (1 + phiy) * si;
        break;
      case SECTION_RESPONSE_T: 
        q(5) += si;
        break;
      default:
        break;
      }
    }
  }

  q[0] += q0[0];
  q[1] += q0[1];
  q[2] += q0[2];
  q[3] += q0[3];
  q[4] += q0[4];

  // Transform to global stiffness
  K = theCoordTransf->getGlobalStiffMatrix(wrapper, q);
  return K;
}

void
CubicFrame3d::getBasicStiff(Matrix& kb, int initial)
{
  // Zero for integral
  kb.Zero();

  double L   = theCoordTransf->getInitialLength();
  double jsx = 1.0 / L;

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    MatrixND<nsr,6> ka;
    ka.zero();

    double xi6  = 6.0 * xi[i];
    double phiz = phizs[i];
    double phiy = phiys[i];

    // Get the section tangent stiffness
    const Matrix& ks =
        (initial) ? theSections[i]->getInitialTangent() 
        : theSections[i]->getSectionTangent();

    // Perform numerical integration
    //kb.addMatrixTripleProduct(1.0, *B, ks, wts(i)/L);
    double wti = wt[i] * jsx;
    for (int j = 0; j < nsr; j++) {
      double tmp;
      switch (scheme[j]) {
      case SECTION_RESPONSE_P:
        for (int k = 0; k < nsr; k++)
          ka(k, 0) += ks(k, j) * wti;
        break;
      case SECTION_RESPONSE_MZ:
        for (int k = 0; k < nsr; k++) {
          tmp = ks(k, j) * wti;
          ka(k, 1) += 1.0 / (1 + phiz) * (xi6 - 4.0 - phiz) * tmp;
          ka(k, 2) += 1.0 / (1 + phiz) * (xi6 - 2.0 + phiz) * tmp;
        }
        break;
      case SECTION_RESPONSE_MY:
        for (int k = 0; k < nsr; k++) {
          tmp = ks(k, j) * wti;
          ka(k, 3) += 1.0 / (1 + phiy) * (xi6 - 4.0 - phiy) * tmp;
          ka(k, 4) += 1.0 / (1 + phiy) * (xi6 - 2.0 + phiy) * tmp;
        }
        break;
      case SECTION_RESPONSE_VY:
        for (int k = 0; k < nsr; k++) {
          tmp = ks(k, j) * wti;
          ka(k, 1) += 0.5 * phiz * L / (1 + phiz) * tmp;
          ka(k, 2) += 0.5 * phiz * L / (1 + phiz) * tmp;
        }
        break;
      case SECTION_RESPONSE_VZ:
        for (int k = 0; k < nsr; k++) {
          tmp = ks(k, j) * wti;
          ka(k, 3) += 0.5 * phiy * L / (1 + phiy) * tmp;
          ka(k, 4) += 0.5 * phiy * L / (1 + phiy) * tmp;
        }
        break;
      case SECTION_RESPONSE_T:
        for (int k = 0; k < nsr; k++)
          ka(k, 5) += ks(k, j) * wti;
        break;
      default: break;
      }
    }

    for (int j = 0; j < nsr; j++) {
      switch (scheme[j]) {
      case SECTION_RESPONSE_P:
        for (int k = 0; k < 6; k++)
          kb(0, k) += ka(j, k);
        break;
      case SECTION_RESPONSE_MZ:
        for (int k = 0; k < 6; k++) {
          double tmp = ka(j, k);
          kb(1, k) += 1.0 / (1 + phiz) * (xi6 - 4.0 - phiz) * tmp;
          kb(2, k) += 1.0 / (1 + phiz) * (xi6 - 2.0 + phiz) * tmp;
        }
        break;
      case SECTION_RESPONSE_MY:
        for (int k = 0; k < 6; k++) {
          double tmp = ka(j, k);
          kb(3, k) += 1.0 / (1 + phiy) * (xi6 - 4.0 - phiy) * tmp;
          kb(4, k) += 1.0 / (1 + phiy) * (xi6 - 2.0 + phiy) * tmp;
        }
        break;
      case SECTION_RESPONSE_VY:
        for (int k = 0; k < 6; k++) {
          double tmp = ka(j, k);
          kb(1, k) += 0.5 * phiz * L / (1 + phiz) * tmp;
          kb(2, k) += 0.5 * phiz * L / (1 + phiz) * tmp;
        }
        break;
      case SECTION_RESPONSE_VZ:
        for (int k = 0; k < 6; k++) {
          double tmp = ka(j, k);
          kb(3, k) += 0.5 * phiy * L / (1 + phiy) * tmp;
          kb(4, k) += 0.5 * phiy * L / (1 + phiy) * tmp;
        }
        break;
      case SECTION_RESPONSE_T:
        for (int k = 0; k < 6; k++)
          kb(5, k) += ka(j, k);
        break;
      default: break;
      }
    }
  }
}

const Matrix&
CubicFrame3d::getInitialStiff()
{
  static Matrix kb(6, 6);

  this->getBasicStiff(kb, 1);

  // Transform to global stiffness
  K = theCoordTransf->getInitialGlobalStiffMatrix(kb);

  return K;
}

const Matrix&
CubicFrame3d::getMass()
{
  K.Zero();

  if (density == 0.0)
    return K;

  double L = theCoordTransf->getInitialLength();

  // lumped mass matrix
  double m = 0.5 * density * L;
  K(0, 0) = K(1, 1) = K(2, 2) = K(6, 6) = K(7, 7) = K(8, 8) = m;

  return K;
}

void
CubicFrame3d::zeroLoad()
{
  Q.Zero();

  q0.zero();
  p0.zero();

  return;
}

int
CubicFrame3d::addLoad(ElementalLoad* theLoad, double loadFactor)
{
  int type;
  const Vector& data = theLoad->getData(type, loadFactor);
  double L           = theCoordTransf->getInitialLength();

  if (type == LOAD_TAG_Beam3dUniformLoad) {
    double wy = data(0) * loadFactor; // Transverse
    double wz = data(1) * loadFactor; // Transverse
    double wx = data(2) * loadFactor; // Axial (+ve from node I to J)

    double Vy = 0.5 * wy * L;
    double Mz = Vy * L / 6.0; // wy*L*L/12
    double Vz = 0.5 * wz * L;
    double My = Vz * L / 6.0; // wz*L*L/12
    double P  = wx * L;

    // Reactions in basic system
    p0[0] -= P;
    p0[1] -= Vy;
    p0[2] -= Vy;
    p0[3] -= Vz;
    p0[4] -= Vz;

    // Fixed end forces in basic system
    q0[0] -= 0.5 * P;
    q0[1] -= Mz;
    q0[2] += Mz;
    q0[3] += My;
    q0[4] -= My;

  } else if (type == LOAD_TAG_Beam3dPointLoad) {
    double Py     = data(0) * loadFactor;
    double Pz     = data(1) * loadFactor;
    double N      = data(2) * loadFactor;
    double aOverL = data(3);

    if (aOverL < 0.0 || aOverL > 1.0)
      return 0;

    double a = aOverL * L;
    double b = L - a;

    // Reactions in basic system
    p0[0] -= N;
    double V1, V2;
    V1 = Py * (1.0 - aOverL);
    V2 = Py * aOverL;
    p0[1] -= V1;
    p0[2] -= V2;
    V1 = Pz * (1.0 - aOverL);
    V2 = Pz * aOverL;
    p0[3] -= V1;
    p0[4] -= V2;

    double L2 = 1.0 / (L * L);
    double a2 = a * a;
    double b2 = b * b;

    // Fixed end forces in basic system
    q0[0] -= N * aOverL;
    double M1, M2;
    M1 = -a * b2 * Py * L2;
    M2 = a2 * b * Py * L2;
    q0[1] += M1;
    q0[2] += M2;
    M1 = -a * b2 * Pz * L2;
    M2 = a2 * b * Pz * L2;
    q0[3] -= M1;
    q0[4] -= M2;
  }

#if 0
  else if (type == LOAD_TAG_Beam3dThermalAction) {

     // load not inside fire load pattern
           //static Vector factors(9);
           //factors.Zero();
           //factors += loadFactor;
          // return this->addLoad(theLoad, factors);
      // This code block is added by LJ and copied from DispBeamColumn2d(Modified edition) for 'FireLoadPattern'--08-May-2012--//[END]
      counterTemperature = 1;
      for(int i=0; i<5; i++){
          residThermal[i]=0;
      }
      double ThermalN, ThermalMz, ThermalMy;
      ThermalN = 0.;
      ThermalMz = 0;
      ThermalMy = 0;

      double xi[maxNumSections];
      beamInt->getSectionLocations(numSections, L, xi);
      double wt[maxNumSections];
      beamInt->getSectionWeights(numSections, L, wt);

        // Zero for integration
        //q.Zero();
      Vector* dataMixV;
      dataMixV = new Vector(data.Size());

      *dataMixV=data;

      // Loop over the integration points
      for (int i = 0; i < numSections; i++) {
        // Get section stress resultant
        const Vector &s = theSections[i]->getTemperatureStress(*dataMixV);

        //apply temp along y
        residThermal[0] = -s(0);
        residThermal[1] = -s(1);
        residThermal[2] = s(1);
        residThermal[3] = -s(2);
        residThermal[4] = s(2);
        //apply temp along z
        //residThermal[0] = -s(0);
        //residThermal[1] = -0.;
        //residThermal[2] = -0.;
        //residThermal[3] = -s(1);
        //residThermal[4] = s(1);
        SectionThermalElong[i]=0;
      }
      AverageThermalElong=0;
    }

    //Added by Liming for implementation of NodalThermalAction
    else if (type == LOAD_TAG_NodalThermalAction) {

        // load not inside fire load pattern
        AverageThermalElong=0.0;
        for(int i=0; i<5; i++){
            residThermal[i]=0;
        }

        NodalThermalAction* theNodalThermal0 = theNodes[0]->getNodalThermalActionPtr();
        NodalThermalAction* theNodalThermal1 = theNodes[1]->getNodalThermalActionPtr();
        int type;
        const Vector &data0 = theNodalThermal0->getData(type);
        const Vector &data1 = theNodalThermal1->getData(type);
        Vector* Loc;Vector* NodalT0;Vector* NodalT1;

        if (data0.Size()==9){
           Loc = new Vector(9);
           NodalT0 = new Vector(9);
           NodalT1 = new Vector(9);

           for(int i =0; i<9;i++){
               if(data0(2*i+1)-data1(2*i+1)>1e-8||data0(2*i+1)-data1(2*i+1)<-1e-8){
                       opserr<<"Warning:The NodalThermalAction in dispBeamColumn2dThermalNUT "<<this->getTag()
                            << "incompatible loc input for datapoint "<< i << endln;
               }
               else {
                   (*Loc)(i)=data0(2*i+1);
                   (*NodalT0)(i)=data0(2*i);
                   (*NodalT1)(i)=data1(2*i);
              }
           }
        }

        //for 9 data POINTS
        else {
          Loc = new Vector(10);
          NodalT0 = new Vector(15);
          NodalT1 = new Vector(15);

          for(int i =0; i<5;i++){
           if (data0(2*i+1)-data1(2*i+1)>1e-8||data0(2*i+1)-data1(2*i+1)<-1e-8){
                  opserr <<"Warning:The NodalThermalAction in dispBeamColumn2dThermalNUT "<<this->getTag()
                         << "incompatible loc input for datapoint "<< i << endln;
           } else {
                 //for loc
                (*Loc)(i)=data0(2*i+1);
                (*Loc)(i+5)=data0(3*i+12);
                //for NodalT0
                (*NodalT0)(i)=data0(2*i);
                (*NodalT0)(i+5)=data0(3*i+10);
                (*NodalT0)(i+10)=data0(3*i+11);
                //for NodalT1
                (*NodalT1)(i)=data1(2*i);
                (*NodalT1)(i+5)=data1(3*i+10);
                (*NodalT1)(i+10)=data1(3*i+11);
           }
          }
        }
        //for 15 data points
        double ThermalN, ThermalMz, ThermalMy;
        ThermalN =  0;
        ThermalMz = 0;
        ThermalMy = 0;

        double xi[maxNumSections];
        beamInt->getSectionLocations(numSections, L, xi);
        double wt[maxNumSections];
        beamInt->getSectionWeights(numSections, L, wt);

        // Loop over the integration points
        for (int i = 0; i < numSections; i++) {
          // Get section stress resultant
          Vector* dataMixV;
          if (NodalT0->Size()==9){
            for(int m=0;m<9;m++){
                (*dataMixV)(2*m)=(*NodalT0)(m)+xi[i]*((*NodalT1)(m) - (*NodalT0)(m)); //Linear temperature interpolation
                (*dataMixV)(2*m+1)=(*Loc)(m);
                (*dataMixV)(18+m)=1000;
            }
          }
          else if (NodalT0->Size()==15) {
             for (int m=0;m<5;m++){
                (*dataMixV)(2*m)=(*NodalT0)(m)+xi[i]*((*NodalT1)(m)-(*NodalT0)(m)); ////5 temps through y, Linear temperature interpolation
                (*dataMixV)(3*m+10)=(*NodalT0)(m+5)+xi[i]*((*NodalT1)(m+5)-(*NodalT0)(m+5));///////5 temps through Z in bottom flange
                (*dataMixV)(3*m+11)=(*NodalT0)(m+10)+xi[i]*((*NodalT1)(m+10)-(*NodalT0)(m+10));/////5 temps through Z in top flange
                (*dataMixV)(2*m+1)=(*Loc)(m);////5 (*Loc)s through y
                (*dataMixV)(3*m+12)=(*Loc)(m+5);///////5 (*Loc)s through Z
             }
          }
          const Vector &s = theSections[i]->getTemperatureStress(*dataMixV);

          SectionThermalElong[i]= theSections[i]->getThermalElong()(0);
          //apply temp along y
          residThermal[0] = -s(0);
          residThermal[1] = -s(1);
          residThermal[2] = s(1);
          residThermal[3] = -s(2);
          residThermal[4] = s(2);
          //apply temp along z
          //residThermal[0] = -s(0);
          //residThermal[1] = -0.;
          //residThermal[2] = -0.;
          //residThermal[3] = -s(1);
          //residThermal[4] = s(1);
          double ThermalEloni  = SectionThermalElong[i]*wt[i];
          AverageThermalElong += ThermalEloni;
        }
        counterTemperature = 1;
    }

    else if(type == LOAD_TAG_ThermalActionWrapper) {
        counterTemperature = 1;
        double ThermalN, ThermalMz, ThermalMy;
        ThermalN =  0;
        ThermalMz = 0;
        ThermalMy = 0;

        AverageThermalElong=0.0;

        double xi[maxNumSections];
        beamInt->getSectionLocations(numSections, L, xi);
        double wt[maxNumSections];
        beamInt->getSectionWeights(numSections, L, wt);

         Vector theNode0Crds = theNodes[0]->getCrds();
         Vector theNode1Crds = theNodes[1]->getCrds();
         int ndm = theNode0Crds.Size();
         Vector theIntCrds = Vector(ndm);

        // Loop over the integration points
         // Get section stress resultant
        for (int i = 0; i < numSections; i++) {
          int order = theSections[i]->getOrder();
          const ID &code = theSections[i]->getType();

          double xi6 = 6.0*xi[i];
          theIntCrds.Zero();
          for (int m = 0; m<ndm; m++){
              theIntCrds(m) =  theNode0Crds(m)+xi[i]*(theNode1Crds(m)- theNode0Crds(m));
          }

          //NodalThermalActions attached with the wrapper would have been updated by pattern;
          Vector dataMixV = ((ThermalActionWrapper*) theLoad)->getIntData(theIntCrds);

          const Vector &s = theSections[i]->getTemperatureStress(dataMixV);
          if (theSections[i]->getClassTag()==SEC_TAG_FiberSection3dThermal)
              SectionThermalElong[i] = ((FiberSection3dThermal*)theSections[i])->getThermalElong()(0);
          else if (theSections[i]->getClassTag()==SEC_TAG_FiberSectionGJThermal)
              SectionThermalElong[i] = ((FiberSectionGJThermal*)theSections[i])->getThermalElong()(0);


          //apply temp along y
          residThermal[0] = -s(0);
          residThermal[1] = -s(1);
          residThermal[2] = s(1);
          residThermal[3] = -s(2);
          residThermal[4] = s(2);

          //apply temp along z
          //residThermal[0] = -s(0);
          //residThermal[1] = -0.;
          //residThermal[2] = -0.;
          //residThermal[3] = -s(1);
          //residThermal[4] = s(1);
          AverageThermalElong += SectionThermalElong[i]*wt[i];
    }
        //end of for loop for section
        //counterTemperature = 1;
  }
#endif

  else {
    opserr << "CubicFrame3d::addLoad() -- load type unknown for element with tag: "
           << this->getTag() << "\n";
    return -1;
  }

  return 0;
}

int
CubicFrame3d::addInertiaLoadToUnbalance(const Vector& accel)
{
  // Check for a quick return
  if (density == 0.0)
    return 0;

  // Get R * accel from the nodes
  const Vector& Raccel1 = theNodes[0]->getRV(accel);
  const Vector& Raccel2 = theNodes[1]->getRV(accel);

  if (6 != Raccel1.Size() || 6 != Raccel2.Size()) {
    opserr << "CubicFrame3d::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
    return -1;
  }

  // want to add ( - fact * M R * accel ) to unbalance
  // take advantage of lumped mass matrix
  double L = theCoordTransf->getInitialLength();
  double m = 0.5 * density * L;

  Q(0) -= m * Raccel1(0);
  Q(1) -= m * Raccel1(1);
  Q(2) -= m * Raccel1(2);
  Q(6) -= m * Raccel2(0);
  Q(7) -= m * Raccel2(1);
  Q(8) -= m * Raccel2(2);

  return 0;
}

const Vector&
CubicFrame3d::getResistingForce()
{
  double L = theCoordTransf->getInitialLength();

  // Zero for integration
  q.Zero();

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {


    double xi6  = 6.0 * xi[i];
    double phiz = phizs[i];
    double phiy = phiys[i];

    // Get section stress resultant
    const VectorND<nsr> s  = theSections[i]->getResultant<nsr,scheme>();

    // Perform numerical integration on internal force
    //q.addMatrixTransposeVector(1.0, *B, s, wts(i));

    for (int j = 0; j < nsr; j++) {
      double si = s[j] * wt[i];
      switch (scheme[j]) {
      case SECTION_RESPONSE_P: q(0) += si; break;
      case SECTION_RESPONSE_MZ:
        q(1) += 1.0 / (1 + phiz) * (xi6 - 4.0 - phiz) * si;
        q(2) += 1.0 / (1 + phiz) * (xi6 - 2.0 + phiz) * si;
        break;
      case SECTION_RESPONSE_MY:
        q(3) += 1.0 / (1 + phiy) * (xi6 - 4.0 - phiy) * si;
        q(4) += 1.0 / (1 + phiy) * (xi6 - 2.0 + phiy) * si;
        break;
      case SECTION_RESPONSE_VY:
        q(1) += 0.5 * phiz * L / (1 + phiz) * si;
        q(2) += 0.5 * phiz * L / (1 + phiz) * si;
        break;
      case SECTION_RESPONSE_VZ:
        q(3) += 0.5 * phiy * L / (1 + phiy) * si;
        q(4) += 0.5 * phiy * L / (1 + phiy) * si;
        break;
      case SECTION_RESPONSE_T: q(5) += si; break;
      default:                 break;
      }
    }
  }

  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];
  q(3) += q0[3];
  q(4) += q0[4];

  // Transform forces
  Vector p0Vec(p0);
  P = theCoordTransf->getGlobalResistingForce(q, p0Vec);

  // Subtract other external nodal loads ... P_res = P_int - P_ext
  if (density != 0)
    P.addVector(1.0, Q, -1.0);

  return P;
}

const Vector&
CubicFrame3d::getResistingForceIncInertia()
{
  P = this->getResistingForce();

  if (density != 0.0) {
    const Vector& accel1 = theNodes[0]->getTrialAccel();
    const Vector& accel2 = theNodes[1]->getTrialAccel();

    // Compute the current resisting force
    this->getResistingForce();

    // take advantage of lumped mass matrix
    double L = theCoordTransf->getInitialLength();
    double m = 0.5 * density * L;

    P(0) += m * accel1(0);
    P(1) += m * accel1(1);
    P(2) += m * accel1(2);
    P(6) += m * accel2(0);
    P(7) += m * accel2(1);
    P(8) += m * accel2(2);

    // add the damping forces if rayleigh damping
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      P.addVector(1.0, this->getRayleighDampingForces(), 1.0);

  } else {

    // add the damping forces if rayleigh damping
    if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      P.addVector(1.0, this->getRayleighDampingForces(), 1.0);
  }

  return P;
}

int
CubicFrame3d::sendSelf(int commitTag, Channel& theChannel)
{
  // place the integer data into an ID

  int dbTag = this->getDbTag();
  int loc = 0;

  const ID& connectedExternalNodes = this->getExternalNodes();

  static Vector data(14);
  data(0)            = this->getTag();
  data(1)            = connectedExternalNodes(0);
  data(2)            = connectedExternalNodes(1);
  data(3)            = numSections;
  data(4)            = theCoordTransf->getClassTag();
  int crdTransfDbTag = theCoordTransf->getDbTag();
  if (crdTransfDbTag == 0) {
    crdTransfDbTag = theChannel.getDbTag();
    if (crdTransfDbTag != 0)
      theCoordTransf->setDbTag(crdTransfDbTag);
  }
  data(5)          = crdTransfDbTag;
  data(6)          = beamInt->getClassTag();
  int beamIntDbTag = beamInt->getDbTag();
  if (beamIntDbTag == 0) {
    beamIntDbTag = theChannel.getDbTag();
    if (beamIntDbTag != 0)
      beamInt->setDbTag(beamIntDbTag);
  }
  data(7) = beamIntDbTag;
  data(8) = density;
  //data(9) = cMass;
  data(10) = alphaM;
  data(11) = betaK;
  data(12) = betaK0;
  data(13) = betaKc;

  if (theChannel.sendVector(dbTag, commitTag, data) < 0) {
    opserr << "CubicFrame3d::sendSelf() - failed to send data Vector\n";
    return -1;
  }

  // send the coordinate transformation
  if (theCoordTransf->sendSelf(commitTag, theChannel) < 0) {
    opserr << "CubicFrame3d::sendSelf() - failed to send crdTranf\n";
    return -1;
  }

  // send the beam integration
  if (beamInt->sendSelf(commitTag, theChannel) < 0) {
    opserr << "CubicFrame3d::sendSelf() - failed to send beamInt\n";
    return -1;
  }

  //
  // send an ID for the sections containing each sections dbTag and classTag
  // if section ha no dbTag get one and assign it
  //

  ID idSections(2 * numSections);
  loc = 0;
  for (int i = 0; i < numSections; i++) {
    int sectClassTag = theSections[i]->getClassTag();
    int sectDbTag    = theSections[i]->getDbTag();
    if (sectDbTag == 0) {
      sectDbTag = theChannel.getDbTag();
      theSections[i]->setDbTag(sectDbTag);
    }

    idSections(loc)     = sectClassTag;
    idSections(loc + 1) = sectDbTag;
    loc += 2;
  }

  if (theChannel.sendID(dbTag, commitTag, idSections) < 0) {
    opserr << "CubicFrame3d::sendSelf() - failed to send ID data\n";
    return -1;
  }

  //
  // send the sections
  //

  for (int j = 0; j < numSections; j++) {
    if (theSections[j]->sendSelf(commitTag, theChannel) < 0) {
      opserr << "CubicFrame3d::sendSelf() - section " << j << "failed to send itself\n";
      return -1;
    }
  }

  return 0;
}

int
CubicFrame3d::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
  //
  // receive the integer data containing tag, numSections and coord transformation info
  //
  int dbTag = this->getDbTag();

  static Vector data(14);

  if (theChannel.recvVector(dbTag, commitTag, data) < 0) {
    opserr << "CubicFrame3d::recvSelf() - failed to recv data Vector\n";
    return -1;
  }

  this->setTag((int)data(0));
  connectedExternalNodes(0) = (int)data(1);
  connectedExternalNodes(1) = (int)data(2);
  int nSect                 = (int)data(3);
  int crdTransfClassTag     = (int)data(4);
  int crdTransfDbTag        = (int)data(5);

  int beamIntClassTag = (int)data(6);
  int beamIntDbTag    = (int)data(7);

  density = data(8);
  //cMass = (int)data(9);

  alphaM = data(10);
  betaK  = data(11);
  betaK0 = data(12);
  betaKc = data(13);

  // create a new crdTransf object if one needed
  if (theCoordTransf == 0 || theCoordTransf->getClassTag() != crdTransfClassTag) {
    if (theCoordTransf != 0)
      delete theCoordTransf;

//  theCoordTransf = theBroker.getNewCrdTransf(crdTransfClassTag);

    if (theCoordTransf == nullptr) {
      opserr << "CubicFrame3d::recvSelf() - " << "failed to obtain a CrdTrans object with classTag"
             << crdTransfClassTag << "\n";
      return -2;
    }
  }

  theCoordTransf->setDbTag(crdTransfDbTag);

  // invoke recvSelf on the crdTransf object
  if (theCoordTransf->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "CubicFrame3d::sendSelf() - failed to recv crdTranf\n";
    return -3;
  }

  // create a new beamInt object if one needed
  if (beamInt == 0 || beamInt->getClassTag() != beamIntClassTag) {
    if (beamInt != 0)
      delete beamInt;

    beamInt = theBroker.getNewBeamIntegration(beamIntClassTag);

    if (beamInt == 0) {
      opserr
          << "CubicFrame3d::recvSelf() - failed to obtain the beam integration object with classTag"
          << beamIntClassTag << "\n";
      exit(-1);
    }
  }

  beamInt->setDbTag(beamIntDbTag);

  // invoke recvSelf on the beamInt object
  if (beamInt->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "CubicFrame3d::sendSelf() - failed to recv beam integration\n";
    return -3;
  }

  //
  // recv an ID for the sections containing each sections dbTag and classTag
  //

  ID idSections(2 * nSect);
  int loc = 0;

  if (theChannel.recvID(dbTag, commitTag, idSections) < 0) {
    opserr << "CubicFrame3d::recvSelf() - failed to recv ID data\n";
    return -1;
  }

  //
  // now receive the sections
  //

  if (numSections != nSect) {

    //
    // we do not have correct number of sections, must delete the old and create
    // new ones before can recvSelf on the sections
    //

    // delete the old
    if (numSections != 0) {
      for (int i = 0; i < numSections; i++)
        delete theSections[i];
      delete[] theSections;
    }

    // create a new array to hold pointers
    theSections = new FrameSection*[nSect];
    if (theSections == 0) {
      opserr << "CubicFrame3d::recvSelf() - out of memory creating sections array of size" << nSect
             << "\n";
      exit(-1);
    }

    // create a section and recvSelf on it
    numSections = nSect;
    loc         = 0;

    for (int i = 0; i < numSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag    = idSections(loc + 1);
      loc += 2;
      // TODO
//    theSections[i] = theBroker.getNewSection(sectClassTag);
      if (theSections[i] == 0) {
        opserr << "CubicFrame3d::recvSelf() - Broker could not create Section of class type"
               << sectClassTag << "\n";
        exit(-1);
      }
      theSections[i]->setDbTag(sectDbTag);
      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
        opserr << "CubicFrame3d::recvSelf() - section " << i << "failed to recv itself\n";
        return -1;
      }
    }

  } else {

    //
    // for each existing section, check it is of correct type
    // (if not delete old & create a new one) then recvSelf on it
    //

    loc = 0;
    for (int i = 0; i < numSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag    = idSections(loc + 1);
      loc += 2;

      // check of correct type
      if (theSections[i]->getClassTag() != sectClassTag) {
        // delete the old section[i] and create a new one
        delete theSections[i];
//      theSections[i] = theBroker.getNewSection(sectClassTag);
        if (theSections[i] == nullptr) {
          opserr << "CubicFrame3d::recvSelf() - Broker could not create Section of class type"
                 << sectClassTag << "\n";
          exit(-1);
        }
      }

      // recvSelf on it
      theSections[i]->setDbTag(sectDbTag);
      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
        opserr << "CubicFrame3d::recvSelf() - section " << i << "failed to recv itself\n";
        return -1;
      }
    }
  }

  return 0;
}

void
CubicFrame3d::Print(OPS_Stream& s, int flag)
{
  const ID& node_tags = this->getExternalNodes();

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_ELEM_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"" << this->getClassType() << "\", ";

    s << "\"nodes\": [" << node_tags(0) << ", " 
                        << node_tags(1) << "]";
    s << ", ";

    s << "\"massperlength\": " << density;
    s << ", ";

    s << "\"sections\": [";
    for (int i = 0; i < numSections - 1; i++)
      s << theSections[i]->getTag() << ", ";
    s << theSections[numSections - 1]->getTag() << "]";
    s << ", ";

    s << "\"integration\": ";
    beamInt->Print(s, flag);
    s << ", ";

    s << "\"crdTransformation\": " << theCoordTransf->getTag() ;
    s << "}";
    return;
  }

  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "\nCubicFrame3d, element id:  " << this->getTag() << "\n";
    s << "\tConnected external nodes:  " << connectedExternalNodes;
    s << "\tCoordTransf: " << theCoordTransf->getTag() << "\n";
    s << "\tmass density:  " << density << "\n";

    double N, Mz1, Mz2, Vy, My1, My2, Vz, T;
    double L        = theCoordTransf->getInitialLength();
    double oneOverL = 1.0 / L;

    N   = q(0);
    Mz1 = q(1);
    Mz2 = q(2);
    Vy  = (Mz1 + Mz2) * oneOverL;
    My1 = q(3);
    My2 = q(4);
    Vz  = -(My1 + My2) * oneOverL;
    T   = q(5);

    s << "\tEnd 1 Forces (P Mz Vy My Vz T): " << -N + p0[0] << ' ' << Mz1 << ' ' << Vy + p0[1]
      << ' ' << My1 << ' ' << Vz + p0[3] << ' ' << -T << "\n";
    s << "\tEnd 2 Forces (P Mz Vy My Vz T): " << N << ' ' << Mz2 << ' ' << -Vy + p0[2] << ' ' << My2
      << ' ' << -Vz + p0[4] << ' ' << T << "\n";
    s << "Number of sections: " << numSections << "\n";
    beamInt->Print(s, flag);

    for (int i = 0; i < numSections; i++) {
      theSections[i]->Print(s, flag);
    }
  }

}


Response*
CubicFrame3d::setResponse(const char** argv, int argc, OPS_Stream& output)
{

  Response* theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType", "CubicFrame3d");
  output.attr("eleTag", this->getTag());
  output.attr("node1", connectedExternalNodes[0]);
  output.attr("node2", connectedExternalNodes[1]);

  //
  // we compare argv[0] for known response types
  //

  // global force -
  if (strcmp(argv[0], "forces") == 0 || strcmp(argv[0], "force") == 0 ||
      strcmp(argv[0], "globalForce") == 0 || strcmp(argv[0], "globalForces") == 0) {

    output.tag("ResponseType", "Px_1");
    output.tag("ResponseType", "Py_1");
    output.tag("ResponseType", "Pz_1");
    output.tag("ResponseType", "Mx_1");
    output.tag("ResponseType", "My_1");
    output.tag("ResponseType", "Mz_1");
    output.tag("ResponseType", "Px_2");
    output.tag("ResponseType", "Py_2");
    output.tag("ResponseType", "Pz_2");
    output.tag("ResponseType", "Mx_2");
    output.tag("ResponseType", "My_2");
    output.tag("ResponseType", "Mz_2");

    theResponse = new ElementResponse(this, 1, Vector(12));

  // Local force
  } else if (strcmp(argv[0], "localForce") == 0 || 
             strcmp(argv[0], "localForces") == 0) {

    output.tag("ResponseType", "N_1");
    output.tag("ResponseType", "Vy_1");
    output.tag("ResponseType", "Vz_1");
    output.tag("ResponseType", "T_1");
    output.tag("ResponseType", "My_1");
    output.tag("ResponseType", "Mz_1");
    output.tag("ResponseType", "N_2");
    output.tag("ResponseType", "Vy_2");
    output.tag("ResponseType", "Vz_2");
    output.tag("ResponseType", "T_2");
    output.tag("ResponseType", "My_2");
    output.tag("ResponseType", "Mz_2");

    theResponse = new ElementResponse(this, 2, Vector(12));

  } else if (strcmp(argv[0], "basicForce") == 0 || strcmp(argv[0], "basicForces") == 0) {
    output.tag("ResponseType", "N");
    output.tag("ResponseType", "M1");
    output.tag("ResponseType", "M2");

    theResponse = new ElementResponse(this, 9, Vector(6));
  } else if (strcmp(argv[0], "basicStiffness") == 0) {
    output.tag("ResponseType", "N");
    output.tag("ResponseType", "M1");
    output.tag("ResponseType", "M2");

    theResponse = new ElementResponse(this, 19, Matrix(6, 6));

  // chord rotation -
  } else if (strcmp(argv[0], "chordRotation") == 0 || strcmp(argv[0], "chordDeformation") == 0 ||
             strcmp(argv[0], "basicDeformation") == 0) {

    output.tag("ResponseType", "eps");
    output.tag("ResponseType", "thetaZ_1");
    output.tag("ResponseType", "thetaZ_2");
    output.tag("ResponseType", "thetaY_1");
    output.tag("ResponseType", "thetaY_2");
    output.tag("ResponseType", "thetaX");

    theResponse = new ElementResponse(this, 3, Vector(6));

  // 4: Plastic rotation
  } else if (strcmp(argv[0], "plasticRotation") == 0 ||
             strcmp(argv[0], "plasticDeformation") == 0) {

    output.tag("ResponseType", "epsP");
    output.tag("ResponseType", "thetaZP_1");
    output.tag("ResponseType", "thetaZP_2");
    output.tag("ResponseType", "thetaYP_1");
    output.tag("ResponseType", "thetaYP_2");
    output.tag("ResponseType", "thetaXP");

    theResponse = new ElementResponse(this, 4, Vector(6));


  } else if (strcmp(argv[0], "RayleighForces") == 0 || 
             strcmp(argv[0], "rayleighForces") == 0) {

    theResponse = new ElementResponse(this, 12, Vector(12));

  // 10-11: Integration
  } else if (strcmp(argv[0], "integrationPoints") == 0)
    theResponse = new ElementResponse(this, 10, Vector(numSections));

  else if (strcmp(argv[0], "integrationWeights") == 0)
    theResponse = new ElementResponse(this, 11, Vector(numSections));

  else if (strcmp(argv[0], "sectionTags") == 0)
    theResponse = new ElementResponse(this, 110, ID(numSections));

  // section response
  else if (strcmp(argv[0], "sectionX") == 0) {
    if (argc > 2) {
      float sectionLoc = atof(argv[1]);

      double xi[maxNumSections];
      double L = theCoordTransf->getInitialLength();
      beamInt->getSectionLocations(numSections, L, xi);

      sectionLoc /= L;

      float minDistance = fabs(xi[0] - sectionLoc);
      int sectionNum    = 0;
      for (int i = 1; i < numSections; i++) {
        if (fabs(xi[i] - sectionLoc) < minDistance) {
          minDistance = fabs(xi[i] - sectionLoc);
          sectionNum  = i;
        }
      }

      output.tag("GaussPointOutput");
      output.attr("number", sectionNum + 1);
      output.attr("eta", xi[sectionNum] * L);

      theResponse = theSections[sectionNum]->setResponse(&argv[2], argc - 2, output);
    }
  }

  else if (strcmp(argv[0], "section") == 0) {
    if (argc > 1) {

      int sectionNum = atoi(argv[1]);

      if (sectionNum > 0 && sectionNum <= numSections && argc > 2) {

        double xi[maxNumSections];
        double L = theCoordTransf->getInitialLength();
        beamInt->getSectionLocations(numSections, L, xi);

        output.tag("GaussPointOutput");
        output.attr("number", sectionNum);
        output.attr("eta", xi[sectionNum - 1] * L);

        theResponse = theSections[sectionNum - 1]->setResponse(&argv[2], argc - 2, output);

        output.endTag();
      } else if (sectionNum == 0) { // argv[1] was not an int, we want all sections,

        CompositeResponse* theCResponse = new CompositeResponse();
        int numResponse                 = 0;
        double xi[maxNumSections];
        double L = theCoordTransf->getInitialLength();
        beamInt->getSectionLocations(numSections, L, xi);

        for (int i = 0; i < numSections; i++) {

          output.tag("GaussPointOutput");
          output.attr("number", i + 1);
          output.attr("eta", xi[i] * L);

          Response* theSectionResponse = theSections[i]->setResponse(&argv[1], argc - 1, output);

          output.endTag();

          if (theSectionResponse != 0) {
            numResponse = theCResponse->addResponse(theSectionResponse);
          }
        }

        if (numResponse == 0) // no valid responses found
          delete theCResponse;
        else
          theResponse = theCResponse;
      }
    }
  }
  // by SAJalali
  else if (strcmp(argv[0], "energy") == 0) {
    theResponse = new ElementResponse(this, 13, 0.0);
  }

  if (theResponse == 0)
    theResponse = theCoordTransf->setResponse(argv, argc, output);

  output.endTag();
  return theResponse;
}

int
CubicFrame3d::getResponse(int responseID, Information& eleInfo)
{
  double N, V, M1, M2, T;
  double L        = theCoordTransf->getInitialLength();
  double oneOverL = 1.0 / L;

  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());

  else if (responseID == 12)
    return eleInfo.setVector(this->getRayleighDampingForces());

  else if (responseID == 2) {
    // Axial
    N    = q(0);
    P(6) = N;
    P(0) = -N + p0[0];

    // Torsion
    T    = q(5);
    P(9) = T;
    P(3) = -T;

    // Moments about z and shears along y
    M1    = q(1);
    M2    = q(2);
    P(5)  = M1;
    P(11) = M2;
    V     = (M1 + M2) * oneOverL;
    P(1)  = V + p0[1];
    P(7)  = -V + p0[2];

    // Moments about y and shears along z
    M1    = q(3);
    M2    = q(4);
    P(4)  = M1;
    P(10) = M2;
    V     = (M1 + M2) * oneOverL;
    P(2)  = -V + p0[3];
    P(8)  = V + p0[4];

    return eleInfo.setVector(P);
  }

  else if (responseID == 9) {
    return eleInfo.setVector(q);
  }

  else if (responseID == 19) {
    static Matrix kb(6, 6);
    this->getBasicStiff(kb);
    return eleInfo.setMatrix(kb);
  }

  // Chord rotation
  else if (responseID == 3)
    return eleInfo.setVector(theCoordTransf->getBasicTrialDisp());

  // Plastic rotation
  else if (responseID == 4) {
    static Vector vp(6);
    static Vector ve(6);
    static Matrix kb(6, 6);
    this->getBasicStiff(kb, 1);
    kb.Solve(q, ve);
    vp = theCoordTransf->getBasicTrialDisp();
    vp -= ve;
    return eleInfo.setVector(vp);
  }

  else if (responseID == 10) {
    double L = theCoordTransf->getInitialLength();
    double pts[maxNumSections];
    beamInt->getSectionLocations(numSections, L, pts);
    Vector locs(numSections);
    for (int i = 0; i < numSections; i++)
      locs(i) = pts[i] * L;
    return eleInfo.setVector(locs);
  }

  else if (responseID == 11) {
    double L = theCoordTransf->getInitialLength();
    double wts[maxNumSections];
    beamInt->getSectionWeights(numSections, L, wts);
    Vector weights(numSections);
    for (int i = 0; i < numSections; i++)
      weights(i) = wts[i] * L;
    return eleInfo.setVector(weights);
  }

  else if (responseID == 110) {
    ID tags(numSections);
    for (int i = 0; i < numSections; i++)
      tags(i) = theSections[i]->getTag();
    return eleInfo.setID(tags);
  }

  //by SAJalali
  else if (responseID == 13) {
    double xi[maxNumSections];
    double L = theCoordTransf->getInitialLength();
    beamInt->getSectionWeights(numSections, L, xi);
    double energy = 0;
    for (int i = 0; i < numSections; i++) {
      energy += theSections[i]->getEnergy() * xi[i] * L;
    }
    return eleInfo.setDouble(energy);
  }

  else
    return -1;
}


int
CubicFrame3d::setParameter(const char** argv, int argc, Parameter& param)
{
  if (argc < 1)
    return -1;

  // don't do anything if MaterialStageParameter calls this element
  if (strcmp(argv[0], "updateMaterialStage") == 0) {
    return -1;
  }

  // If the parameter belongs to the element itself
  if (strcmp(argv[0], "rho") == 0) {
    param.setValue(density);
    return param.addObject(1, this);
  }

  if (strstr(argv[0], "sectionX") != 0) {
    if (argc < 3)
      return -1;

    float sectionLoc = atof(argv[1]);

    double xi[maxNumSections];
    double L = theCoordTransf->getInitialLength();
    beamInt->getSectionLocations(numSections, L, xi);

    sectionLoc /= L;

    float minDistance = fabs(xi[0] - sectionLoc);
    int sectionNum    = 0;
    for (int i = 1; i < numSections; i++) {
      if (fabs(xi[i] - sectionLoc) < minDistance) {
        minDistance = fabs(xi[i] - sectionLoc);
        sectionNum  = i;
      }
    }
    return theSections[sectionNum]->setParameter(&argv[2], argc - 2, param);
  }
  // If the parameter belongs to a section or lower
  if (strstr(argv[0], "section") != 0) {

    if (argc < 3)
      return -1;

    // Get section number
    int sectionNum = atoi(argv[1]);

    if (sectionNum > 0 && sectionNum <= numSections)
      return theSections[sectionNum - 1]->setParameter(&argv[2], argc - 2, param);
    else
      return -1;
  }

  if (strstr(argv[0], "integration") != 0) {

    if (argc < 2)
      return -1;

    return beamInt->setParameter(&argv[1], argc - 1, param);
  }

  // Default, send to every object
  int ok     = 0;
  int result = -1;

  for (int i = 0; i < numSections; i++) {
    ok = theSections[i]->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }

  ok = beamInt->setParameter(argv, argc, param);
  if (ok != -1)
    result = ok;

  return result;
}

int
CubicFrame3d::updateParameter(int parameterID, Information& info)
{
  if (parameterID == 1) {
    density = info.theDouble;
    return 0;
  } else
    return -1;
}


int
CubicFrame3d::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;

  return 0;
}


const Matrix&
CubicFrame3d::getInitialStiffSensitivity(int gradNumber)
{
  thread_local MatrixND<12,12> dK{0.0};
  thread_local Matrix wrapper(dK);
  return wrapper;
}

const Matrix&
CubicFrame3d::getMassSensitivity(int gradNumber)
{
  K.Zero();

  if (density == 0.0 || parameterID != 1)
    return K;

  double L = theCoordTransf->getInitialLength();

  // Lumped mass matrix
  double m = 0.5 * L;
  K(0, 0) = K(1, 1) = K(2, 2) = K(6, 6) = K(7, 7) = K(8, 8) = m;

  return K;
}


const Vector&
CubicFrame3d::getResistingForceSensitivity(int gradNumber)
{
  double L   = theCoordTransf->getInitialLength();
  double jsx = 1.0 / L;


  // Zero for integration
  static Vector dqdh(6);
  dqdh.Zero();

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    int order      = theSections[i]->getOrder();
    const ID& code = theSections[i]->getType();

    double xi6  = 6.0 * xi[i];
    double phiz = phizs[i];
    double phiy = phiys[i];
    double wti  = wt[i];

    // Get section stress resultant gradient
    const Vector& dsdh = theSections[i]->getStressResultantSensitivity(gradNumber, true);

    // Perform numerical integration on internal force gradient
    double sensi;
    for (int j = 0; j < nsr; j++) {
      sensi = dsdh(j) * wti;
      switch (scheme[j]) {
      case SECTION_RESPONSE_P:
        dqdh(0) += sensi;
        break;
      case SECTION_RESPONSE_VY:
        dqdh(1) += 0.5 * phiz * L / (1 + phiz) * sensi;
        dqdh(2) += 0.5 * phiz * L / (1 + phiz) * sensi;
        break;
      case SECTION_RESPONSE_VZ:
        dqdh(3) += 0.5 * phiy * L / (1 + phiy) * sensi;
        dqdh(4) += 0.5 * phiy * L / (1 + phiy) * sensi;
        break;
      case SECTION_RESPONSE_MZ:
        dqdh(1) += 1.0 / (1 + phiz) * (xi6 - 4.0 - phiz) * sensi;
        dqdh(2) += 1.0 / (1 + phiz) * (xi6 - 2.0 + phiz) * sensi;
        break;
      case SECTION_RESPONSE_MY:
        dqdh(3) += 1.0 / (1 + phiy) * (xi6 - 4.0 - phiy) * sensi;
        dqdh(4) += 1.0 / (1 + phiy) * (xi6 - 2.0 + phiy) * sensi;
        break;
      case SECTION_RESPONSE_T: dqdh(5) += sensi; break;
      default:                 break;
      }
    }
  }

  // Transform forces
  static Vector dp0dh(6); // No distributed loads

  P.Zero();



  if (theCoordTransf->isShapeSensitivity()) {

    // Perform numerical integration to obtain basic stiffness matrix
    // Some extra declarations
    static Matrix kbmine(6, 6);
    kbmine.Zero();
    q.Zero();

    for (int i = 0; i < numSections; i++) {
      double tmp;

      int order      = theSections[i]->getOrder();

      double xi6  = 6.0 * xi[i];
      double phiz = phizs[i];
      double phiy = phiys[i];
      double wti  = wt[i];

      // Get the section tangent stiffness and stress resultant
      const MatrixND<nsr,nsr> ks = theSections[i]->getTangent<nsr,scheme>(State::Pres);
      const VectorND<nsr>     s  = theSections[i]->getResultant<nsr,scheme>();

      MatrixND<nsr,6> ka;
      ka.zero();

      for (int j = 0; j < nsr; j++) {
        double si = s(j) * wti;
        switch (scheme[j]) {
        case SECTION_RESPONSE_P:
          q(0) += si;
          for (int k = 0; k < nsr; k++)
            ka(k, 0) += ks(k, j) * wti;
          break;
        case SECTION_RESPONSE_MZ:
          q(1) += 1.0 / (1 + phiz) * (xi6 - 4.0 - phiz) * si;
          q(2) += 1.0 / (1 + phiz) * (xi6 - 2.0 + phiz) * si;
          for (int k = 0; k < nsr; k++) {
            tmp = ks(k, j) * wti;
            ka(k, 1) += 1.0 / (1 + phiz) * (xi6 - 4.0 - phiz) * tmp;
            ka(k, 2) += 1.0 / (1 + phiz) * (xi6 - 2.0 + phiz) * tmp;
          }
          break;
        case SECTION_RESPONSE_MY:
          q(3) += 1.0 / (1 + phiy) * (xi6 - 4.0 - phiy) * si;
          q(4) += 1.0 / (1 + phiy) * (xi6 - 2.0 + phiy) * si;
          for (int k = 0; k < nsr; k++) {
            tmp = ks(k, j) * wti;
            ka(k, 3) += 1.0 / (1 + phiy) * (xi6 - 4.0 - phiy) * tmp;
            ka(k, 4) += 1.0 / (1 + phiy) * (xi6 - 2.0 + phiy) * tmp;
          }
          break;
        case SECTION_RESPONSE_VY:
          q(1) += 0.5 * phiz * L / (1 + phiz) * si;
          q(2) += 0.5 * phiz * L / (1 + phiz) * si;
          for (int k = 0; k < nsr; k++) {
            tmp = ks(k, j) * wti;
            ka(k, 1) += 0.5 * phiz * L / (1 + phiz) * tmp;
            ka(k, 2) += 0.5 * phiz * L / (1 + phiz) * tmp;
          }
          break;
        case SECTION_RESPONSE_VZ:
          q(3) += 0.5 * phiy * L / (1 + phiy) * si;
          q(4) += 0.5 * phiy * L / (1 + phiy) * si;
          for (int k = 0; k < nsr; k++) {
            tmp = ks(k, j) * wti;
            ka(k, 3) += 0.5 * phiy * L / (1 + phiy) * tmp;
            ka(k, 4) += 0.5 * phiy * L / (1 + phiy) * tmp;
          }
          break;
        case SECTION_RESPONSE_T:
          q(5) += si;
          for (int k = 0; k < nsr; k++)
            ka(k, 5) += ks(k, j) * wti;

          break;
        default: break;
        }
      }
      for (int j = 0; j < nsr; j++) {
        switch (scheme[j]) {
        case SECTION_RESPONSE_P:
          for (int k = 0; k < 6; k++) {
            kbmine(0, k) += ka(j, k);
          }
          break;
        case SECTION_RESPONSE_MZ:
          for (int k = 0; k < 6; k++) {
            double tmp = ka(j, k);
            kbmine(1, k) += 1.0 / (1 + phiz) * (xi6 - 4.0 - phiz) * tmp;
            kbmine(2, k) += 1.0 / (1 + phiz) * (xi6 - 2.0 + phiz) * tmp;
          }
          break;
        case SECTION_RESPONSE_MY:
          for (int k = 0; k < 6; k++) {
            tmp = ka(j, k);
            kbmine(3, k) += 1.0 / (1 + phiy) * (xi6 - 4.0 - phiy) * tmp;
            kbmine(4, k) += 1.0 / (1 + phiy) * (xi6 - 2.0 + phiy) * tmp;
          }
          break;
        case SECTION_RESPONSE_VY:
          for (int k = 0; k < 6; k++) {
            tmp = ka(j, k);
            kbmine(1, k) += 0.5 * phiz * L / (1 + phiz) * tmp;
            kbmine(2, k) += 0.5 * phiz * L / (1 + phiz) * tmp;
          }
          break;
        case SECTION_RESPONSE_VZ:
          for (int k = 0; k < 6; k++) {
            tmp = ka(j, k);
            kbmine(3, k) += 0.5 * phiy * L / (1 + phiy) * tmp;
            kbmine(4, k) += 0.5 * phiy * L / (1 + phiy) * tmp;
          }
          break;
        case SECTION_RESPONSE_T:
          for (int k = 0; k < 6; k++) {
            kbmine(5, k) += ka(j, k);
          }
          break;
        default: break;
        }
      }
    }

    const Vector& A_u = theCoordTransf->getBasicTrialDisp();
    double dLdh       = theCoordTransf->getLengthGrad();
    double d1overLdh  = -dLdh / (L * L);
    // a^T k_s dadh v
    dqdh.addMatrixVector(1.0, kbmine, A_u, d1overLdh);



    // k dAdh u
    const Vector& dAdh_u = theCoordTransf->getBasicDisplFixedGrad();
    dqdh.addMatrixVector(1.0, kbmine, dAdh_u, jsx);

    // dAdh^T q
    P += theCoordTransf->getGlobalResistingForceShapeSensitivity(q, dp0dh, gradNumber);
  }

  // A^T (dqdh + k dAdh u)
  P += theCoordTransf->getGlobalResistingForce(dqdh, dp0dh);

  return P;
}


// NEW METHOD
int
CubicFrame3d::commitSensitivity(int gradNumber, int numGrads)
{
  // Get basic deformation and sensitivities
  const Vector& v = theCoordTransf->getBasicTrialDisp();

  static Vector dvdh(6);
  dvdh = theCoordTransf->getBasicDisplTotalGrad(gradNumber);

  double L        = theCoordTransf->getInitialLength();
  double oneOverL = 1.0 / L;

  // Some extra declarations
  double d1oLdh = theCoordTransf->getd1overLdh();

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    int order      = theSections[i]->getOrder();
    const ID& code = theSections[i]->getType();

    VectorND<nsr> e;

    double xi6 = 6.0 * xi[i];
    // Assume the phi values are constant
    double phiz = phizs[i];
    double phiy = phiys[i];

    for (int j = 0; j < nsr; j++) {
      switch (scheme[j]) {
      case SECTION_RESPONSE_P: e(j) = oneOverL * dvdh(0) + d1oLdh * v(0); break;
      case SECTION_RESPONSE_MZ:
        //e(j) = oneOverL*((xi6-4.0)*dvdh(1) + (xi6-2.0)*dvdh(2))
        //  + d1oLdh*((xi6-4.0)*v(1) + (xi6-2.0)*v(2));
        e(j) =
            oneOverL / (1 + phiz) * ((xi6 - 4.0 - phiz) * dvdh(1) + (xi6 - 2.0 + phiz) * dvdh(2)) +
            d1oLdh / (1 + phiz) * ((xi6 - 4.0 - phiz) * v(1) + (xi6 - 2.0 + phiz) * v(2));
        break;
      case SECTION_RESPONSE_MY:
        //e(j) = oneOverL*((xi6-4.0)*dvdh(3) + (xi6-2.0)*dvdh(4))
        //  + d1oLdh*((xi6-4.0)*v(3) + (xi6-2.0)*v(4));
        e(j) =
            oneOverL / (1 + phiy) * ((xi6 - 4.0 - phiy) * dvdh(3) + (xi6 - 2.0 + phiy) * dvdh(4)) +
            d1oLdh / (1 + phiy) * ((xi6 - 4.0 - phiy) * v(3) + (xi6 - 2.0 + phiy) * v(4));
        break;
      case SECTION_RESPONSE_VY:
        e(j) = 0.5 * phiz / (1 + phiz) * dvdh(1) + 0.5 * phiz / (1 + phiz) * dvdh(2);
        break;
      case SECTION_RESPONSE_VZ:
        e(j) = 0.5 * phiy / (1 + phiy) * dvdh(3) + 0.5 * phiy / (1 + phiy) * dvdh(4);
        break;
      case SECTION_RESPONSE_T: e(j) = oneOverL * dvdh(5) + d1oLdh * v(5); break;
      default:                 e(j) = 0.0; break;
      }
    }

    // Set the section deformations
    theSections[i]->commitSensitivity(e, gradNumber, numGrads);
  }

  return 0;
}

