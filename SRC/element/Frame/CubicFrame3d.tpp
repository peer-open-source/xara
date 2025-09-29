//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
// This element is adapted from TimoshenkoBeamColumn3d
//
// Written: CMP, MHS
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
#include <Information.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <CompositeResponse.h>
#include <ElementalLoad.h>
#include <BeamIntegration.h>
#include <Parameter.h>
#include <math.h>
#include <string.h>


#define ELE_TAG_CubicFrame3d 0

namespace OpenSees {

template <bool shear, int nwm>
CubicFrame3d<shear,nwm>::CubicFrame3d(int tag, 
                           std::array<int, 2>& nodes, 
                           std::vector<FrameSection*>& sections,
                           BeamIntegration& bi,
                           CrdTransf& coordTransf, 
                           double r)
 : Element(tag, ELE_TAG_CubicFrame3d),
   numSections(sections.size()),
   theSections(nullptr),
   theCoordTransf(nullptr),
   beamInt(nullptr),
   Q(12),
   q{},
   P{},
   P_wrap(P),
   K{},
   K_wrap(K),
   density(r),
   connectedExternalNodes(2),
   loads(nullptr),
   parameterID(0)
{
  // Allocate arrays of pointers to SectionForceDeformations
  theSections = new FrameSection*[numSections];

  for (int i = 0; i < numSections; i++) {
    theSections[i] = sections[i]->getFrameCopy(scheme);
  }

  beamInt = bi.getCopy();


  theCoordTransf = coordTransf.getCopy3d();

  connectedExternalNodes(0) = nodes[0];
  connectedExternalNodes(1) = nodes[1];


  theNodes[0] = nullptr;
  theNodes[1] = nullptr;
}


template <bool shear, int nwm>
CubicFrame3d<shear,nwm>::CubicFrame3d()
 : Element(0, ELE_TAG_CubicFrame3d),
   numSections(0),
   theSections(nullptr),
   theCoordTransf(nullptr),
   beamInt(nullptr),
   connectedExternalNodes(2),
   Q(12),
   P{},
   P_wrap(P),
   K{},
   K_wrap(K),
   density(0.0),
   loads(nullptr),
   parameterID(0)
{
  theNodes[0] = nullptr;
  theNodes[1] = nullptr;
}

template <bool shear, int nwm>
CubicFrame3d<shear,nwm>::~CubicFrame3d()
{
  if (theSections != nullptr)
    for (int i = 0; i < numSections; i++) {
      if (theSections[i] != nullptr)
        delete theSections[i];
    }

  // Delete the array of pointers to SectionForceDeformation pointer arrays
  delete[] theSections;

  if (theCoordTransf)
    delete theCoordTransf;

  if (beamInt != nullptr)
    delete beamInt;
}

template <bool shear, int nwm>
int
CubicFrame3d<shear,nwm>::getNumExternalNodes() const
{
  return NEN;
}

template <bool shear, int nwm>
const ID&
CubicFrame3d<shear,nwm>::getExternalNodes()
{
  return connectedExternalNodes;
}

template <bool shear, int nwm>
Node**
CubicFrame3d<shear,nwm>::getNodePtrs()
{
  return theNodes;
}

template <bool shear, int nwm>
int
CubicFrame3d<shear,nwm>::getNumDOF()
{
  return NDF*NEN;
}

template <bool shear, int nwm>
void
CubicFrame3d<shear,nwm>::setDomain(Domain* theDomain)
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
    return;
  }

  double L = theCoordTransf->getInitialLength();

  if (L == 0.0)
    return;

  for (int i = 0; i < numSections; i++) {
    const Matrix& ks0 = theSections[i]->getInitialTangent();
    const ID& code    = theSections[i]->getType();

    double EI = 0.0;
    double GA = 0.0;
    for (int k = 0; k < nsr; k++) {
      if (scheme[k] == FrameStress::Mz)
        EI += ks0(k, k);
      if (code(k) == FrameStress::Vy)
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
      if (code(k) == FrameStress::Vz)
        GA += ks0(k, k);
    }
    phiys[i] = 0.0;
    if (GA != 0.0)
      phiys[i] = 12 * EI / (GA * L * L);
  }

  beamInt->getSectionLocations(numSections, L, xi);
  beamInt->getSectionWeights(numSections, L, wt);

  if (theDomain != nullptr)
    this->Element::link(*theDomain);

  this->update();
}

template <bool shear, int nwm>
int
CubicFrame3d<shear,nwm>::commitState()
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

template <bool shear, int nwm>
int
CubicFrame3d<shear,nwm>::revertToLastCommit()
{
  int status = 0;

  // Loop over the integration points and revert to last committed state
  for (int i = 0; i < numSections; i++)
    status += theSections[i]->revertToLastCommit();

  status += theCoordTransf->revertToLastCommit();

  return status;
}

template <bool shear, int nwm>
int
CubicFrame3d<shear,nwm>::revertToStart()
{
  int status = 0;

  // Loop over the integration points and revert states to start
  for (int i = 0; i < numSections; i++)
    status += theSections[i]->revertToStart();

  status += theCoordTransf->revertToStart();

  return status;
}

template <bool shear, int nwm>
int
CubicFrame3d<shear,nwm>::update()
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


    double xi6  = 6.0 * xi[i];
    double phiz = phizs[i];
    double phiy = phiys[i];

    VectorND<nsr> e;
    for (int j = 0; j < nsr; j++) {
      switch (scheme[j]) {
      case FrameStress::N: 
        e(j) = jsx * v(0);
        break;
      case FrameStress::Vy:
        e(j) = 0.5 * phiz/(1 + phiz)*v(1) + 0.5 * phiz/(1 + phiz) * v(2);
        break;
      case FrameStress::Vz:
        e(j) = 0.5 * phiy/(1 + phiy)*v(3) + 0.5 * phiy/(1 + phiy) * v(4);
        break;
      case FrameStress::T: 
        e(j) = jsx * v(5); 
        break;
      case SECTION_RESPONSE_MY:
        e(j) = jsx / (1 + phiy) * ((xi6 - 4.0 - phiy) * v(3) + (xi6 - 2.0 + phiy) * v(4));
        break;
      case FrameStress::Mz:
        e(j) = jsx / (1 + phiz) * ((xi6 - 4.0 - phiz) * v(1) + (xi6 - 2.0 + phiz) * v(2));
        break;
      default:
        e(j) = 0.0; 
        break;
      }
    }

    // Set the section deformations
    err += theSections[i]->template setTrialState<nsr,scheme>(e);
  }

  return err;
}


template <bool shear, int nwm>
const Matrix&
CubicFrame3d<shear,nwm>::getTangentStiff()
{
  static MatrixND<nq,nq> kb;
  static Matrix wrapper(kb);


  double L   = theCoordTransf->getInitialLength();
  double jsx = 1.0 / L;


  kb.zero();
  q.zero();
  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    MatrixND<nsr,nq> ka;
    ka.zero();

    double xi6  = 6.0 * xi[i];
    double phiz = phizs[i];
    double phiy = phiys[i];

    // Get the section tangent stiffness and stress resultant
    MatrixND<nsr,nsr> ks = theSections[i]->template getTangent<nsr,scheme>(State::Pres);
    const VectorND<nsr> s = theSections[i]->template getResultant<nsr,scheme>();
    
    // Perform numerical integration
    double wti = wt[i] * jsx;
    for (int j = 0; j < nsr; j++) {
      double tmp;
      switch (scheme[j]) {
      case FrameStress::N:
        for (int k = 0; k < nsr; k++)
          ka(k, 0) += ks(k, j) * wti;
        break;
      case FrameStress::Mz:
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
      case FrameStress::Vy:
        for (int k = 0; k < nsr; k++) {
          tmp = ks(k, j) * wti;
          ka(k, 1) += 0.5 * phiz * L / (1 + phiz) * tmp;
          ka(k, 2) += 0.5 * phiz * L / (1 + phiz) * tmp;
        }
        break;
      case FrameStress::Vz:
        for (int k = 0; k < nsr; k++) {
          tmp = ks(k, j) * wti;
          ka(k, 3) += 0.5 * phiy * L / (1 + phiy) * tmp;
          ka(k, 4) += 0.5 * phiy * L / (1 + phiy) * tmp;
        }
        break;
      case FrameStress::T:
        for (int k = 0; k < nsr; k++)
          ka(k, 5) += ks(k, j) * wti;
        break;
      default: break;
      }
    }

    for (int j = 0; j < nsr; j++) {
      double tmp;
      switch (scheme[j]) {
      case FrameStress::N:
        for (int k = 0; k < 6; k++)
          kb(0, k) += ka(j, k);
        break;
      case FrameStress::Mz:
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
      case FrameStress::Vy:
        for (int k = 0; k < 6; k++) {
          tmp = ka(j, k);
          kb(1, k) += 0.5 * phiz * L / (1 + phiz) * tmp;
          kb(2, k) += 0.5 * phiz * L / (1 + phiz) * tmp;
        }
        break;
      case FrameStress::Vz:
        for (int k = 0; k < 6; k++) {
          tmp = ka(j, k);
          kb(3, k) += 0.5 * phiy * L / (1 + phiy) * tmp;
          kb(4, k) += 0.5 * phiy * L / (1 + phiy) * tmp;
        }
        break;
      case FrameStress::T:
        for (int k = 0; k < 6; k++)
          kb(5, k) += ka(j, k);
        break;
      default: break;
      }
    }

    // q.addMatrixTransposeVector(1.0, *B, s, wts(i));
    for (int j = 0; j < nsr; j++) {
      double si = s[j] * wt[i];
      switch (scheme[j]) {
      case FrameStress::N: 
        q(0) += si; break;
      case FrameStress::Mz:
        q(1) += 1.0 / (1 + phiz) * (xi6 - 4.0 - phiz) * si;
        q(2) += 1.0 / (1 + phiz) * (xi6 - 2.0 + phiz) * si;
        break;
      case SECTION_RESPONSE_MY:
        q(3) += 1.0 / (1 + phiy) * (xi6 - 4.0 - phiy) * si;
        q(4) += 1.0 / (1 + phiy) * (xi6 - 2.0 + phiy) * si;
        break;
      case FrameStress::Vy:
        q(1) += 0.5 * phiz * L / (1 + phiz) * si;
        q(2) += 0.5 * phiz * L / (1 + phiz) * si;
        break;
      case FrameStress::Vz:
        q(3) += 0.5 * phiy * L / (1 + phiy) * si;
        q(4) += 0.5 * phiy * L / (1 + phiy) * si;
        break;
      case FrameStress::T: 
        q(5) += si;
        break;
      default:
        break;
      }
    }
  }
  
  if (loads != nullptr) {
    loads->addBasicForce(&q[0]);
  }
  // q[0] += q0[0];
  // q[1] += q0[1];
  // q[2] += q0[2];
  // q[3] += q0[3];
  // q[4] += q0[4];

  // Transform to global stiffness
  return theCoordTransf->getGlobalStiffMatrix(wrapper, q);
}


template <bool shear, int nwm>
void
CubicFrame3d<shear,nwm>::getBasicStiff(Matrix& kb, int initial)
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

    //
    // Perform numerical integration
    //
    double wti = wt[i] * jsx;
    for (int j = 0; j < nsr; j++) {
      double tmp;
      switch (scheme[j]) {
      case FrameStress::N:
        for (int k = 0; k < nsr; k++)
          ka(k, 0) += ks(k, j) * wti;
        break;
      case FrameStress::Mz:
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
      case FrameStress::Vy:
        for (int k = 0; k < nsr; k++) {
          tmp = ks(k, j) * wti;
          ka(k, 1) += 0.5 * phiz * L / (1 + phiz) * tmp;
          ka(k, 2) += 0.5 * phiz * L / (1 + phiz) * tmp;
        }
        break;
      case FrameStress::Vz:
        for (int k = 0; k < nsr; k++) {
          tmp = ks(k, j) * wti;
          ka(k, 3) += 0.5 * phiy * L / (1 + phiy) * tmp;
          ka(k, 4) += 0.5 * phiy * L / (1 + phiy) * tmp;
        }
        break;
      case FrameStress::T:
        for (int k = 0; k < nsr; k++)
          ka(k, 5) += ks(k, j) * wti;
        break;
      default: break;
      }
    }

    for (int j = 0; j < nsr; j++) {
      switch (scheme[j]) {
      case FrameStress::N:
        for (int k = 0; k < 6; k++)
          kb(0, k) += ka(j, k);
        break;
      case FrameStress::Mz:
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
      case FrameStress::Vy:
        for (int k = 0; k < 6; k++) {
          double tmp = ka(j, k);
          kb(1, k) += 0.5 * phiz * L / (1 + phiz) * tmp;
          kb(2, k) += 0.5 * phiz * L / (1 + phiz) * tmp;
        }
        break;
      case FrameStress::Vz:
        for (int k = 0; k < 6; k++) {
          double tmp = ka(j, k);
          kb(3, k) += 0.5 * phiy * L / (1 + phiy) * tmp;
          kb(4, k) += 0.5 * phiy * L / (1 + phiy) * tmp;
        }
        break;
      case FrameStress::T:
        for (int k = 0; k < 6; k++)
          kb(5, k) += ka(j, k);
        break;
      default: break;
      }
    }
  }
}

template <bool shear, int nwm>
const Matrix&
CubicFrame3d<shear,nwm>::getInitialStiff()
{
  thread_local Matrix kb(nq, nq);

  this->getBasicStiff(kb, 1);

  return theCoordTransf->getInitialGlobalStiffMatrix(kb);
}

template <bool shear, int nwm>
const Matrix&
CubicFrame3d<shear,nwm>::getMass()
{
  // thread_local MatrixND<NEN*NDF,NEN*NDF> K;
  // thread_local Matrix Wrapper(K);
  K.zero();

  if (density == 0.0)
    return K_wrap;

  double L = theCoordTransf->getInitialLength();

  // lumped mass matrix
  double m = 0.5 * density * L;
  K(0, 0) = K(1, 1) = K(2, 2) = m;
  K(NDF+0, NDF+0) = K(NDF+1, NDF+1) = K(NDF+2, NDF+2) = m;

  return K_wrap;
}

template <bool shear, int nwm>
void
CubicFrame3d<shear,nwm>::zeroLoad()
{
  Q.Zero();
  if (loads != nullptr) {
    loads->zeroLoad();
  }

  // q0.zero();
  // p0.zero();

  return;
}


template <bool shear, int nwm>
int
CubicFrame3d<shear,nwm>::addLoad(ElementalLoad* theLoad, double loadFactor)
{
  if (loads == nullptr) {
    loads = new BasicFrame3d();
    loads->setLength(theCoordTransf->getInitialLength());
  }

  loads->addLoad(theLoad, loadFactor);
  return 0;
}

template <bool shear, int nwm>
int
CubicFrame3d<shear,nwm>::addInertiaLoadToUnbalance(const Vector& accel)
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

template <bool shear, int nwm>
const Vector&
CubicFrame3d<shear,nwm>::getResistingForce()
{
  double L = theCoordTransf->getInitialLength();

  // Zero for integration
  q.zero();

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    double xi6  = 6.0 * xi[i];
    double phiz = phizs[i];
    double phiy = phiys[i];

    // Get section stress resultant
    const VectorND<nsr> s  = theSections[i]->template getResultant<nsr,scheme>();

    // Perform numerical integration on internal force
    //q.addMatrixTransposeVector(1.0, *B, s, wts(i));

    for (int j = 0; j < nsr; j++) {
      double si = s[j] * wt[i];
      switch (scheme[j]) {
      case FrameStress::N: q(0) += si; break;
      case FrameStress::Mz:
        q(1) += 1.0 / (1 + phiz) * (xi6 - 4.0 - phiz) * si;
        q(2) += 1.0 / (1 + phiz) * (xi6 - 2.0 + phiz) * si;
        break;
      case SECTION_RESPONSE_MY:
        q(3) += 1.0 / (1 + phiy) * (xi6 - 4.0 - phiy) * si;
        q(4) += 1.0 / (1 + phiy) * (xi6 - 2.0 + phiy) * si;
        break;
      case FrameStress::Vy:
        q(1) += 0.5 * phiz * L / (1 + phiz) * si;
        q(2) += 0.5 * phiz * L / (1 + phiz) * si;
        break;
      case FrameStress::Vz:
        q(3) += 0.5 * phiy * L / (1 + phiy) * si;
        q(4) += 0.5 * phiy * L / (1 + phiy) * si;
        break;
      case FrameStress::T: q(5) += si; break;
      default:
        break;
      }
    }
  }

  // q(0) += q0[0];
  // q(1) += q0[1];
  // q(2) += q0[2];
  // q(3) += q0[3];
  // q(4) += q0[4];
  if (loads != nullptr) {
    loads->addBasicForce(&q[0]);
    P = theCoordTransf->getGlobalResistingForce(q, loads->getReactions());
  } 
  else {
    static Vector p0(6);
    P = theCoordTransf->getGlobalResistingForce(q, p0);
  }

  // Subtract other external nodal loads
  // P_res = P_int - P_ext
  if (density != 0)
    P.addVector(1.0, Q, -1.0);

  return P_wrap;
}

template <bool shear, int nwm>
const Vector&
CubicFrame3d<shear,nwm>::getResistingForceIncInertia()
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

  return P_wrap;
}


template <bool shear, int nwm>
int
CubicFrame3d<shear,nwm>::sendSelf(int commitTag, Channel& theChannel)
{
  return -1;
}

template <bool shear, int nwm>
int
CubicFrame3d<shear,nwm>::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
  return -1;
}


template <bool shear, int nwm>
void
CubicFrame3d<shear,nwm>::Print(OPS_Stream& s, int flag)
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

    for (int i = 0; i < numSections; i++) {
      theSections[i]->Print(s, flag);
    }
  }

}


template <bool shear, int nwm>
Response*
CubicFrame3d<shear,nwm>::setResponse(const char** argv, int argc, OPS_Stream& output)
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

    theResponse = new ElementResponse(this, 9, Vector(nq));
  } else if (strcmp(argv[0], "basicStiffness") == 0) {
    output.tag("ResponseType", "N");
    output.tag("ResponseType", "M1");
    output.tag("ResponseType", "M2");

    theResponse = new ElementResponse(this, 19, Matrix(nq, nq));

  // chord rotation -
  } else if (strcmp(argv[0], "chordRotation") == 0 || strcmp(argv[0], "chordDeformation") == 0 ||
             strcmp(argv[0], "basicDeformation") == 0) {

    output.tag("ResponseType", "eps");
    output.tag("ResponseType", "thetaZ_1");
    output.tag("ResponseType", "thetaZ_2");
    output.tag("ResponseType", "thetaY_1");
    output.tag("ResponseType", "thetaY_2");
    output.tag("ResponseType", "thetaX");

    theResponse = new ElementResponse(this, 3, Vector(nq));

  // 4: Plastic rotation
  } else if (strcmp(argv[0], "plasticRotation") == 0 ||
             strcmp(argv[0], "plasticDeformation") == 0) {

    output.tag("ResponseType", "epsP");
    output.tag("ResponseType", "thetaZP_1");
    output.tag("ResponseType", "thetaZP_2");
    output.tag("ResponseType", "thetaYP_1");
    output.tag("ResponseType", "thetaYP_2");
    output.tag("ResponseType", "thetaXP");

    theResponse = new ElementResponse(this, 4, Vector(nq));


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

          if (theSectionResponse != nullptr) {
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

  if (theResponse == nullptr)
    theResponse = theCoordTransf->setResponse(argv, argc, output);

  output.endTag();
  return theResponse;
}

template <bool shear, int nwm>
int
CubicFrame3d<shear,nwm>::getResponse(int responseID, Information& eleInfo)
{
  double N, V, M1, M2, T;
  double L        = theCoordTransf->getInitialLength();
  double oneOverL = 1.0 / L;

  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());

  else if (responseID == 12)
    return eleInfo.setVector(this->getRayleighDampingForces());

  else if (responseID == 2) {
    Vector p0(6);
    if (loads != nullptr) {
      p0 = loads->getReactions();
    } else {
      p0.Zero();
    }

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
    P(1)  =  V + p0[1];
    P(7)  = -V + p0[2];

    // Moments about y and shears along z
    M1    = q(3);
    M2    = q(4);
    P(4)  = M1;
    P(10) = M2;
    V     = (M1 + M2) * oneOverL;
    P(2)  = -V + p0[3];
    P(8)  = V + p0[4];

    return eleInfo.setVector(P_wrap);
  }

  else if (responseID == 9) {
    return eleInfo.setVector(q);
  }

  else if (responseID == 19) {
    static Matrix kb(nq, nq);
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


template <bool shear, int nwm>
int
CubicFrame3d<shear,nwm>::setParameter(const char** argv, int argc, Parameter& param)
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

template <bool shear, int nwm>
int
CubicFrame3d<shear,nwm>::updateParameter(int parameterID, Information& info)
{
  if (parameterID == 1) {
    density = info.theDouble;
    return 0;
  } else
    return -1;
}


template <bool shear, int nwm>
int
CubicFrame3d<shear,nwm>::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;
  return 0;
}


template <bool shear, int nwm>
const Matrix&
CubicFrame3d<shear,nwm>::getInitialStiffSensitivity(int gradNumber)
{
  thread_local MatrixND<NEN*NDF,NEN*NDF> dK{};
  thread_local Matrix wrapper(dK);
  return wrapper;
}

template <bool shear, int nwm>
const Matrix&
CubicFrame3d<shear,nwm>::getMassSensitivity(int gradNumber)
{
  K.zero();

  if (density == 0.0 || parameterID != 1)
    return K_wrap;

  double L = theCoordTransf->getInitialLength();

  // Lumped mass matrix
  double m = 0.5 * L;
  K(0, 0) = K(1, 1) = K(2, 2) = m;
  K(6, 6) = K(7, 7) = K(8, 8) = m;

  return K_wrap;
}


template <bool shear, int nwm>
const Vector&
CubicFrame3d<shear,nwm>::getResistingForceSensitivity(int gradNumber)
{
  double L   = theCoordTransf->getInitialLength();
  double jsx = 1.0 / L;


  // Zero for integration
  static Vector dqdh(6);
  dqdh.Zero();

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    double xi6  = 6.0 * xi[i];
    double phiz = phizs[i];
    double phiy = phiys[i];
    double wti  = wt[i];

    // Get section stress resultant gradient
    const VectorND<nsr> ds = theSections[i]->template getResultantGradient<nsr,scheme>(gradNumber, true);

    // Perform numerical integration on internal force gradient
    for (int j = 0; j < nsr; j++) {
      double sensi = ds(j) * wti;
      switch (scheme[j]) {
      case FrameStress::N:
        dqdh(0) += sensi;
        break;
      case FrameStress::Vy:
        dqdh(1) += 0.5 * phiz * L / (1 + phiz) * sensi;
        dqdh(2) += 0.5 * phiz * L / (1 + phiz) * sensi;
        break;
      case FrameStress::Vz:
        dqdh(3) += 0.5 * phiy * L / (1 + phiy) * sensi;
        dqdh(4) += 0.5 * phiy * L / (1 + phiy) * sensi;
        break;
      case FrameStress::Mz:
        dqdh(1) += 1.0 / (1 + phiz) * (xi6 - 4.0 - phiz) * sensi;
        dqdh(2) += 1.0 / (1 + phiz) * (xi6 - 2.0 + phiz) * sensi;
        break;
      case SECTION_RESPONSE_MY:
        dqdh(3) += 1.0 / (1 + phiy) * (xi6 - 4.0 - phiy) * sensi;
        dqdh(4) += 1.0 / (1 + phiy) * (xi6 - 2.0 + phiy) * sensi;
        break;
      case FrameStress::T: 
        dqdh(5) += sensi; 
        break;
      default:
        break;
      }
    }
  }

  // Transform forces
  static Vector dp0dh(nq); // No distributed loads

  P.zero();



  if (theCoordTransf->isShapeSensitivity()) {

    // Perform numerical integration to obtain basic stiffness matrix
    // Some extra declarations
    static Matrix kbmine(nq, nq);
    kbmine.Zero();
    q.zero();

    for (int i = 0; i < numSections; i++) {
      double tmp;

      double xi6  = 6.0 * xi[i];
      double phiz = phizs[i];
      double phiy = phiys[i];
      double wti  = wt[i];

      // Get the section tangent stiffness and stress resultant
      const MatrixND<nsr,nsr> ks = theSections[i]->template getTangent<nsr,scheme>(State::Pres);
      const VectorND<nsr>     s  = theSections[i]->template getResultant<nsr,scheme>();

      MatrixND<nsr,nq> ka;
      ka.zero();

      for (int j = 0; j < nsr; j++) {
        double si = s(j) * wti;
        switch (scheme[j]) {
        case FrameStress::N:
          q(0) += si;
          for (int k = 0; k < nsr; k++)
            ka(k, 0) += ks(k, j) * wti;
          break;
        case FrameStress::Mz:
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
        case FrameStress::Vy:
          q(1) += 0.5 * phiz * L / (1 + phiz) * si;
          q(2) += 0.5 * phiz * L / (1 + phiz) * si;
          for (int k = 0; k < nsr; k++) {
            tmp = ks(k, j) * wti;
            ka(k, 1) += 0.5 * phiz * L / (1 + phiz) * tmp;
            ka(k, 2) += 0.5 * phiz * L / (1 + phiz) * tmp;
          }
          break;
        case FrameStress::Vz:
          q(3) += 0.5 * phiy * L / (1 + phiy) * si;
          q(4) += 0.5 * phiy * L / (1 + phiy) * si;
          for (int k = 0; k < nsr; k++) {
            tmp = ks(k, j) * wti;
            ka(k, 3) += 0.5 * phiy * L / (1 + phiy) * tmp;
            ka(k, 4) += 0.5 * phiy * L / (1 + phiy) * tmp;
          }
          break;
        case FrameStress::T:
          q(5) += si;
          for (int k = 0; k < nsr; k++)
            ka(k, 5) += ks(k, j) * wti;

          break;
        default: break;
        }
      }

      for (int j = 0; j < nsr; j++) {
        switch (scheme[j]) {
        case FrameStress::N:
          for (int k = 0; k < 6; k++) {
            kbmine(0, k) += ka(j, k);
          }
          break;
        case FrameStress::Vy:
          for (int k = 0; k < 6; k++) {
            tmp = ka(j, k);
            kbmine(1, k) += 0.5 * phiz * L / (1 + phiz) * tmp;
            kbmine(2, k) += 0.5 * phiz * L / (1 + phiz) * tmp;
          }
          break;
        case FrameStress::Vz:
          for (int k = 0; k < 6; k++) {
            tmp = ka(j, k);
            kbmine(3, k) += 0.5 * phiy * L / (1 + phiy) * tmp;
            kbmine(4, k) += 0.5 * phiy * L / (1 + phiy) * tmp;
          }
          break;

        case SECTION_RESPONSE_MY:
          for (int k = 0; k < 6; k++) {
            tmp = ka(j, k);
            kbmine(3, k) += 1.0 / (1 + phiy) * (xi6 - 4.0 - phiy) * tmp;
            kbmine(4, k) += 1.0 / (1 + phiy) * (xi6 - 2.0 + phiy) * tmp;
          }
          break;
        case FrameStress::Mz:
          for (int k = 0; k < 6; k++) {
            double tmp = ka(j, k);
            kbmine(1, k) += 1.0 / (1 + phiz) * (xi6 - 4.0 - phiz) * tmp;
            kbmine(2, k) += 1.0 / (1 + phiz) * (xi6 - 2.0 + phiz) * tmp;
          }
          break;
        case FrameStress::T:
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

  return P_wrap;
}


// NEW METHOD
template <bool shear, int nwm>
int
CubicFrame3d<shear,nwm>::commitSensitivity(int gradNumber, int numGrads)
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

    VectorND<nsr> e;

    double xi6 = 6.0 * xi[i];
    // Assume the phi values are constant
    double phiz = phizs[i];
    double phiy = phiys[i];

    for (int j = 0; j < nsr; j++) {
      switch (scheme[j]) {
      case FrameStress::N: e(j) = oneOverL * dvdh(0) + d1oLdh * v(0); break;
      case FrameStress::Mz:
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
      case FrameStress::Vy:
        e(j) = 0.5 * phiz / (1 + phiz) * dvdh(1) + 0.5 * phiz / (1 + phiz) * dvdh(2);
        break;
      case FrameStress::Vz:
        e(j) = 0.5 * phiy / (1 + phiy) * dvdh(3) + 0.5 * phiy / (1 + phiy) * dvdh(4);
        break;
      case FrameStress::T: e(j) = oneOverL * dvdh(5) + d1oLdh * v(5); break;
      default:                 e(j) = 0.0; break;
      }
    }

    // Set the section deformations
    theSections[i]->commitSensitivity(e, gradNumber, numGrads);
  }

  return 0;
}

} // namespace OpenSees
