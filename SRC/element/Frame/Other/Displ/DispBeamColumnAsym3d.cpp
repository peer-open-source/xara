/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
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
// Description: Adapted for analysis of asymmetric sections with introducing
// high-order axial terms for the basic element formulation
// References:
//
// Du, X., & Hajjar, J. (2021). Three-dimensional nonlinear displacement-based beam element
// for members with angle and tee sections. Engineering Structures, 239, 112239.
//
// Written: Xinlong Du and Jerome F. Hajjar, Northeastern University, USA; 
// Created: 2019
// Adapted from: MHS, Feb 2001
//
#include <DispBeamColumnAsym3d.h>
#include <Node.h>
#include <SectionForceDeformation.h>
#include <CrdTransf.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Renderer.h>
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
#include <elementAPI.h>

#include <LobattoBeamIntegration.h>
#include <LegendreBeamIntegration.h>
#include <RadauBeamIntegration.h>
#include <NewtonCotesBeamIntegration.h>
#include <TrapezoidalBeamIntegration.h>
#include <RegularizedHingeIntegration.h>

Matrix DispBeamColumnAsym3d::K(12, 12);
Vector DispBeamColumnAsym3d::P(12);
double DispBeamColumnAsym3d::workArea[200];


DispBeamColumnAsym3d::DispBeamColumnAsym3d(int tag, int nd1, int nd2, int numSec,
                                           SectionForceDeformation **s,
                                           BeamIntegration &bi, CrdTransf &coordTransf,
                                           double yss, double zss, double r,
                                           int cm) //Xinlong
    : Element(tag, ELE_TAG_DispBeamColumnAsym3d), numSections(numSec), 
      theSections(nullptr), crdTransf(nullptr), beamInt(nullptr), 
      connectedExternalNodes(2), 
      Q(12), q(6), ys(yss), zs(zss),
      rho(r), cMass(cm), parameterID(0) //Xinlong
{
  // Allocate arrays of pointers to SectionForceDeformations
  theSections = new SectionForceDeformation *[numSections];


  for (int i = 0; i < numSections; i++) {

    // Get copies of the material model for each integration point
    theSections[i] = s[i]->getCopy();

    // Check allocation
    if (theSections[i] == 0) {
      opserr << "DispBeamColumnAsym3d::DispBeamColumnAsym3d -- failed to get a copy of "
                "section model\n";
      exit(-1);
    }
  }

  beamInt = bi.getCopy();

  crdTransf = coordTransf.getCopy3d();

  // Set connected external node IDs
  connectedExternalNodes(0) = nd1;
  connectedExternalNodes(1) = nd2;

  q0.zero();
  p0.zero();

  theNodes[0] = nullptr;
  theNodes[1] = nullptr;
}

DispBeamColumnAsym3d::DispBeamColumnAsym3d()
    : Element(0, ELE_TAG_DispBeamColumnAsym3d), numSections(0), theSections(0),
      crdTransf(0), beamInt(0), connectedExternalNodes(2), Q(12), q(6), ys(0.0), zs(0.0),
      rho(0.0), cMass(0), parameterID(0) //Xinlong
{
  q0.zero();
  p0.zero();

  theNodes[0] = nullptr;
  theNodes[1] = nullptr;
}

DispBeamColumnAsym3d::~DispBeamColumnAsym3d()
{
  for (int i = 0; i < numSections; i++) {
    if (theSections[i])
      delete theSections[i];
  }

  // Delete the array of pointers to SectionForceDeformation pointer arrays
  if (theSections)
    delete[] theSections;

  if (crdTransf)
    delete crdTransf;

  if (beamInt != nullptr)
    delete beamInt;
}

int
DispBeamColumnAsym3d::getNumExternalNodes() const
{
  return 2;
}

const ID &
DispBeamColumnAsym3d::getExternalNodes()
{
  return connectedExternalNodes;
}

Node **
DispBeamColumnAsym3d::getNodePtrs()
{

  return theNodes;
}

int
DispBeamColumnAsym3d::getNumDOF()
{
  return 12;
}

void
DispBeamColumnAsym3d::setDomain(Domain *theDomain)
{
  // Check Domain is not null - invoked when object removed from a domain
  if (theDomain == 0) {
    theNodes[0] = 0;
    theNodes[1] = 0;
    return;
  }

  int Nd1 = connectedExternalNodes(0);
  int Nd2 = connectedExternalNodes(1);

  theNodes[0] = theDomain->getNode(Nd1);
  theNodes[1] = theDomain->getNode(Nd2);

  if (theNodes[0] == 0 || theNodes[1] == 0) {
    //opserr << "FATAL ERROR DispBeamColumnAsym3d (tag: %d), node not found in domain",
    //	this->getTag());

    return;
  }

  int dofNd1 = theNodes[0]->getNumberDOF();
  int dofNd2 = theNodes[1]->getNumberDOF();

  if (dofNd1 != 6 || dofNd2 != 6) {
    //opserr << "FATAL ERROR DispBeamColumnAsym3d (tag: %d), has differing number of DOFs at its nodes",
    //	this->getTag());

    return;
  }

  if (crdTransf->initialize(theNodes[0], theNodes[1])) {
    // Add some error check
  }

  double L = crdTransf->getInitialLength();

  if (L == 0.0) {
    // Add some error check
  }

  this->DomainComponent::setDomain(theDomain);

  this->update();
}

int
DispBeamColumnAsym3d::commitState()
{
  int retVal = 0;

  // call element commitState to do any base class stuff
  if ((retVal = this->Element::commitState()) != 0) {
    opserr << "DispBeamColumnAsym3d::commitState () - failed in base class";
  }

  // Loop over the integration points and commit the material states
  for (int i = 0; i < numSections; i++)
    retVal += theSections[i]->commitState();

  retVal += crdTransf->commitState();

  return retVal;
}

int
DispBeamColumnAsym3d::revertToLastCommit()
{
  int retVal = 0;

  // Loop over the integration points and revert to last committed state
  for (int i = 0; i < numSections; i++)
    retVal += theSections[i]->revertToLastCommit();

  retVal += crdTransf->revertToLastCommit();

  return retVal;
}

int
DispBeamColumnAsym3d::revertToStart()
{
  int retVal = 0;

  // Loop over the integration points and revert states to start
  for (int i = 0; i < numSections; i++)
    retVal += theSections[i]->revertToStart();

  retVal += crdTransf->revertToStart();

  return retVal;
}

int
DispBeamColumnAsym3d::update()
{
  int err = 0;

  // Update the transformation
  crdTransf->update();

  // Get basic deformations
  const Vector &v = crdTransf->getBasicTrialDisp();

  double L        = crdTransf->getInitialLength();
  double oneOverL = 1.0 / L;

  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    //int order = theSections[i]->getOrder();
    //const ID &code = theSections[i]->getType();

    Vector e(workArea, 5);

    double xi1   = xi[i];
    double dNv1  = 1.0 + 3.0*xi1*xi1 - 4.0*xi1;
    double ddNv1 = 6.0 * xi1 * oneOverL - 4.0 * oneOverL;
    double dNv2  = 3.0 * xi1 * xi1 - 2.0 * xi1;
    double ddNv2 = 6.0 * xi1 * oneOverL - 2.0 * oneOverL;
    double dNw1  = -dNv1;
    double ddNw1 = -ddNv1;
    double dNw2  = -dNv2;
    double ddNw2 = -ddNv2;
    double Nf1   = xi1;

    // generalized strain
    double dx   = oneOverL*v(0);            // u'
    double dy   =  dNw1*v(3) +  dNw2*v(4);  // y'
    double dz   =  dNv1*v(1) +  dNv2*v(2);  // z'
    double phi  = Nf1*v(5);                 // phi
                                            //
    double ddy  = ddNw1*v(3) + ddNw2*v(4);  // y"
    double ddz  = ddNv1*v(1) + ddNv2*v(2);  // z"
    double dphi = oneOverL*v(5);            // phi'
    double s7   = v(1);                     // theta_Iz
    double s8   = v(3);                     // theta_Iy
    double s9   = v(2);                     // theta_Jz
    double ddz0 = v(4);                     // theta_Jy

    // section deformation to be sent to FiberSection
    e(0) = dx
         + (4.0*s7*s7 + 4.0*s8*s8 + 4.0*s9*s9 + 4.0*ddz0*ddz0 -
            2.0*s7*s9 - 2.0*s8*ddz0) / 60.0
         + (zs*dz - ys*dy) * dphi;
    e(1) =  ddz + ddy*phi;
    e(2) = -ddy + ddz*phi;
    e(3) = 0.5 * dphi*dphi;
    e(4) = dphi;

    // Set the section deformations
    err += theSections[i]->setTrialSectionDeformation(e);
  }

  if (err != 0) {
    opserr << "DispBeamColumnAsym3d::update() - failed setTrialSectionDeformations()\n";
    return err;
  }
  return 0;
}

const Matrix &
DispBeamColumnAsym3d::getTangentStiff()
{
  static Matrix kb(6, 6);
  static Matrix N1(5, 11);  //Xinlong
  static Matrix N2(11, 6);  //Xinlong
  static Matrix N3(11, 11); //Xinlong
  static Matrix kbPart1(6, 6);
  static Matrix Gm(11, 11);
  static Matrix kbPart2(6, 6);
  static Matrix Tr(6, 6);
  static Matrix kf1(6, 6);
  static Matrix kf2(6, 6);

  const Vector &v = crdTransf->getBasicTrialDisp();

  // Zero for integral
  kb.Zero();
  q.Zero();

  double L        = crdTransf->getInitialLength();
  double oneOverL = 1.0 / L;

  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  //const Vector &wts = quadRule.getIntegrPointWeights(numSections);
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    N1.Zero();
    N2.Zero();
    N3.Zero();
    kbPart1.Zero();
    kbPart2.Zero();
    Gm.Zero();

    Tr.Zero();
    kf1.Zero();
    kf2.Zero();

    double xi1   = xi[i];
    double dNv1  = 1.0 + 3.0 * xi1 * xi1 - 4.0 * xi1;
    double ddNv1 = 6.0 * xi1 * oneOverL - 4.0 * oneOverL;
    double dNv2  = 3.0 * xi1 * xi1 - 2.0 * xi1;
    double ddNv2 = 6.0 * xi1 * oneOverL - 2.0 * oneOverL;
    double dNw1  = -dNv1;
    double ddNw1 = -ddNv1;
    double dNw2  = -dNv2;
    double ddNw2 = -ddNv2;
    double Nf1   = xi1;

    double dv  =  dNv1*v(1) +  dNv2*v(2); // v'
    double ddv = ddNv1*v(1) + ddNv2*v(2); // v"
    double dw  =  dNw1*v(3) +  dNw2*v(4); // w'
    double ddw = ddNw1*v(3) + ddNw2*v(4); // w"
    double f   = Nf1*v(5);                // phi
    double df  = oneOverL*v(5);           // phi'

    //Matrix Ndeltad1-------------------------------------------------------
    N1(0, 0)  = 1.0;
    N1(0, 1)  = (4.0 * v(1) - v(2)) / 30.0;
    N1(0, 2)  = (4.0 * v(3) - v(4)) / 30.0;
    N1(0, 3)  = (4.0 * v(2) - v(1)) / 30.0;
    N1(0, 4)  = (4.0 * v(4) - v(3)) / 30.0;
    N1(0, 5)  = zs * df;
    N1(0, 6)  = -ys * df;
    N1(0, 10) = zs * dv - ys * dw;
    N1(1, 7)  = 1.0;
    N1(1, 8)  = f;
    N1(1, 9)  = ddw;
    N1(2, 7)  = f;
    N1(2, 8)  = -1.0;
    N1(2, 9)  = ddv;
    N1(3, 10) = df;
    N1(4, 10) = 1.0;

    //Matrix Ndeltad2-------------------------------------------------------
    N2(0, 0)  = oneOverL;
    N2(1, 1)  = 1.0;
    N2(2, 3)  = 1.0;
    N2(3, 2)  = 1.0;
    N2(4, 4)  = 1.0;
    N2(5, 1)  = dNv1;
    N2(5, 2)  = dNv2;
    N2(6, 3)  = dNw1;
    N2(6, 4)  = dNw2;
    N2(7, 1)  = ddNv1;
    N2(7, 2)  = ddNv2;
    N2(8, 3)  = ddNw1;
    N2(8, 4)  = ddNw2;
    N2(9, 5)  = Nf1;
    N2(10, 5) = oneOverL;

    // Transformation matrix - transform axial force form centroid to shear center
    Tr(0, 0) = 1.0;
    Tr(1, 1) = 1.0;
    Tr(2, 2) = 1.0;
    Tr(3, 3) = 1.0;
    Tr(4, 4) = 1.0;
    Tr(5, 5) = 1.0;
    Tr(0, 1) = -ys;
    Tr(0, 2) =  ys;
    Tr(0, 3) =  zs;
    Tr(0, 4) = -zs;

    // Get the section tangent stiffness and stress resultant
    const Matrix &ks = theSections[i]->getSectionTangent();
    const Vector &s  = theSections[i]->getStressResultant();

    // Beam material stiffness matrix
    N3.addMatrixTripleProduct(0.0, N1, ks, 1.0);
    kbPart1.addMatrixTripleProduct(0.0, N2, N3, 1.0);

    // Beam geometric stiffness matrix
    Gm( 1,  1) = Gm(2,  2) = Gm(3, 3) = Gm(4, 4) = s(0) * 4.0 / 30.0; // 4/30*N
    Gm( 1,  3) = Gm(2,  4) = Gm(3, 1) = Gm(4, 2) = -s(0) / 30.0;      // -1/30*N
    Gm( 9,  8) = Gm(8,  9) =  s(1);                                        //Mz
    Gm( 9,  7) = Gm(7,  9) =  s(2);                                        //My
    Gm(10,  5) = Gm(5, 10) =  s(0) * zs;                                 //Nzs
    Gm(10,  6) = Gm(6, 10) = -s(0) * ys;                                //-Nys
    Gm(10, 10)               =  s(3);                                      //W

    kbPart2.addMatrixTripleProduct(0.0, N2, Gm, 1.0);

    //perform transformation - transform axial force form centroid to shear center
    kf1.addMatrixTripleProduct(0.0, Tr, kbPart1, 1.0);
    kf2.addMatrixTripleProduct(0.0, Tr, kbPart2, 1.0);

    //perform numerical integration
    double wti = wt[i];
    for (int j = 0; j < 6; j++) {
      for (int k = 0; k < 6; k++) {
        kb(j, k) += kf1(j, k) * L * wti + kf2(j, k) * L * wti;
      }
    }

    //assemble internal force vector q
    static Vector qProduct1(11);
    static Vector qProduct2(6);
    static Vector qProduct3(6);
    qProduct1.Zero();
    qProduct2.Zero();
    qProduct3.Zero();
    qProduct1.addMatrixTransposeVector(0.0, N1, s, 1.0);
    qProduct2.addMatrixTransposeVector(0.0, N2, qProduct1, 1.0);
    qProduct3.addMatrixTransposeVector(0.0, Tr, qProduct2, 1.0);

    for (int j = 0; j < 6; j++) {
      q(j) += qProduct3(j) * L * wti;
    }
  }

  q[0] += q0[0];
  q[1] += q0[1];
  q[2] += q0[2];
  q[3] += q0[3];
  q[4] += q0[4];

  // Transform to global stiffness
  K = crdTransf->getGlobalStiffMatrix(kb, q);
  return K;
}

const Matrix &
DispBeamColumnAsym3d::getInitialBasicStiff()
{

  static Matrix kb(6, 6);
  static Matrix N1(5, 11);  //Xinlong
  static Matrix N2(11, 6);  //Xinlong
  static Matrix N3(11, 11); //Xinlong
  static Matrix kbPart1(6, 6);
  static Matrix Gm(11, 11);
  static Matrix kbPart2(6, 6);
  static Matrix Tr(6, 6);
  static Matrix kf1(6, 6);
  static Matrix kf2(6, 6);

  // Zero for integral
  kb.Zero();

  const Vector &v = crdTransf->getBasicTrialDisp();

  double L        = crdTransf->getInitialLength();
  double oneOverL = 1.0 / L;

  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  //const Vector &wts = quadRule.getIntegrPointWeights(numSections);
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    N1.Zero();
    N2.Zero();
    N3.Zero();
    kbPart1.Zero();
    kbPart2.Zero();
    Gm.Zero();

    Tr.Zero();
    kf1.Zero();
    kf2.Zero();

    double xi1   = xi[i];
    double dNv1  = 1.0 + 3.0 * xi1 * xi1 - 4.0 * xi1;
    double ddNv1 = 6.0 * xi1 * oneOverL - 4.0 * oneOverL;
    double dNv2  = 3.0 * xi1 * xi1 - 2.0 * xi1;
    double ddNv2 = 6.0 * xi1 * oneOverL - 2.0 * oneOverL;
    double dNw1  = -dNv1;
    double ddNw1 = -ddNv1;
    double dNw2  = -dNv2;
    double ddNw2 = -ddNv2;
    double Nf1   = xi1;

    double dv  = dNv1 * v(1) + dNv2 * v(2);   // v'
    double ddv = ddNv1 * v(1) + ddNv2 * v(2); // v"
    double dw  = dNw1 * v(3) + dNw2 * v(4);   // w'
    double ddw = ddNw1 * v(3) + ddNw2 * v(4); // w"
    double f   = Nf1 * v(5);                  // phi
    double df  = oneOverL * v(5);             // phi'

    //Matrix Ndeltad1-------------------------------------------------------
    N1(0,  0) = 1.0;
    N1(0,  1) = (4.0 * v(1) - v(2)) / 30.0;
    N1(0,  2) = (4.0 * v(3) - v(4)) / 30.0;
    N1(0,  3) = (4.0 * v(2) - v(1)) / 30.0;
    N1(0,  4) = (4.0 * v(4) - v(3)) / 30.0;
    N1(0,  5) =  zs * df;
    N1(0,  6) = -ys * df;
    N1(0, 10) = zs * dv - ys * dw;
    N1(1,  7) = 1.0;
    N1(1,  8) = f;
    N1(1,  9) = ddw;
    N1(2,  7) = f;
    N1(2,  8) = -1.0;
    N1(2,  9) = ddv;
    N1(3, 10) = df;
    N1(4, 10) = 1.0;

    //Matrix Ndeltad2-------------------------------------------------------
    N2( 0, 0)  = oneOverL;
    N2( 1, 1)  = 1.0;
    N2( 2, 3)  = 1.0;
    N2( 3, 2)  = 1.0;
    N2( 4, 4)  = 1.0;
    N2( 5, 1)  = dNv1;
    N2( 5, 2)  = dNv2;
    N2( 6, 3)  = dNw1;
    N2( 6, 4)  = dNw2;
    N2( 7, 1)  = ddNv1;
    N2( 7, 2)  = ddNv2;
    N2( 8, 3)  = ddNw1;
    N2( 8, 4)  = ddNw2;
    N2( 9, 5)  = Nf1;
    N2(10, 5)  = oneOverL;

    //Transformation matrix - transform axial force form centroid to shear center
    Tr(0, 0) = 1.0;
    Tr(1, 1) = 1.0;
    Tr(2, 2) = 1.0;
    Tr(3, 3) = 1.0;
    Tr(4, 4) = 1.0;
    Tr(5, 5) = 1.0;
    Tr(0, 1) = -ys;
    Tr(0, 2) =  ys;
    Tr(0, 3) =  zs;
    Tr(0, 4) = -zs;

    // Get the section tangent stiffness and stress resultant
    const Matrix &ks = theSections[i]->getInitialTangent();
    const Vector &s  = theSections[i]->getStressResultant();

    //calculate kbPart1 - material stiffness matrix
    N3.addMatrixTripleProduct(0.0, N1, ks, 1.0);
    kbPart1.addMatrixTripleProduct(0.0, N2, N3, 1.0);

    //calculate kbPart2 - geometric stiffness matrix
    Gm( 1, 1) = Gm(2,  2) = Gm(3, 3) = Gm(4, 4) =  s(0) * 4.0 / 30.0; //  4/30*N
    Gm( 1, 3) = Gm(2,  4) = Gm(3, 1) = Gm(4, 2) = -s(0) / 30.0;       // -1/30*N
    Gm( 9, 8) = Gm(8,  9) =  s(1);                                    //  Mz
    Gm( 9, 7) = Gm(7,  9) =  s(2);                                    //  My
    Gm(10, 5) = Gm(5, 10) =  s(0) * zs;                               //  Nzs
    Gm(10, 6) = Gm(6, 10) = -s(0) * ys;                               // -Nys
    Gm(10, 10)            =  s(3);                                    //  W

    kbPart2.addMatrixTripleProduct(0.0, N2, Gm, 1.0);

    //perform transformation - transform axial force form centroid to shear center
    kf1.addMatrixTripleProduct(0.0, Tr, kbPart1, 1.0);
    kf2.addMatrixTripleProduct(0.0, Tr, kbPart2, 1.0);

    //perform numerical integration
    double wti = wt[i];
    for (int j = 0; j < 6; j++) {
      for (int k = 0; k < 6; k++) {
        kb(j, k) += kf1(j, k) * L * wti + kf2(j, k) * L * wti;
      }
    }
  }

  return kb;
}

const Matrix &
DispBeamColumnAsym3d::getInitialStiff()
{
  const Matrix &kb = this->getInitialBasicStiff();

  // Transform to global stiffness
  K = crdTransf->getInitialGlobalStiffMatrix(kb);

  return K;
}

const Matrix &
DispBeamColumnAsym3d::getMass()
{
  K.Zero();

  if (rho == 0.0)
    return K;

  double L = crdTransf->getInitialLength();
  if (cMass == 0) {
    // lumped mass matrix
    double m = 0.5 * rho * L;
    K(0, 0) = K(1, 1) = K(2, 2) = K(6, 6) = K(7, 7) = K(8, 8) = m;
  } else {
    // consistent mass matrix
    static Matrix ml(12, 12);
    double m = rho * L / 420.0;
    ml(0, 0) = ml(6, 6) = m * 140.0;
    ml(0, 6) = ml(6, 0) = m * 70.0;
    //ml(3,3) = ml(9,9) = m*(Jx/A)*140.0;  // CURRENTLY NO TORSIONAL MASS
    //ml(3,9) = ml(9,3) = m*(Jx/A)*70.0;   // CURRENTLY NO TORSIONAL MASS

    ml(2, 2) = ml(8, 8) = m * 156.0;
    ml(2, 8) = ml(8, 2) = m * 54.0;
    ml(4, 4) = ml(10, 10) = m * 4.0 * L * L;
    ml(4, 10) = ml(10, 4) = -m * 3.0 * L * L;
    ml(2, 4) = ml(4, 2) = -m * 22.0 * L;
    ml(8, 10) = ml(10, 8) = -ml(2, 4);
    ml(2, 10) = ml(10, 2) = m * 13.0 * L;
    ml(4, 8) = ml(8, 4) = -ml(2, 10);

    ml(1, 1) = ml(7, 7) = m * 156.0;
    ml(1, 7) = ml(7, 1) = m * 54.0;
    ml(5, 5) = ml(11, 11) = m * 4.0 * L * L;
    ml(5, 11) = ml(11, 5) = -m * 3.0 * L * L;
    ml(1, 5) = ml(5, 1) = m * 22.0 * L;
    ml(7, 11) = ml(11, 7) = -ml(1, 5);
    ml(1, 11) = ml(11, 1) = -m * 13.0 * L;
    ml(5, 7) = ml(7, 5) = -ml(1, 11);

    // transform local mass matrix to global system
    K = crdTransf->getGlobalMatrixFromLocal(ml);
  }

  return K;
}

void
DispBeamColumnAsym3d::zeroLoad()
{
  Q.Zero();

  q0.zero();
  p0.zero();

  return;
}

int
DispBeamColumnAsym3d::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);
  double L           = crdTransf->getInitialLength();

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
  } else {
    opserr
        << "DispBeamColumnAsym3d::addLoad() -- load type unknown for element with tag: "
        << this->getTag() << endln;
    return -1;
  }

  return 0;
}

int
DispBeamColumnAsym3d::addInertiaLoadToUnbalance(const Vector &accel)
{
  // Check for a quick return
  if (rho == 0.0)
    return 0;

  // Get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);

  if (6 != Raccel1.Size() || 6 != Raccel2.Size()) {
    opserr << "DispBeamColumnAsym3d::addInertiaLoadToUnbalance matrix and vector sizes "
              "are incompatible\n";
    return -1;
  }

  // want to add ( - fact * M R * accel ) to unbalance
  if (cMass == 0) {
    // take advantage of lumped mass matrix
    double L = crdTransf->getInitialLength();
    double m = 0.5 * rho * L;

    Q(0) -= m * Raccel1(0);
    Q(1) -= m * Raccel1(1);
    Q(2) -= m * Raccel1(2);
    Q(6) -= m * Raccel2(0);
    Q(7) -= m * Raccel2(1);
    Q(8) -= m * Raccel2(2);

  } else {
    // use matrix vector multip. for consistent mass matrix
    static Vector Raccel(12);
    for (int i = 0; i < 6; i++) {
      Raccel(i)     = Raccel1(i);
      Raccel(i + 6) = Raccel2(i);
    }
    Q.addMatrixVector(1.0, this->getMass(), Raccel, -1.0);
  }

  return 0;
}

const Vector &
DispBeamColumnAsym3d::getResistingForce()
{

  static Matrix N1(5, 11); //Xinlong
  static Matrix N2(11, 6); //Xinlong
  static Matrix Tr(6, 6);

  const Vector &v = crdTransf->getBasicTrialDisp();

  // Zero for integral
  q.Zero();

  double L        = crdTransf->getInitialLength();
  double oneOverL = 1.0 / L;

  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  //const Vector &wts = quadRule.getIntegrPointWeights(numSections);
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    N1.Zero();
    N2.Zero();
    Tr.Zero();

    double xi1   = xi[i];
    double dNv1  = 1.0 + 3.0 * xi1 * xi1 - 4.0 * xi1;
    double ddNv1 = 6.0 * xi1 * oneOverL - 4.0 * oneOverL;
    double dNv2  = 3.0 * xi1 * xi1 - 2.0 * xi1;
    double ddNv2 = 6.0 * xi1 * oneOverL - 2.0 * oneOverL;
    double dNw1  = -dNv1;
    double ddNw1 = -ddNv1;
    double dNw2  = -dNv2;
    double ddNw2 = -ddNv2;
    double Nf1   = xi1;

    double dv  = dNv1 * v(1) + dNv2 * v(2);   //v'
    double ddv = ddNv1 * v(1) + ddNv2 * v(2); //v"
    double dw  = dNw1 * v(3) + dNw2 * v(4);   //w'
    double ddw = ddNw1 * v(3) + ddNw2 * v(4); //w"
    double f   = Nf1 * v(5);                  //phi
    double df  = oneOverL * v(5);             //phi'

    //Matrix Ndeltad1-------------------------------------------------------
    N1(0, 0)  = 1.0;
    N1(0, 1)  = (4.0 * v(1) - v(2)) / 30.0;
    N1(0, 2)  = (4.0 * v(3) - v(4)) / 30.0;
    N1(0, 3)  = (4.0 * v(2) - v(1)) / 30.0;
    N1(0, 4)  = (4.0 * v(4) - v(3)) / 30.0;
    N1(0, 5)  = zs * df;
    N1(0, 6)  = -ys * df;
    N1(0, 10) = zs * dv - ys * dw;
    N1(1, 7)  = 1.0;
    N1(1, 8)  = f;
    N1(1, 9)  = ddw;
    N1(2, 7)  = f;
    N1(2, 8)  = -1.0;
    N1(2, 9)  = ddv;
    N1(3, 10) = df;
    N1(4, 10) = 1.0;

    //Matrix Ndeltad2-------------------------------------------------------
    N2(0, 0)  = oneOverL;
    N2(1, 1)  = 1.0;
    N2(2, 3)  = 1.0;
    N2(3, 2)  = 1.0;
    N2(4, 4)  = 1.0;
    N2(5, 1)  = dNv1;
    N2(5, 2)  = dNv2;
    N2(6, 3)  = dNw1;
    N2(6, 4)  = dNw2;
    N2(7, 1)  = ddNv1;
    N2(7, 2)  = ddNv2;
    N2(8, 3)  = ddNw1;
    N2(8, 4)  = ddNw2;
    N2(9, 5)  = Nf1;
    N2(10, 5) = oneOverL;

    //Transformation matrix - transform axial force form centroid to shear center
    Tr(0, 0) = 1.0;
    Tr(1, 1) = 1.0;
    Tr(2, 2) = 1.0;
    Tr(3, 3) = 1.0;
    Tr(4, 4) = 1.0;
    Tr(5, 5) = 1.0;
    Tr(0, 1) = -ys;
    Tr(0, 2) = ys;
    Tr(0, 3) = zs;
    Tr(0, 4) = -zs;

    // Get the section stress resultant
    const Vector &s = theSections[i]->getStressResultant();

    //perform numerical integration
    double wti = wt[i];

    //assemble internal force vector q
    static Vector qProduct1(11);
    static Vector qProduct2(6);
    static Vector qProduct3(6);
    qProduct1.Zero();
    qProduct2.Zero();
    qProduct3.Zero();
    qProduct1.addMatrixTransposeVector(0.0, N1, s, 1.0);
    qProduct2.addMatrixTransposeVector(0.0, N2, qProduct1, 1.0);
    qProduct3.addMatrixTransposeVector(0.0, Tr, qProduct2, 1.0);

    for (int j = 0; j < 6; j++) {
      q(j) += qProduct3(j) * L * wti;
    }
  }

  q(0) += q0[0];
  q(1) += q0[1];
  q(2) += q0[2];
  q(3) += q0[3];
  q(4) += q0[4];

  // Transform forces
  Vector p0Vec(p0);
  P = crdTransf->getGlobalResistingForce(q, p0Vec);

  // Subtract other external nodal loads ... P_res = P_int - P_ext
  if (rho != 0)
    P.addVector(1.0, Q, -1.0);

  return P;
}

const Vector &
DispBeamColumnAsym3d::getResistingForceIncInertia()
{
  P = this->getResistingForce();

  if (rho != 0.0) {
    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();

    if (cMass == 0) {
      // take advantage of lumped mass matrix
      double L = crdTransf->getInitialLength();
      double m = 0.5 * rho * L;

      P(0) += m * accel1(0);
      P(1) += m * accel1(1);
      P(2) += m * accel1(2);
      P(6) += m * accel2(0);
      P(7) += m * accel2(1);
      P(8) += m * accel2(2);
    } else {
      // use matrix vector multip. for consistent mass matrix
      static Vector accel(12);
      for (int i = 0; i < 6; i++) {
        accel(i)     = accel1(i);
        accel(i + 6) = accel2(i);
      }
      P.addMatrixVector(1.0, this->getMass(), accel, 1.0);
    }

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
DispBeamColumnAsym3d::sendSelf(int commitTag, Channel &theChannel)
{
  // place the integer data into an ID

  int dbTag = this->getDbTag();
  int i, j;
  int loc = 0;

  static Vector data(16);
  data(0)            = this->getTag();
  data(1)            = connectedExternalNodes(0);
  data(2)            = connectedExternalNodes(1);
  data(3)            = numSections;
  data(4)            = crdTransf->getClassTag();
  int crdTransfDbTag = crdTransf->getDbTag();
  if (crdTransfDbTag == 0) {
    crdTransfDbTag = theChannel.getDbTag();
    if (crdTransfDbTag != 0)
      crdTransf->setDbTag(crdTransfDbTag);
  }
  data(5)          = crdTransfDbTag;
  data(6)          = beamInt->getClassTag();
  int beamIntDbTag = beamInt->getDbTag();
  if (beamIntDbTag == 0) {
    beamIntDbTag = theChannel.getDbTag();
    if (beamIntDbTag != 0)
      beamInt->setDbTag(beamIntDbTag);
  }
  data(7)  = beamIntDbTag;
  data(8)  = rho;
  data(9)  = cMass;
  data(10) = alphaM;
  data(11) = betaK;
  data(12) = betaK0;
  data(13) = betaKc;
  data(14) = ys;
  data(15) = zs;

  if (theChannel.sendVector(dbTag, commitTag, data) < 0) {
    opserr << "DispBeamColumnAsym3d::sendSelf() - failed to send data Vector\n";
    return -1;
  }

  // send the coordinate transformation
  if (crdTransf->sendSelf(commitTag, theChannel) < 0) {
    opserr << "DispBeamColumnAsym3d::sendSelf() - failed to send crdTranf\n";
    return -1;
  }

  // send the beam integration
  if (beamInt->sendSelf(commitTag, theChannel) < 0) {
    opserr << "DispBeamColumnAsym3d::sendSelf() - failed to send beamInt\n";
    return -1;
  }

  //
  // send an ID for the sections containing each sections dbTag and classTag
  // if section ha no dbTag get one and assign it
  //

  ID idSections(2 * numSections);
  loc = 0;
  for (i = 0; i < numSections; i++) {
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
    opserr << "DispBeamColumnAsym3d::sendSelf() - failed to send ID data\n";
    return -1;
  }

  //
  // send the sections
  //

  for (int j = 0; j < numSections; j++) {
    if (theSections[j]->sendSelf(commitTag, theChannel) < 0) {
      opserr << "DispBeamColumnAsym3d::sendSelf() - section " << j
             << "failed to send itself\n";
      return -1;
    }
  }

  return 0;
}

int
DispBeamColumnAsym3d::recvSelf(int commitTag, Channel &theChannel,
                               FEM_ObjectBroker &theBroker)
{
  //
  // receive the integer data containing tag, numSections and coord transformation info
  //
  int dbTag = this->getDbTag();
  int i;

  static Vector data(16);

  if (theChannel.recvVector(dbTag, commitTag, data) < 0) {
    opserr << "DispBeamColumnAsym3d::recvSelf() - failed to recv data Vector\n";
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

  rho   = data(8);
  cMass = (int)data(9);

  alphaM = data(10);
  betaK  = data(11);
  betaK0 = data(12);
  betaKc = data(13);
  ys     = data(14);
  zs     = data(15);

  // create a new crdTransf object if one needed
  if (crdTransf == 0 || crdTransf->getClassTag() != crdTransfClassTag) {
    if (crdTransf != 0)
      delete crdTransf;

    crdTransf = theBroker.getNewCrdTransf(crdTransfClassTag);

    if (crdTransf == 0) {
      opserr << "DispBeamColumnAsym3d::recvSelf() - "
             << "failed to obtain a CrdTrans object with classTag" << crdTransfClassTag
             << endln;
      return -2;
    }
  }

  crdTransf->setDbTag(crdTransfDbTag);

  // invoke recvSelf on the crdTransf object
  if (crdTransf->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "DispBeamColumnAsym3d::sendSelf() - failed to recv crdTranf\n";
    return -3;
  }

  // create a new beamInt object if one needed
  if (beamInt == 0 || beamInt->getClassTag() != beamIntClassTag) {
    if (beamInt != 0)
      delete beamInt;

    beamInt = theBroker.getNewBeamIntegration(beamIntClassTag);

    if (beamInt == 0) {
      opserr << "DispBeamColumnAsym3d::recvSelf() - failed to obtain the beam "
                "integration object with classTag"
             << beamIntClassTag << endln;
      exit(-1);
    }
  }

  beamInt->setDbTag(beamIntDbTag);

  // invoke recvSelf on the beamInt object
  if (beamInt->recvSelf(commitTag, theChannel, theBroker) < 0) {
    opserr << "DispBeamColumnAsym3d::sendSelf() - failed to recv beam integration\n";
    return -3;
  }

  //
  // recv an ID for the sections containing each sections dbTag and classTag
  //

  ID idSections(2 * nSect);
  int loc = 0;

  if (theChannel.recvID(dbTag, commitTag, idSections) < 0) {
    opserr << "DispBeamColumnAsym3d::recvSelf() - failed to recv ID data\n";
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
    theSections = new SectionForceDeformation *[nSect];
    if (theSections == 0) {
      opserr << "DispBeamColumnAsym3d::recvSelf() - out of memory creating sections "
                "array of size"
             << nSect << endln;
      exit(-1);
    }

    // create a section and recvSelf on it
    numSections = nSect;
    loc         = 0;

    for (i = 0; i < numSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag    = idSections(loc + 1);
      loc += 2;
      theSections[i] = theBroker.getNewSection(sectClassTag);
      if (theSections[i] == 0) {
        opserr << "DispBeamColumnAsym3d::recvSelf() - Broker could not create Section of "
                  "class type"
               << sectClassTag << endln;
        exit(-1);
      }
      theSections[i]->setDbTag(sectDbTag);
      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
        opserr << "DispBeamColumnAsym3d::recvSelf() - section " << i
               << "failed to recv itself\n";
        return -1;
      }
    }

  } else {

    //
    // for each existing section, check it is of correct type
    // (if not delete old & create a new one) then recvSelf on it
    //

    loc = 0;
    for (i = 0; i < numSections; i++) {
      int sectClassTag = idSections(loc);
      int sectDbTag    = idSections(loc + 1);
      loc += 2;

      // check of correct type
      if (theSections[i]->getClassTag() != sectClassTag) {
        // delete the old section[i] and create a new one
        delete theSections[i];
        theSections[i] = theBroker.getNewSection(sectClassTag);
        if (theSections[i] == 0) {
          opserr << "DispBeamColumnAsym3d::recvSelf() - Broker could not create Section "
                    "of class type"
                 << sectClassTag << endln;
          exit(-1);
        }
      }

      // recvSelf on it
      theSections[i]->setDbTag(sectDbTag);
      if (theSections[i]->recvSelf(commitTag, theChannel, theBroker) < 0) {
        opserr << "DispBeamColumnAsym3d::recvSelf() - section " << i
               << "failed to recv itself\n";
        return -1;
      }
    }
  }

  return 0;
}

void
DispBeamColumnAsym3d::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "\nDispBeamColumnAsym3d, element id:  " << this->getTag() << endln;
    s << "\tConnected external nodes:  " << connectedExternalNodes;
    s << "\tCoordTransf: " << crdTransf->getTag() << endln;
    s << "\tmass density:  " << rho << ", cMass: " << cMass << endln;

    double N, Mz1, Mz2, Vy, My1, My2, Vz, T;
    double L        = crdTransf->getInitialLength();
    double oneOverL = 1.0 / L;

    N   = q(0);
    Mz1 = q(1);
    Mz2 = q(2);
    Vy  = (Mz1 + Mz2) * oneOverL;
    My1 = q(3);
    My2 = q(4);
    Vz  = -(My1 + My2) * oneOverL;
    T   = q(5);

    s << "\tEnd 1 Forces (P Mz Vy My Vz T): " << -N + p0[0] << ' ' << Mz1 << ' '
      << Vy + p0[1] << ' ' << My1 << ' ' << Vz + p0[3] << ' ' << -T << endln;
    s << "\tEnd 2 Forces (P Mz Vy My Vz T): " << N << ' ' << Mz2 << ' ' << -Vy + p0[2]
      << ' ' << My2 << ' ' << -Vz + p0[4] << ' ' << T << endln;

    beamInt->Print(s, flag);

    for (int i = 0; i < numSections; i++) {
      //opserr << "Section Type: " << theSections[i]->getClassTag() << endln;
      theSections[i]->Print(s, flag);
    }
    //  if (rho != 0)
    //    opserr << "Mass: \n" << this->getMass();
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"DispBeamColumnAsym3d\", ";
    s << "\"nodes\": [" << connectedExternalNodes(0) << ", " << connectedExternalNodes(1)
      << "], ";
    s << "\"sections\": [";
    for (int i = 0; i < numSections - 1; i++)
      s << "\"" << theSections[i]->getTag() << "\", ";
    s << "\"" << theSections[numSections - 1]->getTag() << "\"], ";
    s << "\"integration\": ";
    beamInt->Print(s, flag);
    s << ", \"massperlength\": " << rho << ", ";
    s << "\"crdTransformation\": \"" << crdTransf->getTag() << "\"}";
  }
}

int
DispBeamColumnAsym3d::displaySelf(Renderer &theViewer, int displayMode, float fact,
                                  const char **modes, int numModes)
{
  static Vector v1(3);
  static Vector v2(3);

  if (displayMode >= 0) {

    theNodes[0]->getDisplayCrds(v1, fact);
    theNodes[1]->getDisplayCrds(v2, fact);

  } else {

    theNodes[0]->getDisplayCrds(v1, 0.);
    theNodes[1]->getDisplayCrds(v2, 0.);

    // add eigenvector values
    int mode             = displayMode * -1;
    const Matrix &eigen1 = theNodes[0]->getEigenvectors();
    const Matrix &eigen2 = theNodes[1]->getEigenvectors();
    if (eigen1.noCols() >= mode) {
      for (int i = 0; i < 3; i++) {
        v1(i) += eigen1(i, mode - 1) * fact;
        v2(i) += eigen2(i, mode - 1) * fact;
      }
    }
  }
  return theViewer.drawLine(v1, v2, 1.0, 1.0, this->getTag());
}

Response *
DispBeamColumnAsym3d::setResponse(const char **argv, int argc, OPS_Stream &output)
{

  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType", "DispBeamColumnAsym3d");
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

    theResponse = new ElementResponse(this, 1, P);

    // local force -
  } else if (strcmp(argv[0], "localForce") == 0 || strcmp(argv[0], "localForces") == 0) {

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

    theResponse = new ElementResponse(this, 2, P);

    // chord rotation -
  } else if (strcmp(argv[0], "chordRotation") == 0 ||
             strcmp(argv[0], "chordDeformation") == 0 ||
             strcmp(argv[0], "basicDeformation") == 0) {

    output.tag("ResponseType", "eps");
    output.tag("ResponseType", "thetaZ_1");
    output.tag("ResponseType", "thetaZ_2");
    output.tag("ResponseType", "thetaY_1");
    output.tag("ResponseType", "thetaY_2");
    output.tag("ResponseType", "thetaX");

    theResponse = new ElementResponse(this, 3, Vector(6));

    // plastic rotation -
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

    theResponse = new ElementResponse(this, 12, P);

  } else if (strcmp(argv[0], "integrationPoints") == 0)
    theResponse = new ElementResponse(this, 10, Vector(numSections));

  else if (strcmp(argv[0], "integrationWeights") == 0)
    theResponse = new ElementResponse(this, 11, Vector(numSections));

  else if (strcmp(argv[0], "sectionTags") == 0)
    theResponse = new ElementResponse(this, 110, ID(numSections));

  // section response -
  else if (strstr(argv[0], "sectionX") != 0) {
    if (argc > 2) {
      float sectionLoc = atof(argv[1]);

      double xi[maxNumSections];
      double L = crdTransf->getInitialLength();
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
        double L = crdTransf->getInitialLength();
        beamInt->getSectionLocations(numSections, L, xi);

        output.tag("GaussPointOutput");
        output.attr("number", sectionNum);
        output.attr("eta", xi[sectionNum - 1] * L);

        theResponse =
            theSections[sectionNum - 1]->setResponse(&argv[2], argc - 2, output);

        output.endTag();
      } else if (sectionNum == 0) { // argv[1] was not an int, we want all sections,

        CompositeResponse *theCResponse = new CompositeResponse();
        int numResponse                 = 0;
        double xi[maxNumSections];
        double L = crdTransf->getInitialLength();
        beamInt->getSectionLocations(numSections, L, xi);

        for (int i = 0; i < numSections; i++) {

          output.tag("GaussPointOutput");
          output.attr("number", i + 1);
          output.attr("eta", xi[i] * L);

          Response *theSectionResponse =
              theSections[i]->setResponse(&argv[1], argc - 1, output);

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

  output.endTag();
  return theResponse;
}

int
DispBeamColumnAsym3d::getResponse(int responseID, Information &eleInfo)
{
  double N, V, M1, M2, T;
  double L        = crdTransf->getInitialLength();
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

  // Chord rotation
  else if (responseID == 3)
    return eleInfo.setVector(crdTransf->getBasicTrialDisp());

  // Plastic rotation
  else if (responseID == 4) {
    static Vector vp(6);
    static Vector ve(6);
    const Matrix &kb = this->getInitialBasicStiff();
    kb.Solve(q, ve);
    vp = crdTransf->getBasicTrialDisp();
    vp -= ve;
    return eleInfo.setVector(vp);
  }

  else if (responseID == 10) {
    double L = crdTransf->getInitialLength();
    double pts[maxNumSections];
    beamInt->getSectionLocations(numSections, L, pts);
    Vector locs(numSections);
    for (int i = 0; i < numSections; i++)
      locs(i) = pts[i] * L;
    return eleInfo.setVector(locs);
  }

  else if (responseID == 11) {
    double L = crdTransf->getInitialLength();
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

  else
    return -1;
}

// AddingSensitivity:BEGIN ///////////////////////////////////
int
DispBeamColumnAsym3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  // don't do anything if MaterialStageParameter calls this element
  if (strcmp(argv[0], "updateMaterialStage") == 0) {
    return -1;
  }

  // If the parameter belongs to the element itself
  if (strcmp(argv[0], "rho") == 0) {
    param.setValue(rho);
    return param.addObject(1, this);
  }

  if (strstr(argv[0], "sectionX") != 0) {
    if (argc < 3)
      return -1;

    float sectionLoc = atof(argv[1]);

    double xi[maxNumSections];
    double L = crdTransf->getInitialLength();
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

  else if (strstr(argv[0], "integration") != 0) {

    if (argc < 2)
      return -1;

    return beamInt->setParameter(&argv[1], argc - 1, param);
  }

  // Default, send to every object
  int ok     = 0;
  int result = 0;

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
DispBeamColumnAsym3d::updateParameter(int parameterID, Information &info)
{
  if (parameterID == 1) {
    rho = info.theDouble;
    return 0;
  } else
    return -1;
}

int
DispBeamColumnAsym3d::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;

  return 0;
}

const Matrix &
DispBeamColumnAsym3d::getKiSensitivity(int gradNumber)
{
  K.Zero();
  return K;
}

const Matrix &
DispBeamColumnAsym3d::getMassSensitivity(int gradNumber)
{
  K.Zero();

  if (rho == 0.0 || parameterID != 1)
    return K;

  double L = crdTransf->getInitialLength();
  if (cMass == 0) {
    // lumped mass matrix
    //double m = 0.5*rho*L;
    double m = 0.5 * L;
    K(0, 0) = K(1, 1) = K(2, 2) = K(6, 6) = K(7, 7) = K(8, 8) = m;
  } else {
    // consistent mass matrix
    static Matrix ml(12, 12);
    //double m = rho*L/420.0;
    double m = L / 420.0;
    ml(0, 0) = ml(6, 6) = m * 140.0;
    ml(0, 6) = ml(6, 0) = m * 70.0;
    //ml(3,3) = ml(9,9) = m*(Jx/A)*140.0;  // CURRENTLY NO TORSIONAL MASS
    //ml(3,9) = ml(9,3) = m*(Jx/A)*70.0;   // CURRENTLY NO TORSIONAL MASS

    ml(2, 2) = ml(8, 8) = m * 156.0;
    ml(2, 8) = ml(8, 2) = m * 54.0;
    ml(4, 4) = ml(10, 10) = m * 4.0 * L * L;
    ml(4, 10) = ml(10, 4) = -m * 3.0 * L * L;
    ml(2, 4) = ml(4, 2) = -m * 22.0 * L;
    ml(8, 10) = ml(10, 8) = -ml(2, 4);
    ml(2, 10) = ml(10, 2) = m * 13.0 * L;
    ml(4, 8) = ml(8, 4) = -ml(2, 10);

    ml(1, 1) = ml(7, 7) = m * 156.0;
    ml(1, 7) = ml(7, 1) = m * 54.0;
    ml(5, 5) = ml(11, 11) = m * 4.0 * L * L;
    ml(5, 11) = ml(11, 5) = -m * 3.0 * L * L;
    ml(1, 5) = ml(5, 1) = m * 22.0 * L;
    ml(7, 11) = ml(11, 7) = -ml(1, 5);
    ml(1, 11) = ml(11, 1) = -m * 13.0 * L;
    ml(5, 7) = ml(7, 5) = -ml(1, 11);

    // transform local mass matrix to global system
    K = crdTransf->getGlobalMatrixFromLocal(ml);
  }

  return K;
}

const Vector &
DispBeamColumnAsym3d::getResistingForceSensitivity(int gradNumber)
{
  double L        = crdTransf->getInitialLength();
  double oneOverL = 1.0 / L;

  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  //const Vector &wts = quadRule.getIntegrPointWeights(numSections);
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);
  double wt[maxNumSections];
  beamInt->getSectionWeights(numSections, L, wt);

  // Zero for integration
  static Vector dqdh(6);
  dqdh.Zero();

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    int order      = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();

    //double xi6 = 6.0*pts(i,0);
    double xi6 = 6.0 * xi[i];
    //double wti = wts(i);
    double wti = wt[i];

    // Get section stress resultant gradient
    const Vector &dsdh = theSections[i]->getStressResultantSensitivity(gradNumber, true);

    // Perform numerical integration on internal force gradient
    double sensi;
    for (int j = 0; j < order; j++) {
      sensi = dsdh(j) * wti;
      switch (code(j)) {
      case SECTION_RESPONSE_P:
        dqdh(0) += sensi;
        break;
      case SECTION_RESPONSE_MZ:
        dqdh(1) += (xi6 - 4.0) * sensi;
        dqdh(2) += (xi6 - 2.0) * sensi;
        break;
      case SECTION_RESPONSE_MY:
        dqdh(3) += (xi6 - 4.0) * sensi;
        dqdh(4) += (xi6 - 2.0) * sensi;
        break;
      case SECTION_RESPONSE_T:
        dqdh(5) += sensi;
        break;
      default:
        break;
      }
    }
  }

  // Transform forces
  static Vector dp0dh(6); // No distributed loads

  P.Zero();

  //////////////////////////////////////////////////////////////

  if (crdTransf->isShapeSensitivity()) {

    // Perform numerical integration to obtain basic stiffness matrix
    // Some extra declarations
    static Matrix kbmine(6, 6);
    kbmine.Zero();
    q.Zero();

    double tmp;

    int j, k;

    for (int i= 0; i < numSections; i++) {

      int order      = theSections[i]->getOrder();
      const ID &code = theSections[i]->getType();

      //double xi6 = 6.0*pts(i,0);
      double xi6 = 6.0 * xi[i];
      //double wti = wts(i);
      double wti = wt[i];

      const Vector &s  = theSections[i]->getStressResultant();
      const Matrix &ks = theSections[i]->getSectionTangent();

      Matrix ka(workArea, order, 6);
      ka.Zero();

      for (int j = 0; j < order; j++) {
        double si = s(j) * wti;
        switch (code(j)) {
        case SECTION_RESPONSE_P:
          q(0) += si;
          for (k = 0; k < order; k++)
            ka(k, 0) += ks(k, j) * wti;

          break;
        case SECTION_RESPONSE_MZ:
          q(1) += (xi6 - 4.0) * si;
          q(2) += (xi6 - 2.0) * si;
          for (k = 0; k < order; k++) {
            double tmp = ks(k, j) * wti;
            ka(k, 1) += (xi6 - 4.0) * tmp;
            ka(k, 2) += (xi6 - 2.0) * tmp;
          }
          break;
        case SECTION_RESPONSE_MY:
          q(3) += (xi6 - 4.0) * si;
          q(4) += (xi6 - 2.0) * si;
          for (k = 0; k < order; k++) {
            double tmp = ks(k, j) * wti;
            ka(k, 3) += (xi6 - 4.0) * tmp;
            ka(k, 4) += (xi6 - 2.0) * tmp;
          }
          break;
        case SECTION_RESPONSE_T:
          q(5) += si;
          for (k = 0; k < order; k++) {
            ka(k, 5) += ks(k, j) * wti;
          }
          break;
        default:
          break;
        }
      }
      for (int j = 0; j < order; j++) {
        switch (code(j)) {
        case SECTION_RESPONSE_P:
          for (k = 0; k < 6; k++) {
            kbmine(0, k) += ka(j, k);
          }
          break;
        case SECTION_RESPONSE_MZ:
          for (k = 0; k < 6; k++) {
            tmp = ka(j, k);
            kbmine(1, k) += (xi6 - 4.0) * tmp;
            kbmine(2, k) += (xi6 - 2.0) * tmp;
          }
          break;
        case SECTION_RESPONSE_MY:
          for (k = 0; k < 6; k++) {
            tmp = ka(j, k);
            kbmine(3, k) += (xi6 - 4.0) * tmp;
            kbmine(4, k) += (xi6 - 2.0) * tmp;
          }
          break;
        case SECTION_RESPONSE_T:
          for (k = 0; k < 6; k++) {
            kbmine(5, k) += ka(j, k);
          }
          break;
        default:
          break;
        }
      }
    }

    const Vector &A_u = crdTransf->getBasicTrialDisp();
    double dLdh       = crdTransf->getLengthGrad();
    double d1overLdh  = -dLdh / (L * L);
    // a^T k_s dadh v
    dqdh.addMatrixVector(1.0, kbmine, A_u, d1overLdh);

    // k dAdh u
    const Vector &dAdh_u = crdTransf->getBasicDisplFixedGrad();
    dqdh.addMatrixVector(1.0, kbmine, dAdh_u, oneOverL);

    // dAdh^T q
    P += crdTransf->getGlobalResistingForceShapeSensitivity(q, dp0dh, gradNumber);
  }

  // A^T (dqdh + k dAdh u)
  P += crdTransf->getGlobalResistingForce(dqdh, dp0dh);

  return P;
}

// NEW METHOD
int
DispBeamColumnAsym3d::commitSensitivity(int gradNumber, int numGrads)
{
  // Get basic deformation and sensitivities
  const Vector &v = crdTransf->getBasicTrialDisp();

  static Vector dvdh(6);
  dvdh = crdTransf->getBasicDisplTotalGrad(gradNumber);

  double L        = crdTransf->getInitialLength();
  double oneOverL = 1.0 / L;
  //const Matrix &pts = quadRule.getIntegrPointCoords(numSections);
  double xi[maxNumSections];
  beamInt->getSectionLocations(numSections, L, xi);

  // Some extra declarations
  double d1oLdh = crdTransf->getd1overLdh();

  // Loop over the integration points
  for (int i = 0; i < numSections; i++) {

    int order      = theSections[i]->getOrder();
    const ID &code = theSections[i]->getType();

    Vector e(workArea, order);

    double xi6 = 6.0 * xi[i];

    for (int j = 0; j < order; j++) {
      switch (code(j)) {
      case SECTION_RESPONSE_P:
        e(j) = oneOverL * dvdh(0) + d1oLdh * v(0);
        break;
      case SECTION_RESPONSE_MZ:
        e(j) = oneOverL * ((xi6 - 4.0) * dvdh(1) + (xi6 - 2.0) * dvdh(2)) +
               d1oLdh * ((xi6 - 4.0) * v(1) + (xi6 - 2.0) * v(2));
        break;
      case SECTION_RESPONSE_MY:
        e(j) = oneOverL * ((xi6 - 4.0) * dvdh(3) + (xi6 - 2.0) * dvdh(4)) +
               d1oLdh * ((xi6 - 4.0) * v(3) + (xi6 - 2.0) * v(4));
        break;
      case SECTION_RESPONSE_T:
        e(j) = oneOverL * dvdh(5) + d1oLdh * v(5);
        break;
      default:
        e(j) = 0.0;
        break;
      }
    }

    // Set the section deformations
    theSections[i]->commitSensitivity(e, gradNumber, numGrads);
  }

  return 0;
}

// AddingSensitivity:END /////////////////////////////////////////////
