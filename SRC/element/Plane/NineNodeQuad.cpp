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
// based on FourNodeQuad by MHS
// Written: Seweryn Kokot, Opole University of Technology, Poland
// Created: Aug 2020
//
// Description: This file contains the class definition for NineNodeQuad.
//
#include "NineNodeQuad.h"
#include <Node.h>
#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <VectorND.h>
#include <ID.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Parameter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>

using OpenSees::VectorND;

double NineNodeQuad::matrixData[(NEN*2)*(NEN*2)];
Matrix NineNodeQuad::K(matrixData, 2*NEN, 2*NEN);
Vector NineNodeQuad::P(2*NEN);
double NineNodeQuad::shp[3][NEN];

NineNodeQuad::NineNodeQuad(int tag, 
                           const std::array<int,9>& nodes,
                           NDMaterial &m,
                           double thickness,
                           double p, 
                           double rho, 
                           double b1, double b2)

:Element (tag, ELE_TAG_NineNodeQuad),
  theMaterial(0), connectedExternalNodes(NEN),
  Q(2*NEN), applyLoad(0), pressureLoad(2*NEN), 
  thickness(thickness), pressure(p), rho(rho), Ki(0)
{
    // Body forces
    b[0] = b1;
    b[1] = b2;

    // Allocate arrays of pointers to NDMaterials
    theMaterial = new NDMaterial *[nip];

    // Get copies of the material model for each integration point
    for (int i = 0; i < nip; i++) {
      theMaterial[i] = m.getCopy();
    }

    // Set connected external node IDs
    for (int i=0; i<NEN; i++) {
      connectedExternalNodes(i) = nodes[i];
      theNodes[i] = nullptr;
    }
}

NineNodeQuad::NineNodeQuad()
:Element (0,ELE_TAG_NineNodeQuad),
  theMaterial(0), connectedExternalNodes(NEN),
 Q(2*NEN), applyLoad(0), pressureLoad(2*NEN), 
 thickness(0.0), pressure(0.0), Ki(0)
{
  for (int i=0; i<NEN; i++)
    theNodes[i] = nullptr;
}

NineNodeQuad::~NineNodeQuad()
{
  for (int i = 0; i < nip; i++) {
    if (theMaterial[i])
      delete theMaterial[i];
  }

  // Delete the array of pointers to NDMaterial pointer arrays
  if (theMaterial)
    delete [] theMaterial;

  if (Ki != nullptr)
    delete Ki;
}

int
NineNodeQuad::getNumExternalNodes() const
{
  return NEN;
}

const ID&
NineNodeQuad::getExternalNodes()
{
  return connectedExternalNodes;
}


Node **
NineNodeQuad::getNodePtrs()
{
  return &theNodes[0];
}

int
NineNodeQuad::getNumDOF()
{
  return NEN*2;
}

void
NineNodeQuad::setDomain(Domain *theDomain)
{
    // Check Domain is not null - invoked when object removed from a domain
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
      if (dofs != 2 && dofs != 3) {
        opserr << "WARNING " << this->getClassType() 
               << " element with tag " << this->getTag() 
               << " does not have 2 or 3 DOFs at node " 
               << theNodes[i]->getTag() << "\n";
        return;
      }
    }

    this->DomainComponent::setDomain(theDomain);

    // Compute consistent nodal loads due to pressure
    this->setPressureLoadAtNodes();
}

int
NineNodeQuad::commitState()
{
    int retVal = 0;

    // call element commitState to do any base class stuff
    if ((retVal = this->Element::commitState()) != 0) {
      opserr << "NineNodeQuad::commitState () - failed in base class";
    }

    // Loop over the integration points and commit the material states
    for (int i = 0; i < nip; i++)
      retVal += theMaterial[i]->commitState();

    return retVal;
}

int
NineNodeQuad::revertToLastCommit()
{
    int retVal = 0;

    // Loop over the integration points and revert to last committed state
    for (int i = 0; i < nip; i++)
        retVal += theMaterial[i]->revertToLastCommit();

    return retVal;
}

int
NineNodeQuad::revertToStart()
{
    int retVal = 0;

    // Loop over the integration points and revert states to start
    for (int i = 0; i < nip; i++)
        retVal += theMaterial[i]->revertToStart();

    return retVal;
}


int
NineNodeQuad::update()
{
  // Collect displacements at each node into a local array
  double u[NDM][NEN];

  for (int i=0; i<NEN; i++) {
      const Vector &displ = theNodes[i]->getTrialDisp();
      for (int j=0; j<NDM; j++) {
          u[j][i] = displ[j];
      }
  }

  int ret = 0;

  // Loop over the integration points
  for (int i = 0; i < nip; i++) {

      // Determine Jacobian for this integration point
      this->shapeFunction(pts[i][0], pts[i][1]);

      // Interpolate strains
      //eps = B*u;
      VectorND<3> eps{};
      for (int beta = 0; beta < NEN; beta++) {
          eps[0] += shp[0][beta]*u[0][beta];
          eps[1] += shp[1][beta]*u[1][beta];
          eps[2] += shp[0][beta]*u[1][beta] + shp[1][beta]*u[0][beta];
      }

      // Set the material strain
      ret += theMaterial[i]->setTrialStrain(eps);
  }

  return ret;
}


const Matrix&
NineNodeQuad::getTangentStiff()
{

  K.Zero();

  double DB[3][2];

  // Loop over the integration points
  for (int i = 0; i < nip; i++) {

    // Determine Jacobian for this integration point
    double dvol = this->shapeFunction(pts[i][0], pts[i][1]);
    dvol *= (thickness*wts[i]);

    // Get the material tangent
    const Matrix &D = theMaterial[i]->getTangent();

    // Perform numerical integration
    //K = K + (B^ D * B) * intWt(i)*intWt(j) * detJ;
    //K.addMatrixTripleProduct(1.0, B, D, intWt(i)*intWt(j)*detJ);

    double D00 = D(0,0); double D01 = D(0,1); double D02 = D(0,2);
    double D10 = D(1,0); double D11 = D(1,1); double D12 = D(1,2);
    double D20 = D(2,0); double D21 = D(2,1); double D22 = D(2,2);

    //      for (int beta = 0, ib = 0, colIb =0, colIbP1 = 8;
    //   beta < 4;
    //   beta++, ib += 2, colIb += 16, colIbP1 += 16) {

    for (int alpha = 0, ia = 0; alpha < NEN; alpha++, ia += 2) {
      for (int beta = 0, ib = 0; beta < NEN; beta++, ib += 2) {

        DB[0][0] = dvol * (D00 * shp[0][beta] + D02 * shp[1][beta]);
        DB[1][0] = dvol * (D10 * shp[0][beta] + D12 * shp[1][beta]);
        DB[2][0] = dvol * (D20 * shp[0][beta] + D22 * shp[1][beta]);
        DB[0][1] = dvol * (D01 * shp[1][beta] + D02 * shp[0][beta]);
        DB[1][1] = dvol * (D11 * shp[1][beta] + D12 * shp[0][beta]);
        DB[2][1] = dvol * (D21 * shp[1][beta] + D22 * shp[0][beta]);

        K(ia,ib) += shp[0][alpha]*DB[0][0] + shp[1][alpha]*DB[2][0];
        K(ia,ib+1) += shp[0][alpha]*DB[0][1] + shp[1][alpha]*DB[2][1];
        K(ia+1,ib) += shp[1][alpha]*DB[1][0] + shp[0][alpha]*DB[2][0];
        K(ia+1,ib+1) += shp[1][alpha]*DB[1][1] + shp[0][alpha]*DB[2][1];
        //          matrixData[colIb   +   ia] += shp[0][alpha]*DB[0][0] + shp[1][alpha]*DB[2][0];
        //matrixData[colIbP1 +   ia] += shp[0][alpha]*DB[0][1] + shp[1][alpha]*DB[2][1];
        //matrixData[colIb   + ia+1] += shp[1][alpha]*DB[1][0] + shp[0][alpha]*DB[2][0];
        //matrixData[colIbP1 + ia+1] += shp[1][alpha]*DB[1][1] + shp[0][alpha]*DB[2][1];

      }
    }
  }

  return K;
}


const Matrix&
NineNodeQuad::getInitialStiff()
{
  if (Ki != 0)
    return *Ki;

  K.Zero();

  double dvol;
  double DB[3][2];

  // Loop over the integration points
  for (int i = 0; i < nip; i++) {

    // Determine Jacobian for this integration point
    dvol = this->shapeFunction(pts[i][0], pts[i][1]);
    dvol *= (thickness*wts[i]);

    // Get the material tangent
    const Matrix &D = theMaterial[i]->getInitialTangent();

    double D00 = D(0,0); double D01 = D(0,1); double D02 = D(0,2);
    double D10 = D(1,0); double D11 = D(1,1); double D12 = D(1,2);
    double D20 = D(2,0); double D21 = D(2,1); double D22 = D(2,2);

    // Perform numerical integration
    //K = K + (B^ D * B) * intWt(i)*intWt(j) * detJ;
    //K.addMatrixTripleProduct(1.0, B, D, intWt(i)*intWt(j)*detJ);
    for (int beta = 0, ib = 0, colIb =0, colIbP1 = 8;
     beta < NEN;
     beta++, ib += 2, colIb += 16, colIbP1 += 16) {

      for (int alpha = 0, ia = 0; alpha < NEN; alpha++, ia += 2) {

        DB[0][0] = dvol * (D00 * shp[0][beta] + D02 * shp[1][beta]);
        DB[1][0] = dvol * (D10 * shp[0][beta] + D12 * shp[1][beta]);
        DB[2][0] = dvol * (D20 * shp[0][beta] + D22 * shp[1][beta]);
        DB[0][1] = dvol * (D01 * shp[1][beta] + D02 * shp[0][beta]);
        DB[1][1] = dvol * (D11 * shp[1][beta] + D12 * shp[0][beta]);
        DB[2][1] = dvol * (D21 * shp[1][beta] + D22 * shp[0][beta]);

        matrixData[colIb   +   ia] += shp[0][alpha]*DB[0][0] + shp[1][alpha]*DB[2][0];
        matrixData[colIbP1 +   ia] += shp[0][alpha]*DB[0][1] + shp[1][alpha]*DB[2][1];
        matrixData[colIb   + ia+1] += shp[1][alpha]*DB[1][0] + shp[0][alpha]*DB[2][0];
        matrixData[colIbP1 + ia+1] += shp[1][alpha]*DB[1][1] + shp[0][alpha]*DB[2][1];
      }
    }
  }

  Ki = new Matrix(K);
  return K;
}

const Matrix&
NineNodeQuad::getMass()
{
    K.Zero();

    int i;
    static double rhoi[nip];
    double sum = 0.0;
    for (int i = 0; i < nip; i++) {
      if (rho == 0)
        rhoi[i] = theMaterial[i]->getRho();
      else
        rhoi[i] = rho;
      sum += rhoi[i];
    }

    if (sum == 0.0)
      return K;

    double rhodvol, Nrho;

    // Compute a lumped mass matrix
    for (int i = 0; i < nip; i++) {
        // Determine Jacobian for this integration point
        rhodvol = this->shapeFunction(pts[i][0], pts[i][1]);

        // Element plus material density ... MAY WANT TO REMOVE ELEMENT DENSITY
        rhodvol *= (rhoi[i]*thickness*wts[i]);

        for (int alpha = 0, ia = 0; alpha < NEN; alpha++, ia++) {
            Nrho = shp[2][alpha]*rhodvol;
            K(ia,ia) += Nrho;
            ia++;
            K(ia,ia) += Nrho;
        }
    }

    return K;
}

void
NineNodeQuad::zeroLoad()
{
  Q.Zero();

  applyLoad = 0;

  appliedB[0] = 0.0;
  appliedB[1] = 0.0;

  return;
}

int
NineNodeQuad::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    // Added option for applying body forces in load pattern: C.McGann, U.Washington
    int type;
    const Vector &data = theLoad->getData(type, loadFactor);

    if (type == LOAD_TAG_SelfWeight) {
        applyLoad = 1;
        appliedB[0] += loadFactor*data(0)*b[0];
        appliedB[1] += loadFactor*data(1)*b[1];
        return 0;
    } else {
        opserr << "NineNodeQuad::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
        return -1;
    }

    return -1;
}

int
NineNodeQuad::addInertiaLoadToUnbalance(const Vector &accel)
{
  int i;
  static double rhoi[nip];
  double sum = 0.0;
  for (int i = 0; i < nip; i++) {
    rhoi[i] = theMaterial[i]->getRho();
    sum += rhoi[i];
  }

  if (sum == 0.0)
    return 0;

  // Get R * accel from the nodes
  const Vector &Raccel1 = theNodes[0]->getRV(accel);
  const Vector &Raccel2 = theNodes[1]->getRV(accel);
  const Vector &Raccel3 = theNodes[2]->getRV(accel);
  const Vector &Raccel4 = theNodes[3]->getRV(accel);
  const Vector &Raccel5 = theNodes[4]->getRV(accel);
  const Vector &Raccel6 = theNodes[5]->getRV(accel);
  const Vector &Raccel7 = theNodes[6]->getRV(accel);
  const Vector &Raccel8 = theNodes[7]->getRV(accel);
  const Vector &Raccel9 = theNodes[8]->getRV(accel);

  if (2 != Raccel1.Size() || 2 != Raccel2.Size() || 2 != Raccel3.Size() ||
      2 != Raccel4.Size() || 2 != Raccel5.Size() || 2 != Raccel6.Size() ||
      2 != Raccel7.Size() || 2 != Raccel8.Size() || 2 != Raccel9.Size()) {
    opserr << "NineNodeQuad::addInertiaLoadToUnbalance matrix and vector sizes are incompatible\n";
    return -1;
  }

  static double ra[2*NEN];

  ra[0] = Raccel1(0);
  ra[1] = Raccel1(1);
  ra[2] = Raccel2(0);
  ra[3] = Raccel2(1);
  ra[4] = Raccel3(0);
  ra[5] = Raccel3(1);
  ra[6] = Raccel4(0);
  ra[7] = Raccel4(1);
  ra[8] = Raccel5(0);
  ra[9] = Raccel5(1);
  ra[10] = Raccel6(0);
  ra[11] = Raccel6(1);
  ra[12] = Raccel7(0);
  ra[13] = Raccel7(1);
  ra[14] = Raccel8(0);
  ra[15] = Raccel8(1);
  ra[16] = Raccel9(0);
  ra[17] = Raccel9(1);

  // Compute mass matrix
  const Matrix& M = this->getMass();

  // Want to add ( - fact * M R * accel ) to unbalance
  // Take advantage of lumped mass matrix
  for (int i = 0; i < NDF*NEN; i++)
    Q[i] += -M(i,i)*ra[i];

  return 0;
}

const Vector&
NineNodeQuad::getResistingForce()
{
    P.Zero();

    double dvol;

    // Loop over the integration points
    for (int i = 0; i < nip; i++) {

        // Determine Jacobian for this integration point
        dvol = this->shapeFunction(pts[i][0], pts[i][1]);
        dvol *= (thickness*wts[i]);

        // Get material stress response
        const Vector &sigma = theMaterial[i]->getStress();

        // Perform numerical integration on internal force
        //P = P + (B^ sigma) * intWt(i)*intWt(j) * detJ;
        //P.addMatrixTransposeVector(1.0, B, sigma, intWt(i)*intWt(j)*detJ);
        for (int alpha = 0, ia = 0; alpha < NEN; alpha++, ia += 2) {

            P(ia) += dvol*(shp[0][alpha]*sigma(0) + shp[1][alpha]*sigma(2));

            P(ia+1) += dvol*(shp[1][alpha]*sigma(1) + shp[0][alpha]*sigma(2));

            // Subtract equiv. body forces from the nodes
            //P = P - (N^ b) * intWt(i)*intWt(j) * detJ;
            //P.addMatrixTransposeVector(1.0, N, b, -intWt(i)*intWt(j)*detJ);
            if (applyLoad == 0) {
                P(ia)   -= dvol*(shp[2][alpha]*b[0]);
                P(ia+1) -= dvol*(shp[2][alpha]*b[1]);
            } else {
                P(ia)   -= dvol*(shp[2][alpha]*appliedB[0]);
                P(ia+1) -= dvol*(shp[2][alpha]*appliedB[1]);
            }
        }
    }

    // Subtract pressure loading from resisting force
    if (pressure != 0.0) {
        //P = P - pressureLoad;
        P.addVector(1.0, pressureLoad, -1.0);
    }

    // Subtract other external nodal loads ... P_res = P_int - P_ext
    //P = P - Q;
    P.addVector(1.0, Q, -1.0);

    return P;
}

const Vector&
NineNodeQuad::getResistingForceIncInertia()
{
    double rhoi[nip];
    double sum = 0.0;
    for (int i = 0; i < nip; i++) {
      rhoi[i] = theMaterial[i]->getRho();
      sum += rhoi[i];
    }

    // if no mass terms .. just add damping terms
    if (sum == 0.0) {
      this->getResistingForce();

      // add the damping forces if rayleigh damping
      if (betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
        P += this->getRayleighDampingForces();

      return P;
    }

    double a[NDF*NEN];
    for (int i=0; i<NEN; i++) {
        const Vector &accel = theNodes[i]->getTrialAccel();
        for (int j=0; j<NDF; j++)
          a[i*NDF+j] = accel[j];
    }

    // Compute the current resisting force
    this->getResistingForce();

    // Compute the mass matrix
    const Matrix& M = this->getMass();

    // Take advantage of lumped mass matrix
    for (int i = 0; i < NDF*NEN; i++)
        P(i) += M(i,i)*a[i];

    // add the damping forces if rayleigh damping
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
      P += this->getRayleighDampingForces();

    return P;
}

int
NineNodeQuad::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();

  // Quad packs its data into a Vector and sends this to theChannel
  // along with its dbTag and the commitTag passed in the arguments
  static Vector data(9);
  data(0) = this->getTag();
  data(1) = thickness;
  data(2) = b[0];
  data(3) = b[1];
  data(4) = pressure;

  data(5) = alphaM;
  data(6) = betaK;
  data(7) = betaK0;
  data(8) = betaKc;

  res += theChannel.sendVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING NineNodeQuad::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return res;
  }


  // Now quad sends the ids of its materials
  int matDbTag;

  static ID idData(2*nip+NEN);

  for (int i = 0; i < nip; i++) {
    idData(i) = theMaterial[i]->getClassTag();
    matDbTag = theMaterial[i]->getDbTag();
    // NOTE: we do have to ensure that the material has a database
    // tag if we are sending to a database channel.
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
            if (matDbTag != 0)
              theMaterial[i]->setDbTag(matDbTag);
    }
    idData(i+nip) = matDbTag;
  }

  for (int i = 0; i < NEN; i++)
    idData(2*nip+i) = connectedExternalNodes(i);

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING NineNodeQuad::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  // Finally, quad asks its material objects to send themselves
  for (int i = 0; i < nip; i++) {
    res += theMaterial[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING NineNodeQuad::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }

  return res;
}

int
NineNodeQuad::recvSelf(int commitTag, Channel &theChannel,
                       FEM_ObjectBroker &theBroker)
{
  int res = 0;

  int dataTag = this->getDbTag();

  // Quad creates a Vector, receives the Vector and then sets the
  // internal data with the data in the Vector
  static Vector data(9);
  res += theChannel.recvVector(dataTag, commitTag, data);
  if (res < 0) {
    opserr << "WARNING NineNodeQuad::recvSelf() - failed to receive Vector\n";
    return res;
  }

  this->setTag((int)data(0));
  thickness = data(1);
  b[0] = data(2);
  b[1] = data(3);
  pressure = data(4);

  alphaM = data(5);
  betaK = data(6);
  betaK0 = data(7);
  betaKc = data(8);

  static ID idData(2*nip+NEN);
  // Quad now receives the tags of its nine external nodes
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING NineNodeQuad::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  for( int i = 0; i < NEN; i++)
    connectedExternalNodes(i) = idData(2*nip+i);

  if (theMaterial == 0) {
    // Allocate new materials
    theMaterial = new NDMaterial *[nip];
    if (theMaterial == 0) {
      opserr << "NineNodeQuad::recvSelf() - Could not allocate NDMaterial* array\n";
      return -1;
    }
    for (int i = 0; i < nip; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+nip);
      // Allocate new material with the sent class tag
      theMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
      if (theMaterial[i] == 0) {
        opserr << "NineNodeQuad::recvSelf() - Broker could not create NDMaterial of class type " << matClassTag << endln;
        return -1;
      }
      // Now receive materials into the newly allocated space
      theMaterial[i]->setDbTag(matDbTag);
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "NineNodeQuad::recvSelf() - material " << i << "failed to recv itself\n";
        return res;
      }
    }
  }

  // materials exist , ensure materials of correct type and recvSelf on them
  else {
    for (int i = 0; i < nip; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+nip);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (theMaterial[i]->getClassTag() != matClassTag) {
        delete theMaterial[i];
        theMaterial[i] = theBroker.getNewNDMaterial(matClassTag);
        if (theMaterial[i] == 0) {
          opserr << "NineNodeQuad::recvSelf() - material " << i << "failed to create\n";
          return -1;
        }
      }
      // Receive the material
      theMaterial[i]->setDbTag(matDbTag);
      res += theMaterial[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "NineNodeQuad::recvSelf() - material " << i << "failed to recv itself\n";
        return res;
      }
    }
  }

  return res;
}

void
NineNodeQuad::Print(OPS_Stream &s, int flag)
{
  if (flag == 2) {

    s << "#NineNodeQuad\n";

    int i;
    const int numNodes = NEN;
    const int nstress = nip ;

    for (i=0; i<numNodes; i++) {
      const Vector &nodeCrd = theNodes[i]->getCrds();
      // const Vector &nodeDisp = theNodes[i]->getDisp();
      s << "#NODE " << nodeCrd(0) << " " << nodeCrd(1) << " " << endln;
     }

    // spit out the section location & invoke print on the scetion
    const int numMaterials = nip;

    static Vector avgStress(nstress);
    static Vector avgStrain(nstress);
    avgStress.Zero();
    avgStrain.Zero();
    for (i=0; i<numMaterials; i++) {
      avgStress += theMaterial[i]->getStress();
      avgStrain += theMaterial[i]->getStrain();
    }
    avgStress /= numMaterials;
    avgStrain /= numMaterials;

    s << "#AVERAGE_STRESS ";
    for (i=0; i<nstress; i++)
      s << avgStress(i) << " " ;
    s << endln;

    s << "#AVERAGE_STRAIN ";
    for (i=0; i<nstress; i++)
      s << avgStrain(i) << " " ;
    s << endln;
  }

  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "\nNineNodeQuad, element id:  " << this->getTag() << endln;
    s << "\tConnected external nodes:  " << connectedExternalNodes;
    s << "\tthickness:  " << thickness << endln;
    s << "\tsurface pressure:  " << pressure << endln;
    s << "\tmass density:  " << rho << endln;
    s << "\tbody forces:  " << b[0] << " " << b[1] << endln;
    theMaterial[0]->Print(s,flag);
    s << "\tStress (xx yy xy)" << endln;
    for (int i = 0; i < nip; i++)
        s << "\t\tGauss point " << i+1 << ": " << theMaterial[i]->getStress();
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    const ID& node_tags = this->getExternalNodes();
    s << OPS_PRINT_JSON_ELEM_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"" << this->getClassType() << "\", ";

    s << "\"nodes\": [";
    for (int i=0; i < NEN-1; i++)
        s << node_tags(i) << ", ";
    s << node_tags(NEN-1) << "]";
    s << ", ";

    s << "\"thickness\": " << thickness << ", ";
    s << "\"surfacePressure\": " << pressure << ", ";
    s << "\"masspervolume\": " << rho << ", ";
    s << "\"bodyForces\": [" << b[0] << ", " << b[1] << "], ";
    s << "\"material\": [";
    for (int i = 0; i < nip - 1; i++)
      s << theMaterial[i]->getTag() << ", ";
    s << theMaterial[nip - 1]->getTag() << "]";

    s << "}";
    return;
  }
}


Response*
NineNodeQuad::setResponse(const char **argv, int argc,
              OPS_Stream &output)
{
  Response *theResponse =0;

  output.tag("ElementOutput");
  output.attr("eleType","NineNodeQuad");
  output.attr("eleTag",this->getTag());
  output.attr("node1",connectedExternalNodes[0]);
  output.attr("node2",connectedExternalNodes[1]);
  output.attr("node3",connectedExternalNodes[2]);
  output.attr("node4",connectedExternalNodes[3]);
  output.attr("node5",connectedExternalNodes[4]);
  output.attr("node6",connectedExternalNodes[5]);
  output.attr("node7",connectedExternalNodes[6]);
  output.attr("node8",connectedExternalNodes[7]);
  output.attr("node9",connectedExternalNodes[8]);

  char dataOut[10];
  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0) {

    for (int i=1; i<=nip; i++) {
      sprintf(dataOut,"P1_%d",i);
      output.tag("ResponseType",dataOut);
      sprintf(dataOut,"P2_%d",i);
      output.tag("ResponseType",dataOut);
    }

    theResponse =  new ElementResponse(this, 1, P);
  }

  else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"integrPoint") == 0) {

    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= nip) {

      output.tag("GaussPoint");
      output.attr("number",pointNum);
      output.attr("eta",pts[pointNum-1][0]);
      output.attr("neta",pts[pointNum-1][1]);

      theResponse =  theMaterial[pointNum-1]->setResponse(&argv[2], argc-2, output);

      output.endTag();

    }
  }
  else if ((strcmp(argv[0],"stresses") ==0) || (strcmp(argv[0],"stress") ==0)) {
    for (int i=0; i<nip; i++) {
      output.tag("GaussPoint");
      output.attr("number",i+1);
      output.attr("eta",pts[i][0]);
      output.attr("neta",pts[i][1]);

      output.tag("NdMaterialOutput");
      output.attr("classType", theMaterial[i]->getClassTag());
      output.attr("tag", theMaterial[i]->getTag());

      output.tag("ResponseType","sigma11");
      output.tag("ResponseType","sigma22");
      output.tag("ResponseType","sigma12");

      output.endTag(); // GaussPoint
      output.endTag(); // NdMaterialOutput
      }
    theResponse =  new ElementResponse(this, 3, Vector(3*nip));
  }

  else if ((strcmp(argv[0],"stressesAtNodes") ==0) || (strcmp(argv[0],"stressAtNodes") ==0)) {
    for (int i=0; i<NEN; i++) {
      output.tag("NodalPoint");
      output.attr("number",i+1);
      // output.attr("eta",pts[i][0]);
      // output.attr("neta",pts[i][1]);

      // output.tag("NdMaterialOutput");
      // output.attr("classType", theMaterial[i]->getClassTag());
      // output.attr("tag", theMaterial[i]->getTag());

      output.tag("ResponseType","sigma11");
      output.tag("ResponseType","sigma22");
      output.tag("ResponseType","sigma12");

      output.endTag(); // GaussPoint
      // output.endTag(); // NdMaterialOutput
      }
    theResponse =  new ElementResponse(this, 11, Vector(3*NEN));
  }

  else if ((strcmp(argv[0],"strain") ==0) || 
           (strcmp(argv[0],"strains") ==0)) {
    for (int i=0; i<nip; i++) {
      output.tag("GaussPoint");
      output.attr("number",i+1);
      output.attr("eta", pts[i][0]);
      output.attr("neta",pts[i][1]);

      output.tag("NdMaterialOutput");
      output.attr("classType", theMaterial[i]->getClassTag());
      output.attr("tag", theMaterial[i]->getTag());

      output.tag("ResponseType","eta11");
      output.tag("ResponseType","eta22");
      output.tag("ResponseType","eta12");

      output.endTag(); // GaussPoint
      output.endTag(); // NdMaterialOutput
      }
    theResponse =  new ElementResponse(this, 4, Vector(3*nip));
  }

  output.endTag(); // ElementOutput

  return theResponse;
}

int
NineNodeQuad::getResponse(int responseID, Information &eleInfo)
{
  if (responseID == 1) {

    return eleInfo.setVector(this->getResistingForce());

  } else if (responseID == 3) {

    // Loop over the integration points
    static Vector stresses(3*nip);
    int cnt = 0;
    for (int i = 0; i < nip; i++) {

      // Get material stress response
      const Vector &sigma = theMaterial[i]->getStress();
      stresses(cnt) = sigma(0);
      stresses(cnt+1) = sigma(1);
      stresses(cnt+2) = sigma(2);
      cnt += 3;
    }

    return eleInfo.setVector(stresses);

  } else if (responseID == 11) {

    // extrapolate stress from Gauss points to element nodes
    static Vector stressGP(3*nip);
    static Vector stressAtNodes(3*NEN);
    stressAtNodes.Zero();
    int cnt = 0;
    // first get stress components (xx, yy, xy) at Gauss points
    for (int i = 0; i < nip; i++) {
      // Get material stress response
      const Vector &sigma = theMaterial[i]->getStress();
      stressGP(cnt) = sigma(0);
      stressGP(cnt+1) = sigma(1);
      stressGP(cnt+2) = sigma(2);
      cnt += 3;
    }

    const double We[NEN][nip] =  {{2.1869398183909485, 0.2777777777777778, 0.0352824038312731, 0.2777777777777778, -0.9858870384674904, -0.1252240726436203, -0.1252240726436203, -0.9858870384674904, 0.4444444444444444},
                                  {0.2777777777777778, 2.1869398183909485, 0.2777777777777778, 0.0352824038312731, -0.9858870384674904, -0.9858870384674904, -0.1252240726436203, -0.1252240726436203, 0.4444444444444444},
                                  {0.0352824038312731, 0.2777777777777778, 2.1869398183909485, 0.2777777777777778, -0.1252240726436203, -0.9858870384674904, -0.9858870384674904, -0.1252240726436203, 0.4444444444444444},
                                  {0.2777777777777778, 0.0352824038312731, 0.2777777777777778, 2.1869398183909485, -0.1252240726436203, -0.1252240726436203, -0.9858870384674904, -0.9858870384674904, 0.4444444444444444},
                                  {0.0, 0.0, 0.0, 0.0, 1.478830557701236,  0.0, 0.1878361089654305, 0.0, -0.6666666666666667},
                                  {0.0, 0.0, 0.0, 0.0, 0.0, 1.478830557701236,  0.0, 0.1878361089654305, -0.6666666666666667},
                                  {0.0, 0.0, 0.0, 0.0, 0.1878361089654305, 0.0, 1.478830557701236,  0.0, -0.6666666666666667},
                                  {0.0, 0.0, 0.0, 0.0, 0.0, 0.1878361089654305, 0.0, 1.478830557701236, -0.6666666666666667},
                                  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0}};

    for (int i = 0; i < NEN; i++) {
      for (int k = 0; k < 3; k++) {
        int p = 3*i + k;
        for (int j = 0; j < nip; j++) {
          int l = 3*j + k;
          stressAtNodes(p) += We[i][j] * stressGP(l);
          // opserr << "stressAtNodes(" << p << ") = We[" << i << "][" << j << "] * stressGP(" << l << ") = " << We[i][j] << " * " << stressGP(l) << " = " << stressAtNodes(p) <<  "\n";
        }
      }
    }

    return eleInfo.setVector(stressAtNodes);

  } else if (responseID == 4) {

    // Loop over the integration points
    static Vector stresses(3*nip);
    int cnt = 0;
    for (int i = 0; i < nip; i++) {

      // Get material stress response
      const Vector &sigma = theMaterial[i]->getStrain();
      stresses(cnt) = sigma(0);
      stresses(cnt+1) = sigma(1);
      stresses(cnt+2) = sigma(2);
      cnt += 3;
    }

    return eleInfo.setVector(stresses);

  } else

    return -1;
}

int
NineNodeQuad::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  int res = -1;

  // quad pressure loading
  if (strcmp(argv[0],"pressure") == 0) {
    return param.addObject(2, this);
  }
  // a material parameter
  else if ((strstr(argv[0],"material") != 0) && (strcmp(argv[0],"materialState") != 0)) {

    if (argc < 3)
      return -1;

    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= nip)
      return theMaterial[pointNum-1]->setParameter(&argv[2], argc-2, param);
    else
      return -1;
  }

  // otherwise it could be just a forall material parameter
  else {
    int matRes = res;
    for (int i=0; i<nip; i++) {

      matRes =  theMaterial[i]->setParameter(argv, argc, param);

      if (matRes != -1)
        res = matRes;
    }
  }

  return res;
}

int
NineNodeQuad::updateParameter(int parameterID, Information &info)
{
    int res = -1;
        int matRes = res;
  switch (parameterID) {
    case -1:
      return -1;

    case 1:

        for (int i = 0; i<nip; i++) {
        matRes = theMaterial[i]->updateParameter(parameterID, info);
        }
        if (matRes != -1) {
            res = matRes;
        }
        return res;

    case 2:
        pressure = info.theDouble;
        this->setPressureLoadAtNodes();    // update consistent nodal loads
        return 0;

    default:
      /*
      if (parameterID >= 100) { // material parameter
        int pointNum = parameterID/100;
        if (pointNum > 0 && pointNum <= 4)
          return theMaterial[pointNum-1]->updateParameter(parameterID-100*pointNum, info);
        else
          return -1;
      } else // unknown
      */
        return -1;
  }
}

double NineNodeQuad::shapeFunction(double s, double t)
{
    const Vector &nd1Crds = theNodes[0]->getCrds();
    const Vector &nd2Crds = theNodes[1]->getCrds();
    const Vector &nd3Crds = theNodes[2]->getCrds();
    const Vector &nd4Crds = theNodes[3]->getCrds();
    const Vector &nd5Crds = theNodes[4]->getCrds();
    const Vector &nd6Crds = theNodes[5]->getCrds();
    const Vector &nd7Crds = theNodes[6]->getCrds();
    const Vector &nd8Crds = theNodes[7]->getCrds();
    const Vector &nd9Crds = theNodes[8]->getCrds();

    double N1 =  (1-s)*(1-t)*s*t/4;
    double N2 = -(1+s)*(1-t)*s*t/4;
    double N3 =  (1+s)*(1+t)*s*t/4;
    double N4 = -(1-s)*(1+t)*s*t/4;

    double N5 = -(1-s*s)*(1-t)*t/2;
    double N6 =  (1+s)*(1-t*t)*s/2;
    double N7 =  (1-s*s)*(1+t)*t/2;
    double N8 = -(1-s)*(1-t*t)*s/2;

    double N9 = (1-s*s)*(1-t*t);

    // alternative
    // double N9 = (1-s*s)*(1-t*t);
    // double N1 = -(1-s)*(1-t)*(1+s+t)/4 + N9/4;
    // double N2 = -(1+s)*(1-t)*(1-s+t)/4 + N9/4;
    // double N3 = -(1+s)*(1+t)*(1-s-t)/4 + N9/4;
    // double N4 = -(1-s)*(1+t)*(1+s-t)/4 + N9/4;

    // double N5 = (1-s*s)*(1-t)/2 - N9/2;
    // double N6 = (1+s)*(1-t*t)/2 - N9/2;
    // double N7 = (1-s*s)*(1+t)/2 - N9/2;
    // double N8 = (1-s)*(1-t*t)/2 - N9/2;

    shp[2][0] = N1;
    shp[2][1] = N2;
    shp[2][2] = N3;
    shp[2][3] = N4;
    shp[2][4] = N5;
    shp[2][5] = N6;
    shp[2][6] = N7;
    shp[2][7] = N8;
    shp[2][8] = N9;

    // derivatives
    double N91 = -2*s*(1-t*t);
    double N92 = -2*t*(1-s*s);

    double N11 =  t*(1-t)*(1-2*s)/4;
    double N12 =  s*(1-s)*(1-2*t)/4;
    double N21 = -t*(1-t)*(1+2*s)/4;
    double N22 = -s*(1+s)*(1-2*t)/4;
    double N31 =  t*(1+t)*(1+2*s)/4;
    double N32 =  s*(1+s)*(1+2*t)/4;
    double N41 = -t*(1+t)*(1-2*s)/4;
    double N42 = -s*(1-s)*(1+2*t)/4;
    double N51 =  s*t*(1-t);
    double N52 = -(1-s*s)*(1-2*t)/2;
    double N61 =  (1-t*t)*(1+2*s)/2;
    double N62 = -s*t*(1+s);
    double N71 = -s*t*(1+t);
    double N72 =  (1-s*s)*(1+2*t)/2;
    double N81 =  -(1-t*t)*(1-2*s)/2;
    double N82 =  s*t*(1-s);

    // remains the same 2x2 matrix Cook Malkus Plesha 6.6, p. 180
    double J[2][2];

    J[0][0] = nd1Crds(0)*N11 + nd2Crds(0)*N21 + nd3Crds(0)*N31 + nd4Crds(0)*N41 +
              nd5Crds(0)*N51 + nd6Crds(0)*N61 + nd7Crds(0)*N71 + nd8Crds(0)*N81 + nd9Crds(0)*N91;

    J[0][1] = nd1Crds(0)*N12 + nd2Crds(0)*N22 + nd3Crds(0)*N32 + nd4Crds(0)*N42 +
              nd5Crds(0)*N52 + nd6Crds(0)*N62 + nd7Crds(0)*N72 + nd8Crds(0)*N82 + nd9Crds(0)*N92;

    J[1][0] = nd1Crds(1)*N11 + nd2Crds(1)*N21 + nd3Crds(1)*N31 + nd4Crds(1)*N41 +
              nd5Crds(1)*N51 + nd6Crds(1)*N61 + nd7Crds(1)*N71 + nd8Crds(1)*N81 + nd9Crds(1)*N91;

    J[1][1] = nd1Crds(1)*N12 + nd2Crds(1)*N22 + nd3Crds(1)*N32 + nd4Crds(1)*N42 +
              nd5Crds(1)*N52 + nd6Crds(1)*N62 + nd7Crds(1)*N72 + nd8Crds(1)*N82 + nd9Crds(1)*N92;

    double detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];

    double oneOverdetJ = 1/detJ;
    double L[2][2];

    // L = inv(J)
    L[0][0] =  J[1][1]*oneOverdetJ;
    L[1][0] = -J[0][1]*oneOverdetJ;
    L[0][1] = -J[1][0]*oneOverdetJ;
    L[1][1] =  J[0][0]*oneOverdetJ;

    double L00 = L[0][0];
    double L10 = L[1][0];
    double L01 = L[0][1];
    double L11 = L[1][1];

    shp[0][0] = L00*N11 + L01*N12;
    shp[0][1] = L00*N21 + L01*N22;
    shp[0][2] = L00*N31 + L01*N32;
    shp[0][3] = L00*N41 + L01*N42;
    shp[0][4] = L00*N51 + L01*N52;
    shp[0][5] = L00*N61 + L01*N62;
    shp[0][6] = L00*N71 + L01*N72;
    shp[0][7] = L00*N81 + L01*N82;
    shp[0][8] = L00*N91 + L01*N92;

    shp[1][0] = L10*N11 + L11*N12;
    shp[1][1] = L10*N21 + L11*N22;
    shp[1][2] = L10*N31 + L11*N32;
    shp[1][3] = L10*N41 + L11*N42;
    shp[1][4] = L10*N51 + L11*N52;
    shp[1][5] = L10*N61 + L11*N62;
    shp[1][6] = L10*N71 + L11*N72;
    shp[1][7] = L10*N81 + L11*N82;
    shp[1][8] = L10*N91 + L11*N92;

    return detJ;
}

void
NineNodeQuad::setPressureLoadAtNodes()
{
    pressureLoad.Zero();

    if (pressure == 0.0)
        return;

    const Vector &node1 = theNodes[0]->getCrds();
    const Vector &node2 = theNodes[1]->getCrds();
    const Vector &node3 = theNodes[2]->getCrds();
    const Vector &node4 = theNodes[3]->getCrds();
    const Vector &node5 = theNodes[4]->getCrds();
    const Vector &node6 = theNodes[5]->getCrds();
    const Vector &node7 = theNodes[6]->getCrds();
    const Vector &node8 = theNodes[7]->getCrds();
    // center node has no pressure commponents
    // const Vector &node9 = theNodes[8]->getCrds();

    double x1 = node1(0);
    double y1 = node1(1);
    double x2 = node2(0);
    double y2 = node2(1);
    double x3 = node3(0);
    double y3 = node3(1);
    double x4 = node4(0);
    double y4 = node4(1);
    double x5 = node5(0);
    double y5 = node5(1);
    double x6 = node6(0);
    double y6 = node6(1);
    double x7 = node7(0);
    double y7 = node7(1);
    double x8 = node8(0);
    double y8 = node8(1);
    // double x9 = node9(0);
    // double y9 = node9(1);

    double dx15 = x5-x1;
    double dy15 = y5-y1;
    double dx52 = x2-x5;
    double dy52 = y2-y5;
    double dx26 = x6-x2;
    double dy26 = y6-y2;
    double dx63 = x3-x6;
    double dy63 = y3-y6;
    double dx37 = x7-x3;
    double dy37 = y7-y3;
    double dx74 = x4-x7;
    double dy74 = y4-y7;
    double dx48 = x8-x4;
    double dy48 = y8-y4;
    double dx81 = x1-x8;
    double dy81 = y1-y8;

    double fac1 = 0.3333333333333333;
    double fac2 = 0.6666666666666667;

    // Contribution from side 15
    pressureLoad(0) += pressure*fac1*dy15;
    pressureLoad(8) += pressure*fac2*dy15;
    pressureLoad(1) += pressure*fac1*-dx15;
    pressureLoad(9) += pressure*fac2*-dx15;

    // Contribution from side 52
    pressureLoad(8) += pressure*fac2*dy52;
    pressureLoad(2) += pressure*fac1*dy52;
    pressureLoad(9) += pressure*fac2*-dx52;
    pressureLoad(3) += pressure*fac1*-dx52;

    // Contribution from side 26
    pressureLoad( 2) += pressure*fac1*dy26;
    pressureLoad(10) += pressure*fac2*dy26;
    pressureLoad( 3) += pressure*fac1*-dx26;
    pressureLoad(11) += pressure*fac2*-dx26;

    // Contribution from side 63
    pressureLoad(10) += pressure*fac2*dy63;
    pressureLoad(4) += pressure*fac1*dy63;
    pressureLoad(11) += pressure*fac2*-dx63;
    pressureLoad(5) += pressure*fac1*-dx63;

    // Contribution from side 37
    pressureLoad( 4) += pressure*fac1*dy37;
    pressureLoad(12) += pressure*fac2*dy37;
    pressureLoad( 5) += pressure*fac1*-dx37;
    pressureLoad(13) += pressure*fac2*-dx37;

    // Contribution from side 74
    pressureLoad(12) += pressure*fac2*dy74;
    pressureLoad( 6) += pressure*fac1*dy74;
    pressureLoad(13) += pressure*fac2*-dx74;
    pressureLoad( 7) += pressure*fac1*-dx74;

    // Contribution from side 48
    pressureLoad( 6) += pressure*fac1*dy48;
    pressureLoad(14) += pressure*fac2*dy48;
    pressureLoad( 7) += pressure*fac1*-dx48;
    pressureLoad(15) += pressure*fac2*-dx48;

    // Contribution from side 81
    pressureLoad(14) += pressure*fac2*dy81;
    pressureLoad( 0) += pressure*fac1*dy81;
    pressureLoad(15) += pressure*fac2*-dx81;
    pressureLoad( 1) += pressure*fac1*-dx81;

    //pressureLoad = pressureLoad*thickness;
}
