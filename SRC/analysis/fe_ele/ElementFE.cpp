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
// Purpose: This file contains the code for implementing the methods
// of the ElementFE class interface.
//
// Written: cmp
// Created: Jan 2026
// Revision: A
//
#include <ElementFE.h>
#include <stdlib.h>
#include <assert.h>

#include <Element.h>
#include <Domain.h>
#include <Node.h>
#include <DOF_Group.h>
#include <Integrator.h>
#include <Subdomain.h>
#include <AnalysisModel.h>
#include <Matrix.h>
#include <Vector.h>

#define MAX_NUM_DOF 64

// static variables initialisation
Matrix **ElementFE::theMatrices; // pointers to class wide matrices
Vector **ElementFE::theVectors;  // pointers to class widde vectors
int ElementFE::numFEs(0);        // number of objects


// constructor that take the corresponding model element.
ElementFE::ElementFE(int tag, Element *ele)
  : FE_Element(tag, (ele->getExternalNodes()).Size(), ele->getNumDOF())
  , myID(ele->getNumDOF()),
   numDOF(ele->getNumDOF()), myEle(ele),
   theResidual(nullptr), theTangent(nullptr), theIntegrator(nullptr)
{
  assert(numDOF > 0);

  // get elements domain & check it is valid
  Domain *theDomain = ele->getDomain();
  assert(theDomain != nullptr);

  // keep a pointer to all DOF_Groups
  int numGroups = ele->getNumExternalNodes();
  const ID &nodes = ele->getExternalNodes();

  for (int i=0; i<numGroups; i++) {
    Node *nodePtr =theDomain->getNode(nodes(i));
    assert(nodePtr != nullptr);

    DOF_Group *dofGrpPtr = nodePtr->getDOF_GroupPtr();
    assert(dofGrpPtr != nullptr);
    myDOF_Groups(i) = dofGrpPtr->getTag();
  }

  // if this is the first ElementFE we now
  // create the arrays used to store pointers to class wide
  // matrix and vector objects used to return tangent and residual
  if (numFEs == 0) {
    theMatrices = new Matrix *[MAX_NUM_DOF+1];
    theVectors  = new Vector *[MAX_NUM_DOF+1];

    for (int i=0; i<MAX_NUM_DOF; i++) {
      theMatrices[i] = nullptr;
      theVectors[i]  = nullptr;
    }
  }

  // if Elements are not subdomains, set up pointers to
  // objects to return tangent Matrix and residual Vector.
  if (numDOF <= MAX_NUM_DOF) {
    // use class wide objects
    if (theVectors[numDOF] == nullptr) {
        theVectors[numDOF] = new Vector(numDOF);
        theMatrices[numDOF] = new Matrix(numDOF,numDOF);
        theResidual = theVectors[numDOF];
        theTangent  = theMatrices[numDOF];

    } else {
        theResidual = theVectors[numDOF];
        theTangent  = theMatrices[numDOF];
    }

  }

  numFEs++;
}


ElementFE::~ElementFE()
{
  // decrement number of ElementFEs
  numFEs--;

  // delete tangent and residual if created specially
  if (numDOF > MAX_NUM_DOF) {
      if (theTangent != nullptr)
        delete theTangent;
      if (theResidual != nullptr) 
        delete theResidual;
  }

  // if this is the last ElementFE, clean up the
  // storage for the matrix and vector objects
  if (numFEs == 0) {
      for (int i=0; i<MAX_NUM_DOF; i++) {
          if (theVectors[i] != nullptr)
              delete theVectors[i];
          if (theMatrices[i] != nullptr)
              delete theMatrices[i];
      }
      delete [] theMatrices;
      delete [] theVectors;
  }
}


const ID &
ElementFE::getID() const
{
  return myID;
}


// void setID(int index, int value);
//        Method to set the corresponding index of the ID to value.

int
ElementFE::setID(AnalysisModel &theModel)
{
  int current = 0;


  int numGrps = myDOF_Groups.Size();
  for (int i=0; i<numGrps; i++) {
    int tag = myDOF_Groups(i);

    DOF_Group *dofPtr = theModel.getDOF_GroupPtr(tag);
    assert(dofPtr != nullptr);

    const ID &theDOFid = dofPtr->getID();

    for (int j=0; j<theDOFid.Size(); j++) {
      assert(current < numDOF);
      myID(current++) = theDOFid(j);
    }

  }
  return 0;
}


const Matrix &
ElementFE::getTangent(Integrator *theNewIntegrator)
{
  theIntegrator = theNewIntegrator;

  assert(myEle != nullptr);

  if (theNewIntegrator != nullptr)
    theNewIntegrator->formEleTangent(this);

  return *theTangent;

}


void
ElementFE::zeroTangent()
{
  assert(myEle != nullptr);
  theTangent->Zero();
}


void
ElementFE::addKtToTang(double fact)
{
  assert (myEle != nullptr);

  // check for a quick return
  if (fact == 0.0)
      return;
  else
      theTangent->addMatrix(myEle->getTangentStiff(),fact);
}


void
ElementFE::addCtoTang(double fact)
{
  assert (myEle != nullptr);

  // check for a quick return
  if (fact == 0.0)
    return;
  else
    theTangent->addMatrix(myEle->getDamp(),fact);
}


void
ElementFE::addMtoTang(double fact)
{
  if (myEle != nullptr) {

    // check for a quick return
    if (fact == 0.0)
      return;
    else
      theTangent->addMatrix(myEle->getMass(),fact);
  }
}


void
ElementFE::addKiToTang(double fact)
{
  // check for a quick return
  if (fact == 0.0)
    return;

  else // if (myEle->isSubdomain() == false)
    theTangent->addMatrix(myEle->getInitialStiff(), fact);
}


void
ElementFE::addKpToTang(double fact, int numP)
{
  if (myEle != nullptr) {
    // check for a quick return
    if (fact == 0.0)
      return;

    else if (myEle->isSubdomain() == false) {
      const Matrix *thePrevMat = myEle->getPreviousK(numP);
      if (thePrevMat != nullptr)
        theTangent->addMatrix(*thePrevMat, fact);

    }
  }
}

int
ElementFE::storePreviousK(int numP)
{
  int res = 0;
  if (myEle != nullptr)
    res = myEle->storePreviousK(numP);

  return res;
}

//
// RESIDUAL
//
const Vector &
ElementFE::getResidual(Integrator *theNewIntegrator)
{
  theIntegrator = theNewIntegrator;

  if (theIntegrator == nullptr)
    return *theResidual;

  assert(myEle != nullptr);

  theNewIntegrator->formEleResidual(this);
  return *theResidual;
}

void
ElementFE::zeroResidual()
{
  assert(myEle != nullptr);
  assert(myEle->isSubdomain() == false);

  theResidual->Zero();
}


void
ElementFE::addRtoResidual(double fact)
{
  assert(myEle != nullptr);
  assert(myEle->isSubdomain() == false);

  // check for a quick return
  if (fact == 0.0)
    return;

  else {
    const Vector &eleResisting = myEle->getResistingForce();
    theResidual->addVector(1.0, eleResisting, -fact);
  }
}


void
ElementFE::addRIncInertiaToResidual(double fact)
{
  assert(myEle != nullptr);

  // check for a quick return
  if (fact == 0.0)
    return;

  else {
    const Vector &eleResisting = myEle->getResistingForceIncInertia();
    theResidual->addVector(1.0, eleResisting, -fact);
  }
}


const Vector &
ElementFE::getTangForce(const Vector &disp, double fact)
{
  assert(myEle != nullptr);

  // zero out the force vector
  theResidual->Zero();

  // check for a quick return
  if (fact == 0.0)
    return *theResidual;

  // get the components we need out of the vector
  // and place in a temporary vector
  Vector tmp(numDOF);
  for (int i=0; i<numDOF; i++) {
    int dof = myID(i);
    if (dof >= 0)
      tmp(i) = disp(myID(i));
    else
      tmp(i) = 0.0;
  }

  // form the tangent again and then add the force
  theIntegrator->formEleTangent(this);
  theResidual->addMatrixVector(*theTangent,tmp,fact);

  return *theResidual;
}



const Vector &
ElementFE::getK_Force(const Vector &disp, double fact)
{
  assert(myEle != nullptr);

  // zero out the force vector
  theResidual->Zero();

  // check for a quick return
  if (fact == 0.0)
      return *theResidual;

  // get the components we need out of the vector
  // and place in a temporary vector
  Vector tmp(numDOF);
  for (int i=0; i<numDOF; i++) {
    int dof = myID(i);
    if (dof >= 0)
      tmp(i) = disp(myID(i));
    else
      tmp(i) = 0.0;
  }

  theResidual->addMatrixVector(1.0, myEle->getTangentStiff(), tmp, fact);

  return *theResidual;
}


const Vector &
ElementFE::getKi_Force(const Vector &disp, double fact)
{
  assert(myEle != nullptr);

  // zero out the force vector
  theResidual->Zero();

  // check for a quick return
  if (fact == 0.0)
    return *theResidual;

  // get the components we need out of the vector
  // and place in a temporary vector
  Vector tmp(numDOF);
  for (int i=0; i<numDOF; i++) {
    int dof = myID(i);
    if (dof >= 0)
      tmp(i) = disp(myID(i));
    else
      tmp(i) = 0.0;
  }

  theResidual->addMatrixVector(myEle->getInitialStiff(), tmp, fact);

  return *theResidual;
}

const Vector &
ElementFE::getM_Force(const Vector &disp, double fact)
{
  assert(myEle != nullptr);

  // zero out the force vector
  theResidual->Zero();

  // check for a quick return
  if (fact == 0.0)
      return *theResidual;

  // get the components we need out of the vector
  // and place in a temporary vector
  Vector tmp(numDOF);
  for (int i=0; i<numDOF; i++) {
    int dof = myID(i);
    if (dof >= 0)
      tmp(i) = disp(myID(i));
    else
      tmp(i) = 0.0;
  }

  theResidual->addMatrixVector(myEle->getMass(), tmp, fact);

  return *theResidual;
}

const Vector &
ElementFE::getC_Force(const Vector &disp, double fact)
{
  assert(myEle != nullptr);

  // zero out the force vector
  theResidual->Zero();

  // check for a quick return
  if (fact == 0.0)
      return *theResidual;

  // get the components we need out of the vector
  // and place in a temporary vector
  Vector tmp(numDOF);
  for (int i=0; i<numDOF; i++) {
    int dof = myID(i);
    if (dof >= 0)
      tmp(i) = disp(myID(i));
    else
      tmp(i) = 0.0;
  }

  theResidual->addMatrixVector(myEle->getDamp(), tmp, fact);

  return *theResidual;
}


Integrator *
ElementFE::getLastIntegrator()
{
  return theIntegrator;
}


const Vector &
ElementFE::getLastResponse()
{
  assert(myEle != nullptr);

  if (theIntegrator != nullptr) {
    if (theIntegrator->getLastResponse(*theResidual,myID) < 0) {
      opserr << "WARNING ElementFE::getLastResponse()";
      opserr << " - the Integrator had problems with getLastResponse()\n";
    }
  }
  else {
    theResidual->Zero();
    opserr << "WARNING  ElementFE::getLastResponse()";
    opserr << " No Integrator yet passed\n";
  }

  Vector &result = *theResidual;
  return result;
}

void
ElementFE::addM_Force(const Vector &accel, double fact)
{
  assert(myEle != nullptr);
  assert(myEle->isSubdomain() == false);

  // check for a quick return
  if (fact == 0.0)
    return;

  // get the components we need out of the vector
  // and place in a temporary vector
  Vector tmp(numDOF);
  for (int i=0; i<numDOF; i++) {
    int loc = myID(i);
    if (loc >= 0)
      tmp(i) = accel(loc);
    else
      tmp(i) = 0.0;
  }

  theResidual->addMatrixVector(1.0, myEle->getMass(), tmp, fact);
}

void
ElementFE::addD_Force(const Vector &accel, double fact)
{
  assert(myEle != nullptr);
  assert(myEle->isSubdomain() == false);

  // check for a quick return
  if (fact == 0.0)
    return;

  // get the components we need out of the vector
  // and place in a temporary vector
  Vector tmp(numDOF);
  for (int i=0; i<numDOF; i++) {
    int loc = myID(i);
    if (loc >= 0)
        tmp(i) = accel(loc);
    else
        tmp(i) = 0.0;
  }

  theResidual->addMatrixVector(myEle->getDamp(), tmp, fact);
}

void
ElementFE::addK_Force(const Vector &disp, double fact)
{
  assert(myEle != nullptr);
  assert(myEle->isSubdomain() == false);

  // check for a quick return
  if (fact == 0.0)
    return;

  // get the components we need out of the vector
  // and place in a temporary vector
  Vector tmp(numDOF);
  for (int i=0; i<numDOF; i++) {
    int loc = myID(i);
    if (loc >= 0)
        tmp(i) = disp(loc);
    else
        tmp(i) = 0.0;
  }

  theResidual->addMatrixVector(1.0, myEle->getTangentStiff(), tmp, fact);
}



void
ElementFE::addLocalM_Force(const Vector &accel, double fact)
{
  assert(myEle != nullptr);
  assert(myEle->isSubdomain() == false);

  // check for a quick return
  if (fact == 0.0)
    return;

  theResidual->addMatrixVector(myEle->getMass(), accel, fact);
}

void
ElementFE::addLocalD_Force(const Vector &accel, double fact)
{
  assert(myEle != nullptr);
  assert(myEle->isSubdomain() == false);

  // check for a quick return
  if (fact == 0.0)
      return;

  if (theResidual->addMatrixVector(myEle->getDamp(), accel, fact) < 0){
    opserr << "WARNING ElementFE::addLocalD_Force() - ";
    opserr << "- addMatrixVector returned error\n";
  }
}


Element *
ElementFE::getElement()
{
  return myEle;
}


// AddingSensitivity:BEGIN /////////////////////////////////
void
ElementFE::addResistingForceSensitivity(int gradNumber, double fact)
{
  theResidual->addVector(1.0, myEle->getResistingForceSensitivity(gradNumber), -fact);
}

void
ElementFE::addM_ForceSensitivity(int gradNumber, const Vector &vect, double fact)
{
  // Get the components we need out of the vector
  // and place in a temporary vector
  Vector tmp(numDOF);
  for (int i=0; i<numDOF; i++) {
    int loc = myID(i);
    if (loc >= 0) {
      tmp(i) = vect(loc);
    } else {
      tmp(i) = 0.0;
    }
  }
  theResidual->addMatrixVector(myEle->getMassSensitivity(gradNumber),tmp,fact);
}

void
ElementFE::addD_ForceSensitivity(int gradNumber, const Vector &vect, double fact)
{
  assert(myEle != nullptr);

  // check for a quick return
  if (fact == 0.0)
    return;

  if (myEle->isSubdomain() == false) {
    // get the components we need out of the vector
    // and place in a temporary vector
    Vector tmp(numDOF);
    for (int i=0; i<numDOF; i++) {
      int loc = myID(i);
      if (loc >= 0)
        tmp(i) = vect(loc);
      else
        tmp(i) = 0.0;
    }
    theResidual->addMatrixVector(myEle->getDampSensitivity(gradNumber), tmp, fact);
  }
  else {
    opserr << "WARNING ElementFE::addD_ForceSensitivity() - ";
    opserr << "- this should not be called on a Subdomain!\n";
  }
}

void
ElementFE::addLocalD_ForceSensitivity(int gradNumber, const Vector &accel, double fact)
{
  // check for a quick return
  if (fact == 0.0)
      return;
  theResidual->addMatrixVector(myEle->getDampSensitivity(gradNumber),  accel, fact);

}


void
ElementFE::addLocalM_ForceSensitivity(int gradNumber, const Vector &accel, double fact)
{
  assert(myEle != nullptr);
  assert(myEle->isSubdomain() == false);

  // check for a quick return
  if (fact == 0.0)
    return;

  if (theResidual->addMatrixVector(myEle->getMassSensitivity(gradNumber), accel, fact) < 0) {
    opserr << "WARNING ElementFE::addLocalD_ForceSensitivity() - ";
    opserr << "- addMatrixVector returned error\n";
  }
}





int
ElementFE::commitSensitivity(int gradNum, int numGrads)
{
  myEle->commitSensitivity(gradNum, numGrads);
  return 0;
}

// AddingSensitivity:END ////////////////////////////////////


int
ElementFE::updateElement()
{
  if (myEle != nullptr)
    return myEle->update();
  return 0;
}
