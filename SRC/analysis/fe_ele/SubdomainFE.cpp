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
// of the SubdomainFE class interface.
//
// Written: cmp
//
#include <SubdomainFE.h>
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
Matrix **SubdomainFE::theMatrices; // pointers to class wide matrices
Vector **SubdomainFE::theVectors;  // pointers to class widde vectors
int SubdomainFE::numFEs(0);           // number of objects


// constructor that take the corresponding model element.
SubdomainFE::SubdomainFE(int tag, Element *ele)
  : FE_Element(tag, (ele->getExternalNodes()).Size(), ele->getNumDOF())
  , myID(ele->getNumDOF()),
   numDOF(ele->getNumDOF()), myEle(ele),
   theResidual(nullptr), theTangent(nullptr), theIntegrator(nullptr)
{
  assert(numDOF > 0);
  opserr << "SubdomainFE::SubdomainFE - called\n";

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

  // if this is the first SubdomainFE we now
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

  {
    // as subdomains have own matrix for tangent and residual don't need
    // to set matrix and vector pointers to these objects
    theResidual = new Vector(numDOF);
      // invoke setFE_ElementPtr() method on Subdomain
    Subdomain *theSub = (Subdomain *)ele;
    theSub->setFE_ElementPtr(this);
  }

  // increment number of SubdomainFEs by 1
  numFEs++;
}





SubdomainFE::~SubdomainFE()
{
  // decrement number of SubdomainFEs
  numFEs--;

  // delete tangent and residual if created specially
  if (numDOF > MAX_NUM_DOF) {
    if (theTangent != nullptr)
      delete theTangent;
    if (theResidual != nullptr) 
      delete theResidual;
  }

  // if this is the last SubdomainFE, clean up the
  // storage for the matrix and vector objects
  if (numFEs == 0) {
    for (int i=0; i<MAX_NUM_DOF; i++) {
      if (theVectors && theVectors[i] != nullptr)
        delete theVectors[i];
      if (theMatrices && theMatrices[i] != nullptr)
        delete theMatrices[i];
    }
    delete [] theMatrices;
    delete [] theVectors;
  }
}


const ID &
SubdomainFE::getID() const
{
  return myID;
}


// void setID(int index, int value);
//        Method to set the corresponding index of the ID to value.

int
SubdomainFE::setID(AnalysisModel &theModel)
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
SubdomainFE::getTangent(Integrator *theNewIntegrator)
{
  theIntegrator = theNewIntegrator;

  {
    Subdomain *theSub = (Subdomain *)myEle;
    theSub->computeTang();
    return theSub->getTang();
  }
}

void
SubdomainFE::zeroTangent()
{
  opserr << "This should not be called on a Subdomain!\n";
}


void
SubdomainFE::addKtToTang(double fact)
{
  opserr << "WARNING FE_Element::addKToTang() - ";
  opserr << "- this should not be called on a Subdomain!\n";
}


void
SubdomainFE::addCtoTang(double fact)
{
#ifdef OLD_FE
  assert (myEle != nullptr);
  assert(myEle->isSubdomain() == false);

  // check for a quick return
  if (fact == 0.0)
    return;
  else
    theTangent->addMatrix(myEle->getDamp(),fact);
#endif
}


void
SubdomainFE::addMtoTang(double fact)
{
#ifdef OLD_FE
  if (myEle != nullptr) {
    assert(myEle->isSubdomain() == false);

    // check for a quick return
    if (fact == 0.0)
      return;
    else
      theTangent->addMatrix(myEle->getMass(),fact);
  }
#endif
}


void
SubdomainFE::addKiToTang(double fact)
{
#ifdef OLD_FE
  if (myEle != nullptr) {
    assert(myEle->isSubdomain() == false);

    // check for a quick return
    if (fact == 0.0)
      return;

    else // if (myEle->isSubdomain() == false)
      theTangent->addMatrix(myEle->getInitialStiff(), fact);
  }
#endif
}


void
SubdomainFE::addKpToTang(double fact, int numP)
{
#ifdef OLD_FE
  if (myEle != nullptr) {
    // check for a quick return
    if (fact == 0.0)
      return;

    else if (myEle->isSubdomain() == false) {
      const Matrix *thePrevMat = myEle->getPreviousK(numP);
      if (thePrevMat != nullptr)
        theTangent->addMatrix(*thePrevMat, fact);

    } else {
      opserr << "WARNING SubdomainFE::addKpToTang() - ";
      opserr << "- this should not be called on a Subdomain!\n";
    }
  }
#endif
}

int
SubdomainFE::storePreviousK(int numP)
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
SubdomainFE::getResidual(Integrator *theNewIntegrator)
{
  theIntegrator = theNewIntegrator;

#ifdef OLD_FE
  if (theIntegrator == nullptr)
    return *theResidual;

  assert(myEle != nullptr);

  {
    Subdomain *theSub = (Subdomain *)myEle;
    theSub->computeResidual();
    return theSub->getResistingForce();
  }
#else 
  return *theResidual;
#endif
}

void
SubdomainFE::zeroResidual()
{
#ifdef OLD_FE
  assert(myEle != nullptr);
  assert(myEle->isSubdomain() == false);

  theResidual->Zero();
#endif
}


void
SubdomainFE::addRtoResidual(double fact)
{
#ifdef OLD_FE
  assert(myEle != nullptr);
  assert(myEle->isSubdomain() == false);

  // check for a quick return
  if (fact == 0.0)
    return;

  else {
    const Vector &eleResisting = myEle->getResistingForce();
    theResidual->addVector(1.0, eleResisting, -fact);
  }
#endif
}


void
SubdomainFE::addRIncInertiaToResidual(double fact)
{
#ifdef OLD_FE
  assert(myEle != nullptr);
  assert(myEle->isSubdomain() == false);

  // check for a quick return
  if (fact == 0.0)
      return;

  else {
    const Vector &eleResisting = myEle->getResistingForceIncInertia();
    theResidual->addVector(1.0, eleResisting, -fact);
  }
#endif
}


const Vector &
SubdomainFE::getTangForce(const Vector &disp, double fact)
{
#ifdef OLD_FE
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

  if (myEle->isSubdomain() == false) {
    // form the tangent again and then add the force
    theIntegrator->formEleTangent(this);
    theResidual->addMatrixVector(*theTangent,tmp,fact);

  } else {
    theResidual->addMatrixVector(((Subdomain *)myEle)->getTang(),tmp,fact);
  }
#endif
  return *theResidual;
}



const Vector &
SubdomainFE::getK_Force(const Vector &disp, double fact)
{
#ifdef OLD_FE
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
#endif
  return *theResidual;
}


const Vector &
SubdomainFE::getKi_Force(const Vector &disp, double fact)
{
#ifdef OLD_FE
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
#endif
  return *theResidual;
}

const Vector &
SubdomainFE::getM_Force(const Vector &disp, double fact)
{
#ifdef OLD_FE
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
#endif
  return *theResidual;
}

const Vector &
SubdomainFE::getC_Force(const Vector &disp, double fact)
{
#ifdef OLD_FE
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
#endif
  return *theResidual;
}


Integrator *
SubdomainFE::getLastIntegrator()
{
  return theIntegrator;
}


const Vector &
SubdomainFE::getLastResponse()
{
#ifndef OLD_FE
  assert(myEle != nullptr);

  if (theIntegrator != nullptr) {
    if (theIntegrator->getLastResponse(*theResidual,myID) < 0) {
      opserr << "WARNING SubdomainFE::getLastResponse()";
      opserr << " - the Integrator had problems with getLastResponse()\n";
    }
  }
  else {
    theResidual->Zero();
    opserr << "WARNING  SubdomainFE::getLastResponse()";
    opserr << " No Integrator yet passed\n";
  }
#endif
  Vector &result = *theResidual;
  return result;
}

void
SubdomainFE::addM_Force(const Vector &accel, double fact)
{
#ifdef OLD_FE
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
#endif
}

void
SubdomainFE::addD_Force(const Vector &accel, double fact)
{
#ifdef OLD_FE
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
#endif
}

void
SubdomainFE::addK_Force(const Vector &disp, double fact)
{
#ifdef OLD_FE
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
#endif
}





Element *
SubdomainFE::getElement()
{
  return myEle;
}


// AddingSensitivity:BEGIN /////////////////////////////////
void
SubdomainFE::addResistingForceSensitivity(int gradNumber, double fact)
{
  theResidual->addVector(1.0, myEle->getResistingForceSensitivity(gradNumber), -fact);
}

void
SubdomainFE::addM_ForceSensitivity(int gradNumber, const Vector &vect, double fact)
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
  if (theResidual->addMatrixVector(myEle->getMassSensitivity(gradNumber),tmp,fact) < 0) {
    opserr << "WARNING SubdomainFE::addM_ForceSensitivity() - ";
    opserr << "- addMatrixVector returned error\n";
  }
}

void
SubdomainFE::addD_ForceSensitivity(int gradNumber, const Vector &vect, double fact)
{
#ifdef OLD_FE
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
    if (theResidual->addMatrixVector(myEle->getDampSensitivity(gradNumber), tmp, fact) < 0){
      opserr << "WARNING SubdomainFE::addD_ForceSensitivity() - ";
      opserr << "- addMatrixVector returned error\n";
    }
  }
  else {
    opserr << "WARNING SubdomainFE::addD_ForceSensitivity() - ";
    opserr << "- this should not be called on a Subdomain!\n";
  }
#endif
}





int
SubdomainFE::commitSensitivity(int gradNum, int numGrads)
{
  myEle->commitSensitivity(gradNum, numGrads);
  return 0;
}

// AddingSensitivity:END ////////////////////////////////////


int
SubdomainFE::updateElement()
{
  if (myEle != nullptr)
    return myEle->update();
  return 0;
}
