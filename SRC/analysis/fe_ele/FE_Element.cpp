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
// Purpose: This file contains the code for implementing the methods
// of the FE_Element class interface.
//
// Written: fmk
// Created: 11/96
// Revision: A
//
#include <FE_Element.h>
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
Matrix **FE_Element::theMatrices; // pointers to class wide matrices
Vector **FE_Element::theVectors;  // pointers to class widde vectors
int FE_Element::numFEs(0);           // number of objects


// constructor that take the corresponding model element.
FE_Element::FE_Element(int tag, Element *ele)
  : TaggedObject(tag)
  , myDOF_Groups((ele->getExternalNodes()).Size())
  , myID(ele->getNumDOF()),
   numDOF(ele->getNumDOF()), theModel(0), myEle(ele),
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

  // if this is the first FE_Element we now
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

  if (ele->isSubdomain() == false) {

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

    } else {
      // create matrices and vectors for each object instance
      theResidual = new Vector(numDOF);
      theTangent  = new Matrix(numDOF, numDOF);
    }

  } else {
    // as subdomains have own matrix for tangent and residual don't need
    // to set matrix and vector pointers to these objects
    theResidual = new Vector(numDOF);
      // invoke setFE_ElementPtr() method on Subdomain
    Subdomain *theSub = (Subdomain *)ele;
    theSub->setFE_ElementPtr(this);
  }

  // increment number of FE_Elements by 1
  numFEs++;
}


FE_Element::FE_Element(int tag, int numDOF_Group, int ndof)
  :TaggedObject(tag),
   myDOF_Groups(numDOF_Group), myID(ndof), numDOF(ndof), theModel(nullptr),
   myEle(nullptr), theResidual(nullptr), theTangent(nullptr), theIntegrator(nullptr)
{
  // this is for a subtype, the subtype must set the myDOF_Groups ID array
  numFEs++;

  // if this is the first FE_Element we now
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

  // as subtypes have no access to the tangent or residual we don't set them
  // this way we can detect if subclass does not provide all methods it should
}



FE_Element::~FE_Element()
{
  // decrement number of FE_Elements
  numFEs--;

  // delete tangent and residual if created specially
  if (numDOF > MAX_NUM_DOF) {
      if (theTangent != nullptr)
        delete theTangent;
      if (theResidual != nullptr) 
        delete theResidual;
  }

  // if this is the last FE_Element, clean up the
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
FE_Element::getDOFtags() const
{
  return myDOF_Groups;
}


const ID &
FE_Element::getID() const
{
  return myID;
}

void
FE_Element::setAnalysisModel(AnalysisModel &theAnalysisModel)
{
  theModel = &theAnalysisModel;
}

// void setID(int index, int value);
//        Method to set the corresponding index of the ID to value.

int
FE_Element::setID()
{
  int current = 0;

  assert(theModel != nullptr);

  int numGrps = myDOF_Groups.Size();
  for (int i=0; i<numGrps; i++) {
    int tag = myDOF_Groups(i);

    DOF_Group *dofPtr = theModel->getDOF_GroupPtr(tag);
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
FE_Element::getTangent(Integrator *theNewIntegrator)
{
  theIntegrator = theNewIntegrator;

  assert(myEle != nullptr);

  if (myEle->isSubdomain() == false) {
    if (theNewIntegrator != nullptr)
      theNewIntegrator->formEleTangent(this);

    return *theTangent;

  } else {
    Subdomain *theSub = (Subdomain *)myEle;
    theSub->computeTang();
    return theSub->getTang();
  }
}


void
FE_Element::zeroTangent()
{
  assert(myEle != nullptr);
  assert(myEle->isSubdomain() == false);
  theTangent->Zero();
}


void
FE_Element::addKtToTang(double fact)
{
  assert (myEle != nullptr);
  assert(myEle->isSubdomain() == false);

  // check for a quick return
  if (fact == 0.0)
      return;
  else
      theTangent->addMatrix(myEle->getTangentStiff(),fact);
}


void
FE_Element::addCtoTang(double fact)
{
  assert (myEle != nullptr);
  assert(myEle->isSubdomain() == false);

  // check for a quick return
  if (fact == 0.0)
    return;
  else
    theTangent->addMatrix(myEle->getDamp(),fact);
}


void
FE_Element::addMtoTang(double fact)
{
  if (myEle != nullptr) {
    assert(myEle->isSubdomain() == false);

    // check for a quick return
    if (fact == 0.0)
      return;
    else
      theTangent->addMatrix(myEle->getMass(),fact);
  }
}


void
FE_Element::addKiToTang(double fact)
{
  if (myEle != nullptr) {
    assert(myEle->isSubdomain() == false);

    // check for a quick return
    if (fact == 0.0)
      return;

    else // if (myEle->isSubdomain() == false)
      theTangent->addMatrix(myEle->getInitialStiff(), fact);
  }
}


void
FE_Element::addKpToTang(double fact, int numP)
{
  if (myEle != nullptr) {
    // check for a quick return
    if (fact == 0.0)
      return;

    else if (myEle->isSubdomain() == false) {
      const Matrix *thePrevMat = myEle->getPreviousK(numP);
      if (thePrevMat != nullptr)
        theTangent->addMatrix(*thePrevMat, fact);

    } else {
      opserr << "WARNING FE_Element::addKpToTang() - ";
      opserr << "- this should not be called on a Subdomain!\n";
    }
  }
}

int
FE_Element::storePreviousK(int numP)
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
FE_Element::getResidual(Integrator *theNewIntegrator)
{
  theIntegrator = theNewIntegrator;

  if (theIntegrator == nullptr)
    return *theResidual;

  assert(myEle != nullptr);

  if (myEle->isSubdomain() == false) {
    theNewIntegrator->formEleResidual(this);
    return *theResidual;

  } else {
    Subdomain *theSub = (Subdomain *)myEle;
    theSub->computeResidual();
    return theSub->getResistingForce();
  }
}

void
FE_Element::zeroResidual()
{
  assert(myEle != nullptr);
  assert(myEle->isSubdomain() == false);

  theResidual->Zero();
}


void
FE_Element::addRtoResidual(double fact)
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
FE_Element::addRIncInertiaToResidual(double fact)
{
  assert(myEle != nullptr);
  assert(myEle->isSubdomain() == false);

  // check for a quick return
  if (fact == 0.0)
      return;

  else {
    const Vector &eleResisting = myEle->getResistingForceIncInertia();
    theResidual->addVector(1.0, eleResisting, -fact);
  }
}


const Vector &
FE_Element::getTangForce(const Vector &disp, double fact)
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

  if (myEle->isSubdomain() == false) {
    // form the tangent again and then add the force
    theIntegrator->formEleTangent(this);
    theResidual->addMatrixVector(*theTangent,tmp,fact);

  } else {
    theResidual->addMatrixVector(((Subdomain *)myEle)->getTang(),tmp,fact);
  }
  return *theResidual;
}



const Vector &
FE_Element::getK_Force(const Vector &disp, double fact)
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
FE_Element::getKi_Force(const Vector &disp, double fact)
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
FE_Element::getM_Force(const Vector &disp, double fact)
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
FE_Element::getC_Force(const Vector &disp, double fact)
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
FE_Element::getLastIntegrator()
{
  return theIntegrator;
}


const Vector &
FE_Element::getLastResponse()
{
  assert(myEle != nullptr);

  if (theIntegrator != nullptr) {
    if (theIntegrator->getLastResponse(*theResidual,myID) < 0) {
      opserr << "WARNING FE_Element::getLastResponse()";
      opserr << " - the Integrator had problems with getLastResponse()\n";
    }
  }
  else {
    theResidual->Zero();
    opserr << "WARNING  FE_Element::getLastResponse()";
    opserr << " No Integrator yet passed\n";
  }

  Vector &result = *theResidual;
  return result;
}

void
FE_Element::addM_Force(const Vector &accel, double fact)
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
FE_Element::addD_Force(const Vector &accel, double fact)
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
FE_Element::addK_Force(const Vector &disp, double fact)
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
FE_Element::addLocalM_Force(const Vector &accel, double fact)
{
  assert(myEle != nullptr);
  assert(myEle->isSubdomain() == false);

  // check for a quick return
  if (fact == 0.0)
    return;

  theResidual->addMatrixVector(myEle->getMass(), accel, fact);
}

void
FE_Element::addLocalD_Force(const Vector &accel, double fact)
{
  assert(myEle != nullptr);
  assert(myEle->isSubdomain() == false);

  // check for a quick return
  if (fact == 0.0)
      return;

  if (theResidual->addMatrixVector(myEle->getDamp(), accel, fact) < 0){
    opserr << "WARNING FE_Element::addLocalD_Force() - ";
    opserr << "- addMatrixVector returned error\n";
  }
}


Element *
FE_Element::getElement()
{
  return myEle;
}


// AddingSensitivity:BEGIN /////////////////////////////////
void
FE_Element::addResistingForceSensitivity(int gradNumber, double fact)
{
  theResidual->addVector(1.0, myEle->getResistingForceSensitivity(gradNumber), -fact);
}

void
FE_Element::addM_ForceSensitivity(int gradNumber, const Vector &vect, double fact)
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
    opserr << "WARNING FE_Element::addM_ForceSensitivity() - ";
    opserr << "- addMatrixVector returned error\n";
  }
}

void
FE_Element::addD_ForceSensitivity(int gradNumber, const Vector &vect, double fact)
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
      if (theResidual->addMatrixVector(myEle->getDampSensitivity(gradNumber), tmp, fact) < 0){
        opserr << "WARNING FE_Element::addD_ForceSensitivity() - ";
        opserr << "- addMatrixVector returned error\n";
      }
    }
    else {
      opserr << "WARNING FE_Element::addD_ForceSensitivity() - ";
      opserr << "- this should not be called on a Subdomain!\n";
    }
}

void
FE_Element::addLocalD_ForceSensitivity(int gradNumber, const Vector &accel, double fact)
{
    if (myEle != nullptr) {
      // check for a quick return
      if (fact == 0.0)
          return;
      if (myEle->isSubdomain() == false) {
          if (theResidual->addMatrixVector(myEle->getDampSensitivity(gradNumber),
                                            accel, fact) < 0){

            opserr << "WARNING FE_Element::addLocalD_ForceSensitivity() - ";
            opserr << "- addMatrixVector returned error\n";
          }
      }
      else {
          opserr << "WARNING FE_Element::addLocalD_ForceSensitivity() - ";
          opserr << "- this should not be called on a Subdomain!\n";
      }
    }
    else {
        opserr << "WARNING FE_Element::addLocalD_ForceSensitivity() - no Element *given ";
        opserr << "- subclasses must provide implementation\n";
    }
}


void
FE_Element::addLocalM_ForceSensitivity(int gradNumber, const Vector &accel, double fact)
{
  assert(myEle != nullptr);
  assert(myEle->isSubdomain() == false);

  // check for a quick return
  if (fact == 0.0)
    return;

  if (theResidual->addMatrixVector(myEle->getMassSensitivity(gradNumber), accel, fact) < 0) {
    opserr << "WARNING FE_Element::addLocalD_ForceSensitivity() - ";
    opserr << "- addMatrixVector returned error\n";
  }
}





int
FE_Element::commitSensitivity(int gradNum, int numGrads)
{
  myEle->commitSensitivity(gradNum, numGrads);
  return 0;
}

// AddingSensitivity:END ////////////////////////////////////


int
FE_Element::updateElement()
{
  if (myEle != nullptr)
    return myEle->update();
  return 0;
}
