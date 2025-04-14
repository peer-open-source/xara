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
// This file contains the implementation of the Node class
//
// Written: fmk
// Created: 11/96
// Revision: A
//
#include <Node.h>
#include <NodeData.h>
#include <stdlib.h>
#include <assert.h>
#include <Element.h>
#include <Vector.h>
#include <Versor.h>
#include <Matrix.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <DOF_Group.h>
#include <string.h>
#include <Information.h>
#include <Parameter.h>

// AddingSensitivity:BEGIN //////////////////////////
#include <Domain.h>
#include <Element.h>
#include <ElementIter.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
// AddingSensitivity:END ////////////////////////////

#include <NodalLoad.h>

#include <OPS_Globals.h>

using namespace OpenSees;
Matrix **Node::theMatrices = nullptr;
int Node::numMatrices = 0;


// for FEM_Object Broker to use
Node::Node(int theClassTag)
:TaggedObject(0),MovableObject(theClassTag),
 numberDOF(0), theDOF_GroupPtr(0),
 Crd(0), commitDisp(0), commitVel(0), commitAccel(0),
 trialDisp(0), trialVel(0), trialAccel(0), unbalLoad(0), incrDisp(0),
 incrDeltaDisp(0),
 disp(0), vel(0), accel(0), dbTag1(0), 
 rotation(nullptr),
 dbTag2(0), dbTag3(0), dbTag4(0),
 R(0), mass(0), unbalLoadWithInertia(0), alphaM(0.0), theEigenvectors(0),
 index(-1), reaction(0)//, displayLocation(0)
{
  // for FEM_ObjectBroker, recvSelf() must be invoked on object

  // Sensitivity
  dispSensitivity = nullptr;
  velSensitivity  = nullptr;
  accSensitivity  = nullptr;
  parameterID     = 0;


  theNodalThermalActionPtr = 0; //Added by Liming for initializing NodalLoadPointer, [SIF]
}


Node::Node(int tag, int theClassTag)
:TaggedObject(tag), MovableObject(theClassTag),
 numberDOF(0), theDOF_GroupPtr(0),
 Crd(0), commitDisp(0), commitVel(0), commitAccel(0),
 trialDisp(0), trialVel(0), trialAccel(0), unbalLoad(0), incrDisp(0),
 incrDeltaDisp(0),
 disp(0), vel(0), accel(0), 
 rotation(nullptr),
 dbTag1(0), dbTag2(0), dbTag3(0), dbTag4(0),
  R(0), mass(0), unbalLoadWithInertia(0), alphaM(0.0), theEigenvectors(0),
 index(-1), reaction(0)//, displayLocation(0)
{
  // for subclasses - they must implement all the methods with
  // their own data structures.

  dispSensitivity = nullptr;
  velSensitivity  = nullptr;
  accSensitivity  = nullptr;
  parameterID = 0;

  theNodalThermalActionPtr = 0;//Added by Liming for initializing NodalLoadPointer, [SIF]
}

Node::Node(int tag, int ndof, double Crd1, Vector *dLoc)
:TaggedObject(tag),MovableObject(NOD_TAG_Node),
 numberDOF(ndof), theDOF_GroupPtr(0),
 Crd(0), commitDisp(0), commitVel(0), commitAccel(0),
 trialDisp(0), trialVel(0), trialAccel(0), unbalLoad(0), incrDisp(0),
 incrDeltaDisp(0),
 disp(0), vel(0), accel(0), 
 rotation(nullptr),
 dbTag1(0), dbTag2(0), dbTag3(0), dbTag4(0),
 R(0), mass(0), unbalLoadWithInertia(0), alphaM(0.0), theEigenvectors(0),
 index(-1), reaction(0)//, displayLocation(0)
{
  this->createDisp();
  dispSensitivity = 0;
  velSensitivity  = 0;
  accSensitivity  = 0;
  parameterID     = 0;

  theNodalThermalActionPtr = 0;//Added by Liming for initializing NodalLoadPointer, [SIF]

  Crd = new Vector(1);
  (*Crd)(0) = Crd1;
}


//  Node(int tag, int ndof, double Crd1, double yCrd);
//      constructor for 2d nodes
Node::Node(int tag, int ndof, double Crd1, double Crd2, Vector *dLoc)
:TaggedObject(tag), MovableObject(NOD_TAG_Node),
 numberDOF(ndof), theDOF_GroupPtr(0),
 Crd(0), commitDisp(0), commitVel(0), commitAccel(0),
 trialDisp(0), trialVel(0), trialAccel(0), unbalLoad(0), incrDisp(0),
 incrDeltaDisp(0),
 disp(0), vel(0), accel(0), 
 rotation(nullptr),
 dbTag1(0), dbTag2(0), dbTag3(0), dbTag4(0),
 R(0), mass(0), unbalLoadWithInertia(0), alphaM(0.0), theEigenvectors(0),
 reaction(0)//, displayLocation(0)
{
  this->createDisp();
  dispSensitivity = 0;
  velSensitivity = 0;
  accSensitivity = 0;
  parameterID = 0;

  theNodalThermalActionPtr = 0;//Added by Liming for initializing NodalLoadPointer, [SIF]

  Crd = new Vector(2);
  (*Crd)(0) = Crd1;
  (*Crd)(1) = Crd2;

  index = -1;
}


//  Node(int tag, int ndof, double Crd1, double Crd2, double zCrd);
//      constructor for 3d nodes

 Node::Node(int tag, int ndof, double Crd1, double Crd2, double Crd3, Vector *dLoc)
:TaggedObject(tag), MovableObject(NOD_TAG_Node),
 numberDOF(ndof), theDOF_GroupPtr(0),
 Crd(0), commitDisp(0), commitVel(0), commitAccel(0),
 trialDisp(0), trialVel(0), trialAccel(0), unbalLoad(0), incrDisp(0),
 incrDeltaDisp(0),
 disp(0), vel(0), accel(0), 
 rotation(nullptr),
 dbTag1(0), dbTag2(0), dbTag3(0), dbTag4(0),
 R(0), mass(0), unbalLoadWithInertia(0), alphaM(0.0), theEigenvectors(0),
 reaction(0)//, displayLocation(0)
{
  this->createDisp();
  dispSensitivity = 0;
  velSensitivity = 0;
  accSensitivity = 0;
  parameterID = 0;

  theNodalThermalActionPtr = 0;//Added by Liming for initializing NodalLoadPointer, [SIF]

  Crd = new Vector(3);
  (*Crd)(0) = Crd1;
  (*Crd)(1) = Crd2;
  (*Crd)(2) = Crd3;

  index = -1;
}


// used for domain decomposition & external nodes
//  copy everything but the mass
//  we should really set the mass to 0.0
Node::Node(const Node &otherNode, bool copyMass)
  :TaggedObject(otherNode.getTag()),MovableObject(otherNode.getClassTag()),
 numberDOF(otherNode.numberDOF), theDOF_GroupPtr(0),
 Crd(0), commitDisp(0), commitVel(0), commitAccel(0),
 trialDisp(0), trialVel(0), trialAccel(0), unbalLoad(0), incrDisp(0),
 incrDeltaDisp(0),
 disp(0), vel(0), accel(0), 
 rotation(nullptr),
 dbTag1(0), dbTag2(0), dbTag3(0), dbTag4(0),
 R(0), mass(0), unbalLoadWithInertia(0), alphaM(0.0), theEigenvectors(0),
   reaction(0)//, displayLocation(0)
{
  this->createDisp();
  dispSensitivity = 0;
  velSensitivity = 0;
  accSensitivity = 0;
  parameterID = 0;

  theNodalThermalActionPtr = 0;//Added by Liming for initializing NodalLoadPointer, [SIF]

  Crd = new Vector(otherNode.getCrds());

  if (otherNode.commitDisp != nullptr) {
    // this->createDisp();
    *commitDisp = *otherNode.commitDisp;
//  for (int i=0; i<4*numberDOF; i++)
//    disp[i] = otherNode.disp[i];
  }

  if (otherNode.commitVel != nullptr) {
    this->createVel();
    for (int i=0; i<2*numberDOF; i++)
      vel[i] = otherNode.vel[i];
  }

  if (otherNode.commitAccel != nullptr) {
    this->createAccel();
    for (int i=0; i<2*numberDOF; i++)
      accel[i] = otherNode.accel[i];
  }


  if (otherNode.unbalLoad != nullptr){
    unbalLoad = new Vector(*(otherNode.unbalLoad));
    unbalLoad->Zero();
  }

  if (otherNode.mass != 0 && copyMass == true) {
    mass = new Matrix(*(otherNode.mass)) ;
    // TODO; zero?
  }

  if (otherNode.R != 0) {
    R = new Matrix(*(otherNode.R));
  }

  index = -1;
}


// ~Node():
//       destructor

Node::~Node()
{
  // delete anything that we created with new
  if (Crd != 0)
    delete Crd;

  if (commitDisp != 0)
    delete commitDisp;

  if (commitVel != 0)
    delete commitVel;

  if (commitAccel != 0)
    delete commitAccel;

  if (trialDisp != 0)
    delete trialDisp;

  if (trialVel != 0)
    delete trialVel;

  if (trialAccel != 0)
    delete trialAccel;

  if (rotation != nullptr)
    delete rotation;

  if (incrDisp != 0)
    delete incrDisp;

  if (incrDeltaDisp != 0)
    delete incrDeltaDisp;

  if (unbalLoad != 0)
    delete unbalLoad;

  if (disp != 0)
    delete [] disp;

  if (vel != 0)
    delete [] vel;

  if (accel != 0)
    delete [] accel;

  if (mass != 0)
    delete mass;

  if (R != 0)
    delete R;

  if (unbalLoadWithInertia != 0)
    delete unbalLoadWithInertia;

  if (theEigenvectors != 0)
    delete theEigenvectors;

  // AddingSensitivity:BEGIN ///////////////////////////////////////
  if (dispSensitivity != 0)
    delete dispSensitivity;
  if (velSensitivity != 0)
    delete velSensitivity;
  if (accSensitivity != 0)
    delete accSensitivity;
  // AddingSensitivity:END /////////////////////////////////////////

  if (reaction != 0)
    delete reaction;

  if (theDOF_GroupPtr != 0)
    theDOF_GroupPtr->resetNodePtr();
}


int
Node::getNumberDOF(void) const
{
  return  numberDOF;
}


void
Node::setDOF_GroupPtr(DOF_Group *theDOF_Grp)
{
  // set the DOF_Group pointer
  theDOF_GroupPtr = theDOF_Grp;
}


DOF_Group *
Node::getDOF_GroupPtr(void)
{
  return theDOF_GroupPtr;
}


const Vector &
Node::getCrds() const
{
  // return the vector of nodal coordinates
  return *Crd;
}




const Vector &
Node::getDisp()
{
  // construct memory and Vectors for trial and committed
  // displacement on first call to this method, getTrialDisp()
  // setTrialDisp() or incrTrialDisp()
#if 0
  if (commitDisp == nullptr)
    this->createDisp();
#endif
  // return the committed disp
  return *commitDisp;
}

const Vector &
Node::getVel()
{
  // construct memory and Vectors for trial and committed
  // velocity on first call to this method, getTrialVel()
  // setTrialVel() or incrTrialVel()
  if (commitVel == nullptr)
    this->createVel();

  // return the velocity
  return *commitVel;
}


const Vector &
Node::getAccel()
{
  // construct memory and Vectors for trial and committed
  // accel on first call to this method, getTrialAccel()
  // setTrialAccel() or incrTrialAccel()
  if (commitAccel == nullptr)
      this->createAccel();

  return *commitAccel;
}


/* *********************************************************************
**
**   Methods to return the trial response quantities similar to committed
**
** *********************************************************************/

const Vector &
Node::getTrialDisp()
{
#if 0
    if (trialDisp == nullptr)
      this->createDisp();
#endif
    return *trialDisp;
}

const Vector &
Node::getIncrDisp(void)
{
#if 0
    if (incrDisp == nullptr)
      this->createDisp();
#endif
    return *incrDisp;
}

const Vector &
Node::getIncrDeltaDisp()
{
#if 0
    if (incrDeltaDisp == 0)
      this->createDisp();
#endif
    return *incrDeltaDisp;
}

Versor
Node::getTrialRotation()
{
  if (rotation == nullptr) {
    if (this->getNumberDOF() < 6)
      return Versor();
    else
      rotation = new Versor{{0.0, 0.0, 0.0}, 1.0};
  }

  return *rotation;
}


const Vector &
Node::getTrialVel()
{
    if (trialVel == nullptr)
      this->createVel();

    return *trialVel;
}



const Vector &
Node::getTrialAccel(void)
{
    if (trialAccel == nullptr)
      this->createAccel();

    return *trialAccel;
}


int
Node::setTrialDisp(double value, int dof)
{
  // check vector arg is of correct size
  assert(dof >= 0 && dof < numberDOF);
#if 0
  if (trialDisp == nullptr)
    this->createDisp();
#endif
  // perform the assignment .. we don't go through Vector interface
  // as we are sure of size and this way is quicker
  double tDisp = value;
  disp[dof+2*numberDOF] = tDisp - disp[dof+numberDOF];
  disp[dof+3*numberDOF] = tDisp - disp[dof];
  disp[dof]             = tDisp;
  return 0;
}

int
Node::setTrialDisp(const Vector &newTrialDisp)
{
  // check vector arg is of correct size
  assert(newTrialDisp.Size() == numberDOF);

#if 0
  if (trialDisp == 0)
      this->createDisp();
#endif

  // perform the assignment .. we don't go through Vector interface
  // as we are sure of size and this way is quicker
  for (int i=0; i<numberDOF; i++) {
      double tDisp = newTrialDisp(i);
      disp[i+2*numberDOF] = tDisp - disp[i+numberDOF];
      disp[i+3*numberDOF] = tDisp - disp[i];
      disp[i] = tDisp;
  }

  return 0;
}

int
Node::setTrialVel(const Vector &newTrialVel)
{
  // check vector arg is of correct size
  assert(newTrialVel.Size() == numberDOF);

  // construct memory and Vectors for trial and committed
  // accel on first call to this method, getTrialVEl(),
  // getVEl(), or incrTrialVel()
  if (trialVel == 0)
    this->createVel();

  // set the trial quantities
  for (int i=0; i<numberDOF; i++)
    vel[i] = newTrialVel(i);

  return 0;
}

int
Node::incrTrialDisp(const Vector &incrDispl)
{
  // check vector arg is of correct size
  assert(incrDispl.Size() == numberDOF);

  // create a copy if no trial exists and add committed
  if (trialDisp == nullptr) {
    this->createDisp();
    for (int i = 0; i<numberDOF; i++) {
      double incrDispI = incrDispl(i);
      disp[i]             = incrDispI;
      disp[i+2*numberDOF] = incrDispI;
      disp[i+3*numberDOF] = incrDispI;
    }
    return 0;
  }

  // otherwise set trial = incr + trial
  for (int i = 0; i<numberDOF; i++) {
      double incrDispI = incrDispl(i);
      disp[i]             += incrDispI;
      disp[i+2*numberDOF] += incrDispI;
      disp[i+3*numberDOF]  = incrDispI;
  }
  if (rotation != nullptr && this->getNumberDOF() >= 6) {
    (*rotation) = (*rotation) * Versor::from_vector(&disp[3+3*numberDOF]);
  }

  return 0;
}


int
Node::setTrialAccel(const Vector &newTrialAccel)
{
    // check vector arg is of correct size
    assert(newTrialAccel.Size() == numberDOF);

    // create a copy if no trial exists
    if (trialAccel == 0)
      this->createAccel();

    // use vector assignment otherwise
    for (int i=0; i<numberDOF; i++)
      accel[i] = newTrialAccel(i);

    return 0;
}



int
Node::incrTrialVel(const Vector &incrVel)
{
  // check vector arg is of correct size
  assert(incrVel.Size() == numberDOF);

  // create Vectors and array if none exist and set trial
  if (trialVel == nullptr) {
    this->createVel();
    for (int i = 0; i<numberDOF; i++)
        vel[i] = incrVel(i);
    return 0;
  }

  // otherwise set trial = incr + trial
  for (int i = 0; i<numberDOF; i++)
    vel[i] += incrVel(i);

  return 0;
}


int
Node::incrTrialAccel(const Vector &incrAccel)
{
  // check vector arg is of correct size
  assert(incrAccel.Size() == numberDOF);

  // create a copy if no trial exists and add committed
  if (trialAccel == 0) {
      this->createAccel();
    for (int i = 0; i<numberDOF; i++)
          accel[i] = incrAccel(i);

      return 0;
  }

  // otherwise set trial = incr + trial
  for (int i = 0; i<numberDOF; i++)
    accel[i] += incrAccel(i);

  return 0;
}


void
Node::zeroUnbalancedLoad(void)
{
  if (unbalLoad != nullptr)
    unbalLoad->Zero();
}

int
Node::addUnbalancedLoad(const Vector &add, double fact)
{
  // check vector arg is of correct size
  if (add.Size() != numberDOF) {
    opserr << "Node::addunbalLoad - load to add of incorrect size ";
    opserr << add.Size() << " should be " <<  numberDOF << endln;
    return -1;
  }
  // if no load yet create it and assign
  if (unbalLoad == 0) {
    unbalLoad = new Vector(add);
    if (fact != 1.0)
          (*unbalLoad) *= fact;
    return 0;
  }

  // add fact*add to the unbalanced load
  unbalLoad->addVector(1.0, add,fact);
  return 0;
}



int
Node::addInertiaLoadToUnbalance(const Vector &accelG, double fact)
{
  // simply return if node has no mass or R matrix
  if (mass == 0 || R == 0)
    return 0;

  // otherwise we must determine MR accelG
  assert(accelG.Size() == R->noCols());

  // if no load yet create it and assign
  if (unbalLoad == nullptr)
      unbalLoad = new Vector(numberDOF);

  // form - fact * M*R*accelG and add it to the unbalanced load
  //(*unbalLoad) -= ((*mass) * (*R) * accelG)*fact;

  Matrix MR(mass->noRows(), R->noCols());
  MR.addMatrixProduct(0.0, *mass, *R, 1.0);
  unbalLoad->addMatrixVector(1.0, MR, accelG, -fact);

  return 0;
}



int
Node::addInertiaLoadSensitivityToUnbalance(const Vector &accelG, double fact, bool somethingRandomInMotions)
{
  // simply return if node has no mass or R matrix
  if (mass == 0 || R == 0)
    return 0;

  // otherwise we must determine MR accelG
  assert(accelG.Size() == R->noCols());

  // if no load yet create it and assign
  if (unbalLoad == nullptr)
      unbalLoad = new Vector(numberDOF);

  // form - fact * M*R*accelG and add it to the unbalanced load
  //(*unbalLoad) -= ((*mass) * (*R) * accelG)*fact;


  Matrix massSens(mass->noRows(),mass->noCols());
  massSens = this->getMassSensitivity();

  Matrix MR(mass->noRows(), R->noCols());

  if (somethingRandomInMotions) {
    MR.addMatrixProduct(0.0, *mass, *R, 1.0);
  }
  else {
    MR.addMatrixProduct(0.0, massSens, *R, 1.0);
  }
  unbalLoad->addMatrixVector(1.0, MR, accelG, -fact);

  return 0;
}



const Vector &
Node::getUnbalancedLoad()
{
  // make sure it was created before we return it
  if (unbalLoad == nullptr)
    unbalLoad = new Vector(numberDOF);

  // return the unbalanced load
  return *unbalLoad;
}



const Vector &
Node::getUnbalancedLoadIncInertia()
{
    // make sure it was created before we return it
    if (unbalLoadWithInertia == nullptr) {
      unbalLoadWithInertia = new Vector(this->getUnbalancedLoad());

    } else
      (*unbalLoadWithInertia) = this->getUnbalancedLoad();

    if (mass != nullptr) {

      const Vector &theAccel = this->getTrialAccel(); // in case accel not created
      unbalLoadWithInertia->addMatrixVector(1.0, *mass, theAccel, -1.0);

      if (alphaM != 0.0) {
      const Vector &theVel = this->getTrialVel(); // in case vel not created
      unbalLoadWithInertia->addMatrixVector(1.0, *mass, theVel, -alphaM);
      }
    }

    return *unbalLoadWithInertia;
}

int
Node::commitState()
{
  // check disp exists, if does set commit = trial, incr = 0.0
  if (trialDisp != 0) {
    for (int i=0; i<numberDOF; i++) {
    disp[i+numberDOF] = disp[i];
      disp[i+2*numberDOF] = 0.0;
      disp[i+3*numberDOF] = 0.0;
    }
  }

  // Check vel exists, if does set commit = trial
  if (trialVel != 0) {
    for (int i=0; i<numberDOF; i++)
    vel[i+numberDOF] = vel[i];
  }

  // Check accel exists, if does set commit = trial
  if (trialAccel != 0) {
    for (int i=0; i<numberDOF; i++)
    accel[i+numberDOF] = accel[i];
  }

  return 0;
}



int
Node::revertToLastCommit()
{
  // check memory was allocated, if it does set trial = last commit, incr = 0
  if (disp != 0) {
    for (int i=0 ; i<numberDOF; i++) {
    disp[i] = disp[i+numberDOF];
    disp[i+2*numberDOF] = 0.0;
    disp[i+3*numberDOF] = 0.0;
    }
  }

  // check vel exists, if does set trial = last commit
  if (vel != 0) {
    for (int i=0 ; i<numberDOF; i++)
    vel[i] = vel[numberDOF+i];
  }

  // check accel exists, if does set trial = last commit
  if (accel != 0) {
    for (int i=0 ; i<numberDOF; i++)
    accel[i] = accel[numberDOF+i];
  }

  // if we get here we are done
  return 0;
}


int
Node::revertToStart()
{
  // check disp exists, if does set all to zero
  if (disp != 0) {
    for (int i=0 ; i<4*numberDOF; i++)
    disp[i] = 0.0;
  }

  // check vel exists, if does set all to zero
  if (vel != 0) {
    for (int i=0 ; i<2*numberDOF; i++)
    vel[i] = 0.0;
  }

  // check accel exists, if does set all to zero
  if (accel != 0) {
    for (int i=0 ; i<2*numberDOF; i++)
    accel[i] = 0.0;
  }

  if (unbalLoad != nullptr)
    (*unbalLoad) *= 0;


  if (dispSensitivity != 0)
    dispSensitivity->Zero();

  if (velSensitivity != 0)
    velSensitivity->Zero();

  if (accSensitivity != 0)
    accSensitivity->Zero();

  if (rotation != nullptr)
    *rotation = Versor{{0.0, 0.0, 0.0}, 1.0};

  return 0;
}


const Matrix &
Node::getMass()
{
    if (index == -1) {
      setGlobalMatrices();
    }

    // make sure it was created before we return it
    if (mass == 0) {
      theMatrices[index]->Zero();
      return *theMatrices[index];
    } else
      return *mass;
}


int
Node::setRayleighDampingFactor(double alpham) {
  alphaM = alpham;
  return 0;
}


const Matrix &
Node::getDamp()
{
    if (index == -1) {
      setGlobalMatrices();
    }

    // make sure it was created before we return it
    if (mass == 0 || alphaM == 0.0) {
      theMatrices[index]->Zero();
      return *theMatrices[index];
    } else {
      Matrix &result = *theMatrices[index];
      result = *mass;
      result *= alphaM;
      return result;
    }
}


const Matrix &
Node::getDampSensitivity(void)
{
    if (index == -1) {
      setGlobalMatrices();
    }

    // make sure it was created before we return it
    if (mass == 0 || alphaM == 0.0) {
      theMatrices[index]->Zero();
      return *theMatrices[index];
    } else {
      Matrix &result = *theMatrices[index];
        result.Zero();
      //result = *mass;
      //result *= alphaM;
      return result;
    }
}


int
Node::setMass(const Matrix &newMass)
{
    // check right size
    assert(newMass.noRows() == numberDOF && newMass.noCols() == numberDOF);


    // create a matrix if no mass yet set
    if (mass == nullptr) {
      mass = new Matrix(newMass);
      return 0;
    }

    // otherwise assign mass
    (*mass) = newMass;

    return 0;
}



int
Node::setNumColR(int numCol)
{
  if (R != 0) {
    if (R->noCols() != numCol) {
      delete R;
      R = new Matrix(numberDOF, numCol);
    }
  } else
    R = new Matrix(numberDOF, numCol);

  R->Zero();
  return 0;
}

int
Node::setR(int row, int col, double Value)
{
  // ensure R had been set
  if (R == nullptr) {
    opserr << "Node:setR() - R has not been initialised\n";
    return -1;
  }

  // ensure row, col in range (matrix assignment will catch this - extra work)
  if (row < 0 || row > numberDOF || col < 0 || col > R->noCols()) {
    opserr << "Node:setR() - row, col index out of range\n";
    return -1;
  }

  // do the assignment
  (*R)(row,col) = Value;

  /*
  // to test uniform excitation pattern with consistent mass matrices:
  // found that the static application of a unit ground displacement
  // needs to also be applied to the constrained DOFs
  Domain *theDomain = this->getDomain();
  SP_ConstraintIter &theSPs = theDomain->getSPs();
  SP_Constraint *theSP;
  // assign zero if there is a homogeneous SP
  while ((theSP = theSPs()) != 0) {
      if (theSP->getNodeTag() == this->getTag() &&
          theSP->getDOF_Number() == row &&
          theSP->isHomogeneous()) {
              (*R)(row,col) = 0.0;
      }
  }
  */

  return 0;
}



const Vector &
Node::getRV(const Vector &V)
{
    // we store the product of RV in unbalLoadWithInertia

    // make sure unbalLoadWithInertia was created, if not create it
    if (unbalLoadWithInertia == nullptr)
      unbalLoadWithInertia = new Vector(numberDOF);

    // see if quick return , i.e. R == 0
    if (R == 0) {
      unbalLoadWithInertia->Zero();
      return *unbalLoadWithInertia;
    }

    // check dimesions of R and V
    if (R->noCols() != V.Size()) {
      opserr << "WARNING Node::getRV() - R and V of incompatible dimesions\n";
      opserr << "R: " << *R << "V: " << V;
      unbalLoadWithInertia->Zero();
      return *unbalLoadWithInertia;
    }

    // determine the product
    unbalLoadWithInertia->addMatrixVector(0.0, *R, V, 1.0);
    return *unbalLoadWithInertia;
}


int
Node::setNumEigenvectors(int numVectorsToStore)
{
  // ensure a positive number of vectors
  if (numVectorsToStore <= 0) {
    opserr << "Node::setNumEigenvectors() - " << numVectorsToStore << " < 0\n";
    return -1;
  }

  // if matrix not yet assigned or not of correct size delete old and create new
  if (theEigenvectors == 0 || theEigenvectors->noCols() != numVectorsToStore) {
    if (theEigenvectors != 0)
      delete theEigenvectors;

    theEigenvectors = new Matrix(numberDOF, numVectorsToStore);

  } else {
    // zero the eigenvector matrix
    theEigenvectors->Zero();
  }

  return 0;
}

int
Node::setEigenvector(int mode, const Vector &eigenVector)
{
  if (theEigenvectors == 0 || theEigenvectors->noCols() < mode) {
    opserr << "Node::setEigenvectors() - mode " << mode << " invalid\n";
    return -1;
  }

  assert(eigenVector.Size() == numberDOF);

  // set the values
  for (int i=0; i<numberDOF; i++)
    (*theEigenvectors)(i, mode-1) = eigenVector(i);

  return 0;
}

const Matrix &
Node::getEigenvectors(void)
{
  // check the eigen vectors have been set
  if (theEigenvectors == 0) {
    opserr << "Node::getEigenvectors() - eigenvectors have not been set\n";
    // TODO: Handle this!
    exit(-1);
  }

  return *theEigenvectors;
}


int
Node::sendSelf(int cTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();

    ID data(14);
    data(0) = this->getTag();
    data(1) = numberDOF;

    // indicate whether vector quantaties have been formed
    if (disp == 0)       data(2) = 1; else data(2) = 0;
    if (vel == 0)        data(3) = 1; else data(3) = 0;
    if (accel == 0)      data(4) = 1; else data(4) = 0;
    if (mass == 0)       data(5) = 1; else data(5) = 0;
    if (unbalLoad  == 0) data(6) = 1; else data(6) = 0;
    if (R == 0)
      data(12) = 1;
    else {
      data(12) = 0;
      data(13) = R->noCols();
    }

    data(7) = Crd->Size();

    if (dbTag1 == 0)
      dbTag1 = theChannel.getDbTag();
    if (dbTag2 == 0)
      dbTag2 = theChannel.getDbTag();
    if (dbTag3 == 0)
      dbTag3 = theChannel.getDbTag();
    if (dbTag4 == 0)
      dbTag4 = theChannel.getDbTag();

    data(8) = dbTag1;
    data(9) = dbTag2;
    data(10) = dbTag3;
    data(11) = dbTag4;

    int res = 0;

    res = theChannel.sendID(dataTag, cTag, data);
    if (res < 0) {
      opserr << " Node::sendSelf() - failed to send ID data\n";
      return res;
    }

    res = theChannel.sendVector(dataTag, cTag, *Crd);
    if (res < 0) {
      opserr << " Node::sendSelf() - failed to send Vecor data\n";
      return res;
    }

    if (commitDisp != 0) {
      res = theChannel.sendVector(dbTag1, cTag, *commitDisp);
      if (res < 0) {
        opserr << " Node::sendSelf() - failed to send Disp data\n";
        return res;
      }
    }

    if (commitVel != 0) {
      res = theChannel.sendVector(dbTag2, cTag, *commitVel);
      if (res < 0) {
        opserr << " Node::sendSelf() - failed to send Vel data\n";
        return res;
      }
    }

    if (commitAccel != 0) {
      res = theChannel.sendVector(dbTag3, cTag, *commitAccel);
      if (res < 0) {
        opserr << " Node::sendSelf() - failed to send Accel data\n";
        return res;
      }
    }

    if (mass != 0) {
      res = theChannel.sendMatrix(dataTag, cTag, *mass);
      if (res < 0) {
        opserr << " Node::sendSelf() - failed to send Mass data\n";
        return res;
      }
    }

    if (R != nullptr) {
      res = theChannel.sendMatrix(dbTag2, cTag, *R);
      if (res < 0) {
        opserr << " Node::sendSelf() - failed to send R data\n";
        return res;
      }
    }

    if (unbalLoad  != 0) {
      res = theChannel.sendVector(dbTag4, cTag, *unbalLoad);
      if (res < 0) {
        opserr << " Node::sendSelf() - failed to send Load data\n";
        return res;
      }
    }

    // if get here successful
    return 0;
}

int
Node::recvSelf(int cTag, Channel &theChannel,
             FEM_ObjectBroker &theBroker)
{
    int res = 0;
    int dataTag = this->getDbTag();


    ID data(14);
    res = theChannel.recvID(dataTag, cTag, data);
    if (res < 0) {
      opserr << "Node::recvSelf() - failed to receive ID data\n";
      return res;
    }

    this->setTag(data(0));
    numberDOF = data(1);
    int numberCrd = data(7);

    dbTag1 = data(8);
    dbTag2 = data(9);
    dbTag3 = data(10);
    dbTag4 = data(11);

    // create a Vector to hold coordinates IF one needed
    if (Crd == nullptr)
      Crd = new Vector(numberCrd);

    if (theChannel.recvVector(dataTag, cTag, *Crd) < 0) {
      opserr << "Node::recvSelf() - failed to receive the Coordinate vector\n";
      return -2;
    }

    if (commitDisp == nullptr)
      this->createDisp();

    if (data(2) == 0) {
      // create the disp vectors if node is a total blank
//    if (commitDisp == 0)
//       this->createDisp();

      // recv the committed disp
      if (theChannel.recvVector(dbTag1, cTag, *commitDisp) < 0) {
      opserr << "Node::recvSelf - failed to receive Disp data\n";
      return res;
      }

      // set the trial quantities equal to committed
      for (int i=0; i<numberDOF; i++)
      disp[i] = disp[i+numberDOF];  // set trial equal committed

    } else if (commitDisp != nullptr) {
      // if going back to initial we will just zero the vectors
      commitDisp->Zero();
      trialDisp->Zero();
    }


    if (data(3) == 0) {
      // create the vel vectors if node is a total blank
      if (commitVel == nullptr)
      this->createVel();

      // recv the committed vel
      if (theChannel.recvVector(dbTag2, cTag, *commitVel) < 0) {
      opserr << "Node::recvSelf - failed to receive Velocity data\n";
      return -3;
      }

      // set the trial quantity
      for (int i=0; i<numberDOF; i++)
      vel[i] = vel[i+numberDOF];  // set trial equal committed
    }

    if (data(4) == 0) {
      // create the vel vectors if node is a total blank
      if (commitAccel == 0)
      this->createAccel();

      // recv the committed accel
      if (theChannel.recvVector(dbTag3, cTag, *commitAccel) < 0) {
      opserr << "Node::recvSelf - failed to receive Acceleration data\n";
      return -4;
      }

      // set the trial values
      for (int i=0; i<numberDOF; i++)
      accel[i] = accel[i+numberDOF];  // set trial equal committed
    }

    if (data(5) == 0) {
      // make some room and read in the vector
      if (mass == 0) {
      mass = new Matrix(numberDOF,numberDOF);
      }
      if (theChannel.recvMatrix(dataTag, cTag, *mass) < 0) {
      opserr << "Node::recvSelf() - failed to receive Mass data\n";
      return -6;
      }
    }

    if (data(12) == 0) {
      // create a matrix for R
      int noCols = data(13);
      if (R == nullptr) {
        R = new Matrix(numberDOF, noCols);
      }
      // now recv the R matrix
      if (theChannel.recvMatrix(dbTag2, cTag, *R) < 0) {
        opserr << "Node::recvSelf() - failed to receive R data\n";
        return res;
      }
    }


    if (data(6) == 0) {
      // create a vector for the load
      if (unbalLoad == 0) {
      unbalLoad = new Vector(numberDOF);
      if (unbalLoad == 0) {
        opserr << "Node::recvData -- ran out of memory\n";
        return -10;
      }
      }
      if (theChannel.recvVector(dbTag4, cTag, *unbalLoad) < 0) {
      opserr << "Node::recvSelf() - failed to receive Load data\n";
      return res;
      }
    }


  index = -1;
  if (numMatrices != 0) {
    for (int i=0; i<numMatrices; i++)
      if (theMatrices[i]->noRows() == numberDOF) {
      index = i;
      i = numMatrices;
      }
  }
  if (index == -1) {
    Matrix **nextMatrices = new Matrix *[numMatrices+1];

    for (int j=0; j<numMatrices; j++)
      nextMatrices[j] = theMatrices[j];

    nextMatrices[numMatrices] = new Matrix(numberDOF, numberDOF);

    if (numMatrices != 0)
      delete [] theMatrices;
    index = numMatrices;
    numMatrices++;
    theMatrices = nextMatrices;
  }


  return 0;
}


// createDisp(), createVel() and createAccel():
// private methods to create the arrays to hold the disp, vel and acceleration
// values and the Vector objects for the committed and trial quantaties.

int
Node::createDisp(void)
{
  // trial , committed, incr = (committed-trial)
  // Use {} to allocate zero-initialized space for the data
  disp          = new double[4*numberDOF]{};
  trialDisp     = new Vector(disp, numberDOF);
  commitDisp    = new Vector(&disp[numberDOF], numberDOF);
  incrDisp      = new Vector(&disp[2*numberDOF], numberDOF);
  incrDeltaDisp = new Vector(&disp[3*numberDOF], numberDOF);
  return 0;
}


int
Node::createVel(void)
{
  // Use {} to allocate zero-initialized space for the data
  vel       = new double[2*numberDOF]{};
  commitVel = new Vector(&vel[numberDOF], numberDOF);
  trialVel  = new Vector(vel, numberDOF);
  return 0;
}

int
Node::createAccel(void)
{
  // Use {} to allocate zero-initialized space for the data
  accel       = new double[2*numberDOF]{};
  commitAccel = new Vector(&accel[numberDOF], numberDOF);
  trialAccel  = new Vector(accel, numberDOF);
  return 0;
}

void
Node::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_CURRENTSTATE) { // print out everything
        s << "\n  Node: " << this->getTag() << endln;
        s << "\tCoordinates  : " << *Crd;
        if (commitDisp != 0)
            s << "\tDisps: " << *trialDisp;
        if (commitVel != 0)
            s << "\tVelocities   : " << *trialVel;
        if (commitAccel != 0)
            s << "\tcommitAccels: " << *trialAccel;
        if (unbalLoad != 0)
            s << "\tunbalanced Load: " << *unbalLoad;
        if (reaction != 0)
            s << "\treaction: " << *reaction;
        if (mass != 0) {
            s << "\tMass : " << *mass;
            s << "\tRayleigh Factor: alphaM: " << alphaM << endln;
            s << "\tRayleigh Forces: " << *this->getResponse(NodeData::RayleighForces);
        }
        if (theEigenvectors != 0)
            s << "\t Eigenvectors: " << *theEigenvectors;
        if (theDOF_GroupPtr != 0)
            s << "\tID : " << theDOF_GroupPtr->getID();
        s << "\n";
    }

    else if (flag == 1) { // print out: nodeId displacements
        s << this->getTag() << "  " << *commitDisp;
    }

    else if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << OPS_PRINT_JSON_NODE_INDENT << "{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"ndf\": " << numberDOF << ", ";
        s << "\"crd\": [";
        int numCrd = Crd->Size();
        for (int i = 0; i < numCrd - 1; i++)
            s << (*Crd)(i) << ", ";
        s << (*Crd)(numCrd - 1) << "]";
        if (mass != 0) {
            s << ", \"mass\": [";
            for (int i = 0; i < numberDOF - 1; i++)
                s << (*mass)(i, i) << ", ";
            s << (*mass)(numberDOF - 1, numberDOF - 1) << "]";
        }
        s << "}";
    }
}


// AddingSensitivity:BEGIN ///////////////////////////////////////

Matrix
Node::getMassSensitivity(void)
{
  if (index == -1)
      setGlobalMatrices();

  if (mass == 0) {
    theMatrices[index]->Zero();
    return *theMatrices[index];

  } else {
    Matrix massSens(mass->noRows(),mass->noCols());
    if ( (parameterID == 1) || (parameterID == 2) || (parameterID == 3) ) {
            massSens(parameterID-1,parameterID-1) = 1.0;
    }
    if (parameterID == 7) {
      massSens(0,0) = 1.0;
      massSens(1,1) = 1.0;
    }
    if (parameterID == 8) {
      massSens(0,0) = 1.0;
      massSens(1,1) = 1.0;
      massSens(2,2) = 1.0;
    }
    return massSens;
  }
}


int
Node::getCrdsSensitivity(void)
{
  if ( (parameterID == 4) || (parameterID == 5) || (parameterID == 6) ) {
    return (parameterID-3);
  }
  else {
    return 0;
  }
}


int
Node::setParameter(const char **argv, int argc, Parameter &param)
{
  // The following parameterID map is being used:
  // 1: nodal mass in direction 1
  // 2: nodal mass in direction 2
  // 3: nodal mass in direction 3
  // 4: coordinate in direction 1
  // 5: coordinate in direction 2
  // 6: coordinate in direction 3

  if (argc < 2)
    return -1;

  if ((strstr(argv[0],"mass") != 0) || (strstr(argv[0],"-mass") != 0)) {
    int direction = 0; // atoi(argv[1]);
    if ((strcmp(argv[1],"x") == 0)||(strcmp(argv[1],"X") == 0)||(strcmp(argv[1],"1") == 0)) {
      direction = 1;
      if (mass != 0)
      param.setValue((*mass)(0,0));
    }
    else if ((strcmp(argv[1],"y") == 0)||(strcmp(argv[1],"Y") == 0)||(strcmp(argv[1],"2") == 0)) {
      direction = 2;
      if (mass != 0)
      param.setValue((*mass)(1,1));
    }
    else if ((strcmp(argv[1],"z") == 0)||(strcmp(argv[1],"Z") == 0)||(strcmp(argv[1],"3") == 0)) {
      direction = 3;
      if (mass != 0)
      param.setValue((*mass)(2,2));
    }
    else if ((strcmp(argv[1],"xy") == 0)||(strcmp(argv[1],"XY") == 0)) {
      direction = 7;
      if (mass != 0)
      param.setValue((*mass)(0,0));
    }
    else if ((strcmp(argv[1],"xyz") == 0)||(strcmp(argv[1],"XYZ") == 0)) {
      direction = 8;
      if (mass != 0)
      param.setValue((*mass)(0,0));
    }

    if ((direction >= 1 && direction <= 3) || direction == 7 || direction == 8)
      return param.addObject(direction, this);
  }
  else if (strstr(argv[0],"coord") != 0) {
    int direction = atoi(argv[1]);
    if (direction >= 1 && direction <= 3) {
      if (Crd != 0)
      param.setValue((*Crd)(direction-1));
      return param.addObject(direction+3, this);
    }
  }
  else
    opserr << "WARNING: Could not set parameter in Node. " << endln;

  return -1;
}



int
Node::updateParameter(int pparameterID, Information &info)
{
  if (pparameterID >= 1 && pparameterID <= 3)
    (*mass)(pparameterID-1,pparameterID-1) = info.theDouble;

  else if (pparameterID == 7) {
    (*mass)(0,0) = info.theDouble;
    (*mass)(1,1) = info.theDouble;
  } else if (pparameterID == 8) {
    (*mass)(0,0) = info.theDouble;
    (*mass)(1,1) = info.theDouble;
    (*mass)(2,2) = info.theDouble;
  }

  else if (pparameterID >= 4 && pparameterID <= 6) {

    if ( (*Crd)(pparameterID-4) != info.theDouble) {

      // Set the new coordinate value
      (*Crd)(pparameterID-4) = info.theDouble;

      // Need to "setDomain" to make the change take effect.
      Domain *theDomain = this->getDomain();
      ElementIter &theElements = theDomain->getElements();
      Element *theElement;
      while ((theElement = theElements()) != 0) {
      theElement->setDomain(theDomain);
      }
    }
    else {
      // No change in nodal coordinate
    }
  }

  return -1;
}


int
Node::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;

  return 0;
}

int
Node::saveDispSensitivity(const Vector &v, int gradIndex, int numGrads)
{
  // If the sensitivity matrices are not already created:
  if (dispSensitivity == 0) {
    dispSensitivity = new Matrix( numberDOF, numGrads );
  }

  if (dispSensitivity->noRows() != numberDOF ||
      dispSensitivity->noCols() != numGrads) {
    delete dispSensitivity;
    dispSensitivity = new Matrix( numberDOF, numGrads );
  }

  //opserr << "Node::saveDispSens " << dispSensitivity->noRows() << ' ' << dispSensitivity->noCols() << endln;
  for (int i=0; i<numberDOF; i++ )
    (*dispSensitivity)(i,gradIndex) = v(i);

  return 0;
}

int
Node::saveVelSensitivity(const Vector &vdot, int gradIndex, int numGrads)
{
  // If the sensitivity matrices are not already created:
  if (velSensitivity == 0) {
    velSensitivity = new Matrix( numberDOF, numGrads );
  }

  for (int i=0; i<numberDOF; i++ )
    (*velSensitivity)(i,gradIndex) = vdot(i);

  return 0;
}

int
Node::saveAccelSensitivity(const Vector &vdotdot, int gradIndex, int numGrads)
{
  // If the sensitivity matrices are not already created:
  if (accSensitivity == 0) {
    accSensitivity = new Matrix( numberDOF, numGrads );
  }

  for (int i=0; i<numberDOF; i++ )
    (*accSensitivity)(i,gradIndex) = vdotdot(i);

  return 0;
}

double
Node::getDispSensitivity(int dof, int gradIndex)
{
  if (dispSensitivity != 0)
    return (*dispSensitivity)(dof-1,gradIndex);
  else
    return 0.0;
}

double
Node::getVelSensitivity(int dof, int gradIndex)
{
  if (velSensitivity != 0)
    return (*velSensitivity)(dof-1,gradIndex);
  else
    return 0.0;
}

double
Node::getAccSensitivity(int dof, int gradIndex)
{
  if (accSensitivity != 0)
    return (*accSensitivity)(dof-1,gradIndex);
  else
    return 0.0;
}
// AddingSensitivity:END /////////////////////////////////////////



const Vector &
Node::getReaction() {
  if (reaction == 0)
    reaction = new Vector(numberDOF);

  // TODO; zero this?

  return *reaction;
}

int
Node::addReactionForce(const Vector &add, double factor){

  // create rection vector if have not done so already
  if (reaction == 0)
    reaction = new Vector(numberDOF);

  // check vector of appropraie size
  if (add.Size() != numberDOF) {
    opserr << "WARNING Node::addReactionForce() - vector not of correct size\n";
    return -1;
  }

  if (factor == 1.0)
    *reaction += add;
  else if (factor == -1.0)
    *reaction -= add;
  else
    *reaction = add * factor;

  return 0;
}

int
Node::resetReactionForce(int flag)
{
  // TODO: Use NodeData for flag
  // create rection vector if have not done so already
  if (reaction == nullptr)
    reaction = new Vector(numberDOF);

  reaction->Zero();

  // add unbalance, the negative of applied forces hence the -=
  if (flag == 0) {
    *reaction -= this->getUnbalancedLoad();
  } if (flag == 1) {
    *reaction -= this->getUnbalancedLoadIncInertia();
  } else {
    if (mass != 0 && alphaM != 0) {
      if (alphaM != 0.0) {
      const Vector &theVel = this->getTrialVel(); // in case vel not created
      reaction->addMatrixVector(1.0, *mass, theVel, alphaM);
      }
    }
  }
  return 0;
}

int
Node::fillResponse(NodeData responseType, Vector& result, int offset)
{
  const Vector *resp = getResponse(responseType);

  if (resp == nullptr)
    return -1;

  if (resp->Size() + offset > result.Size())
    result.resize(offset + resp->Size());
  for (int i = 0; i < resp->Size(); i++)
    result(i+offset) = (*resp)(i);
  return resp->Size();
}

const Vector *
Node::getResponse(NodeData responseType)
{
  const Vector *result = NULL;
  if (responseType == NodeData::Disp)
    result  = &(this->getDisp());

  else if (responseType == NodeData::Vel)
    return &(this->getVel());

  else if (responseType == NodeData::Accel)
    return &(this->getAccel());

  else if (responseType == NodeData::IncrDisp)
    return &(this->getIncrDisp());

  else if (responseType == NodeData::IncrDeltaDisp)
    return &(this->getIncrDeltaDisp());

  else if (responseType == NodeData::Reaction)
    return &(this->getReaction());

  else if (responseType == NodeData::UnbalancedLoad)
    return &(this->getUnbalancedLoad());

  else if (responseType == NodeData::UnbalanceInclInertia)
    return &(this->getUnbalancedLoadIncInertia());

  else if (responseType == NodeData::RayleighForces) {
    if (unbalLoadWithInertia == 0) {
      unbalLoadWithInertia = new Vector(this->getUnbalancedLoad());
    }
    if (alphaM != 0.0 && mass != 0) {
      const Vector &theVel = this->getTrialVel(); // in case vel not created
      unbalLoadWithInertia->addMatrixVector(0.0, *mass, theVel, -alphaM);
    } else
      unbalLoadWithInertia->Zero();

    return unbalLoadWithInertia;

  } 
  else
    return NULL;

  return result;
}

void
Node::setCrds(double Crd1)
{
  if (Crd != 0 && Crd->Size() >= 1)
    (*Crd)(0) = Crd1;

  // Need to "setDomain" to make the change take effect.
  Domain *theDomain = this->getDomain();
  ElementIter &theElements = theDomain->getElements();
  Element *theElement;
  while ((theElement = theElements()) != nullptr)
    theElement->setDomain(theDomain);

}

void
Node::setCrds(double Crd1, double Crd2)
{
  if (Crd != 0 && Crd->Size() >= 2) {
    (*Crd)(0) = Crd1;
    (*Crd)(1) = Crd2;

    // Need to "setDomain" to make the change take effect.
    Domain *theDomain = this->getDomain();
    ElementIter &theElements = theDomain->getElements();
    Element *theElement;
    while ((theElement = theElements()) != nullptr)
      theElement->setDomain(theDomain);
  }
}

void
Node::setCrds(double Crd1, double Crd2, double Crd3)
{
  if (Crd != 0 && Crd->Size() >= 3) {
    (*Crd)(0) = Crd1;
    (*Crd)(1) = Crd2;
    (*Crd)(2) = Crd3;

    // Need to "setDomain" to make the change take effect.
    Domain *theDomain = this->getDomain();
    ElementIter &theElements = theDomain->getElements();
    Element *theElement;
    while ((theElement = theElements()) != nullptr)
      theElement->setDomain(theDomain);
  }
}

void
Node::setCrds(const Vector &newCrds)
{
  if (Crd != 0 && Crd->Size() == newCrds.Size()) {
    (*Crd) = newCrds;

      return;

    // Need to "setDomain" to make the change take effect.
    Domain *theDomain = this->getDomain();
    ElementIter &theElements = theDomain->getElements();
    Element *theElement;
    while ((theElement = theElements()) != nullptr)
      theElement->setDomain(theDomain);
  }
}

#if 1
int
Node::getDisplayRots(Vector& res, double fact, int mode)
{
  int ndm = Crd->Size();
  int resSize = res.Size();
  int nRotDOFs = numberDOF - ndm;
  if (resSize < nRotDOFs)
      return -1;

  if (mode < 0) {
    int eigenMode = -mode;
    for (int i = ndm; i < resSize; i++)
        res(i) = (*theEigenvectors)(i, eigenMode - 1) * fact;

  } else {
    for (int i = ndm; i < resSize; i++)
        res(i) = (*commitDisp)(i) * fact;
  }

  // zero rest
  for (int i = nRotDOFs; i < resSize; i++)
    res(i) = 0;

  return 0;
}

int
Node::getDisplayCrds(Vector &res, double fact, int mode)
{
  return -1;
}

int
Node::setDisplayCrds(const Vector &theCrds)
{
  return -1;
}
#endif


//Add Pointer to NodalThermalAction id applicable------begin-----L.Jiang, [SIF]
NodalThermalAction*
Node::getNodalThermalActionPtr(void)
{
    return theNodalThermalActionPtr;
}
void
Node::setNodalThermalActionPtr(NodalThermalAction* theAction)
{
    theNodalThermalActionPtr = theAction;
}
//Add Pointer to NodalThermalAction id applicable-----end------L.Jiang, {SIF]

int
Node::setGlobalMatrices()
{
    if (index == -1) {
      for (int i=0; i<numMatrices; i++) {
          if (theMatrices[i]->noRows() == numberDOF) {
            index = i;
            i = numMatrices;
          }
      }
    }
    if (index == -1) {
      Matrix **nextMatrices = new Matrix *[numMatrices+1];

      for (int j=0; j<numMatrices; j++)
          nextMatrices[j] = theMatrices[j];

      Matrix *theMatrix = new Matrix(numberDOF, numberDOF);

      nextMatrices[numMatrices] = theMatrix;
      if (numMatrices != 0)
          delete [] theMatrices;
      index = numMatrices;
      numMatrices++;
      theMatrices = nextMatrices;
    }

    return 0;
}
