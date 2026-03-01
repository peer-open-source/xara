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
// of the DOF_Group class interface.
//
// File: ~/analysis/dof_grp/DOF_Group.C
//
// Written: fmk 
// Created: 11/96
//
#include <DOF_Group.h>
#include <stdlib.h>
#include <assert.h>
#include <Node.h>
#include <Vector.h>
#include <Matrix.h>
#include <TransientIntegrator.h>

#define MAX_NUM_DOF 256


// static variables initialisation
Matrix **DOF_Group::theMatrices; // array of pointers to class wide matrices
Vector **DOF_Group::theVectors;  // array of pointers to class widde vectors
int DOF_Group::numDOFs(0);       // number of objects


DOF_Group::DOF_Group(int tag, Node *node)
: TaggedObject(tag),
  unbalance(0), tangent(0), myNode(node), 
  myID(node->getNumberDOF()), 
  numDOF(node->getNumberDOF())
{
  // get number of DOF & verify valid
  int numDOF = node->getNumberDOF();
  if (numDOF <= 0) {
    opserr << "DOF_Group::DOF_Group(Node *) ";
    opserr << " node must have at least 1 dof " << *node;
    exit(-1);
  }

  // initially set all the IDs to be -2
  for (int i=0; i<numDOF; i++)
    myID(i) = -2;

  // if this is the first DOF_Group we now
  // create the arrays used to store pointers to class wide
  // matrix and vector objects used to return tangent and residual
  if (numDOFs == 0) {
    theMatrices = new Matrix *[MAX_NUM_DOF+1];
    theVectors  = new Vector *[MAX_NUM_DOF+1];
    
    for (int i=0; i<MAX_NUM_DOF; i++) {
      theMatrices[i] = 0;
      theVectors[i] = 0;
    }
  }    
  
  // set the pointers for the tangent and residual
  if (numDOF <= MAX_NUM_DOF) {
    // use class wide objects
    if (theVectors[numDOF] == 0) {
      // have to create matrix and vector of size as none yet created
      theVectors[numDOF] = new Vector(numDOF);
      theMatrices[numDOF] = new Matrix(numDOF,numDOF);
      unbalance = theVectors[numDOF];
      tangent = theMatrices[numDOF];

    } else {
      unbalance = theVectors[numDOF];
      tangent = theMatrices[numDOF];
    }
  } else {
    // create matrices and vectors for each object instance
    unbalance = new Vector(numDOF);
    tangent = new Matrix(numDOF, numDOF);
  }

  numDOFs++;
}


DOF_Group::DOF_Group(int tag, int ndof)
:TaggedObject(tag),
 unbalance(0), tangent(0), myNode(0), 
 myID(ndof), 
 numDOF(ndof)
{
  // get number of DOF & verify valid
  int numDOF = ndof;
  assert(numDOF > 0);

  // initially set all the IDs to be -2
  for (int i=0; i<numDOF; i++)
    myID(i) = -2;
  
  // if this is the first DOF_Group we now
  // create the arrays used to store pointers to class wide
  // matrix and vector objects used to return tangent and residual
  if (numDOFs == 0) {
    theMatrices = new Matrix *[MAX_NUM_DOF+1];
    theVectors  = new Vector *[MAX_NUM_DOF+1];
    
    for (int i=0; i<MAX_NUM_DOF; i++) {
        theMatrices[i] = 0;
        theVectors[i]  = 0;
    }
  }

  // set the pointers for the tangent and residual
  if (numDOF <= MAX_NUM_DOF) {
    // use class wide objects
    if (theVectors[numDOF] == 0) {
        // have to create matrix and vector of size as none yet created
        theVectors[numDOF] = new Vector(numDOF);
        theMatrices[numDOF] = new Matrix(numDOF,numDOF);
        unbalance = theVectors[numDOF];
        tangent   = theMatrices[numDOF];
    } else {
        unbalance = theVectors[numDOF];
        tangent   = theMatrices[numDOF];
    }
  } else {
    // create matrices and vectors for each object instance
    unbalance = new Vector(numDOF);
    tangent   = new Matrix(numDOF, numDOF);
  }
  numDOFs++;
}


DOF_Group::~DOF_Group()
{
  numDOFs--;

  int numDOF = unbalance->Size();

  // set the pointer in the associated Node to 0, to stop
  // segmentation fault if node tries to use this object after destroyed
  if (myNode != 0) 
    myNode->setDOF_GroupPtr(0);

  // delete tangent and residual if created specially
  if (numDOF > MAX_NUM_DOF) {
      if (tangent != 0) delete tangent;
      if (unbalance != 0) delete unbalance;
  }

  // if this is the last FE_Element, clean up the
  // storage for the matrix and vector objects
  if (numDOFs == 0) {
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

// void setID(int index, int value);
//	Method to set the corresponding index of the ID to value.

void
DOF_Group::setID(int index, int value)
{
  assert((index >=0) && (index < numDOF));
  myID(index) = value;
}

// void setID(const ID &);
//	Method to set the ID to be same as that given.

void
DOF_Group::setID(const ID &copy)
{
  myID = copy;
}
 

const ID &
DOF_Group::getID() const
{
  return myID;
}

int
DOF_Group::doneID()
{
  return 0;
}

int
DOF_Group::getNumDOF() const
{
  return numDOF;
}

int
DOF_Group::getNodeTag() const
{
  if (myNode != nullptr)
    return myNode->getTag();
  else
    return -1;
}

int
DOF_Group::getNumFreeDOF() const
{
  int numFreeDOF = numDOF;
  for (int i=0; i<numDOF; i++)
    if (myID(i) == -1 || myID(i) == -4)
      numFreeDOF--;
  
  return numFreeDOF;
}

int
DOF_Group::getNumConstrainedDOF() const
{   
  int numConstr = 0;
  for (int i=0; i<numDOF; i++)
    if (myID(i) < 0)
      numConstr++;    

  return numConstr;
}    



const Matrix &
DOF_Group::getTangent(Integrator *theIntegrator) 
{	
  if (theIntegrator != nullptr)
      theIntegrator->formNodTangent(this);    
  return *tangent;
}

void  
DOF_Group::zeroTangent()
{
  tangent->Zero();
}


void  
DOF_Group::addMtoTang(double fact)
{
  // if there is no associated node, subclass should 
  // implement this method (ie addMtoTang())
  assert(myNode != nullptr);
  tangent->addMatrix(1.0, myNode->getMass(), fact);
}


void  
DOF_Group::addCtoTang(double fact)
{
  // if there is no associated node, subclass should 
  // implement this method (ie addCtoTang())
  assert(myNode != nullptr);
  tangent->addMatrix(1.0, myNode->getDamp(), fact);
}



void
DOF_Group::zeroUnbalance() 
{
  unbalance->Zero();
}


const Vector &
DOF_Group::getUnbalance(Integrator *theIntegrator)
{
  if (theIntegrator != nullptr)
    theIntegrator->formNodUnbalance(this);

  return *unbalance;
}


void
DOF_Group::addPtoUnbalance(double fact)
{
  // if there is no associated node, subclass should 
  // implement this method (ie addPtoUnbalance())
  assert(myNode != nullptr);
  unbalance->addVector(1.0, myNode->getUnbalancedLoad(), fact);
}


void
DOF_Group::addPIncInertiaToUnbalance(double fact)
{
  // if there is no associated node, subclass should 
  // implement this method (ie addPIncInertiaToUnbalance())
  assert(myNode != nullptr);
  unbalance->addVector(1.0, myNode->getUnbalancedLoadIncInertia(), fact);
}


void  
DOF_Group::addM_Force(const Vector &Udotdot, double fact)
{
  assert(myNode != nullptr);
//  if (myNode == 0) {
//      opserr << "DOF_Group::addM_Force() - no Node associated";	
//      opserr << " subclass should not call this method \n";	    
//      return;
//  }

    Vector accel(numDOF);
    // get accel for the unconstrained dof
    for (int i=0; i<numDOF; i++) {
      int loc = myID(i);
      if (loc >= 0)
          accel(i) = Udotdot(loc); 
      else accel(i) = 0.0;
    }

    unbalance->addMatrixVector(1.0, myNode->getMass(), accel, fact);
}



const Vector &
DOF_Group::getTangForce(const Vector &Udotdot, double fact)
{
  opserr << "DOF_Group::getTangForce() - not yet implemented";
  return *unbalance;
}


const Vector &
DOF_Group::getM_Force(const Vector &Udotdot, double fact)
{
  if (myNode == 0) {
    opserr << "DOF_Group::getM_Force() - no Node associated";	
    opserr << " subclass should not call this method \n";	    
    return *unbalance;
  }

  Vector accel(numDOF);
  // get accel for the unconstrained dof
  for (int i=0; i<numDOF; i++) {
    int loc = myID(i);
    if (loc >= 0)
        accel(i) = Udotdot(loc); 
    else accel(i) = 0.0;
  }
  
  unbalance->Zero();
  unbalance->addMatrixVector(myNode->getMass(), accel, fact);
  
  return *unbalance;
}



const Vector &
DOF_Group::getC_Force(const Vector &Udotdot, double fact)
{
  assert(myNode != nullptr);

  Vector accel(numDOF);
  // get accel for the unconstrained dof
  for (int i=0; i<numDOF; i++) {
      int loc = myID(i);
      if (loc >= 0)
          accel(i) = Udotdot(loc); 
      else accel(i) = 0.0;
  }
      
  unbalance->addMatrixVector(0.0, myNode->getDamp(), accel, fact);
  return *unbalance;
}


const Vector & 
DOF_Group::getCommittedDisp()
{
  assert(myNode != nullptr);
  return myNode->getDisp();
}


const Vector & 
DOF_Group::getCommittedVel()
{
  assert(myNode != nullptr);
  return myNode->getVel();
}


const Vector & 
DOF_Group::getCommittedAccel()
{
  assert(myNode != nullptr);
  return myNode->getAccel();
}

// void setNodeDisp(const Vector &u);
//	Method to set the corresponding nodes displacements to the
//	values in u, components identified by myID;
void
DOF_Group::setNodeDisp(const Vector &u)
{
  assert(myNode != nullptr); 
  Vector &disp = *unbalance;
  disp = myNode->getTrialDisp();
  int i;
  
  // get disp for my dof out of vector u
  for (i=0; i<numDOF; i++) {
    int loc = myID(i);
    if (loc >= 0)
      disp(i) = u(loc);  
  }

  myNode->setTrialDisp(disp);
}
	
	
// void setNodeVel(const Vector &udot);
//	Method to set the corresponding nodes velocities to the
//	values in udot, components identified by myID;

void
DOF_Group::setNodeVel(const Vector &udot)
{
  assert(myNode != nullptr);

  Vector &vel = *unbalance;
  vel = myNode->getTrialVel();
  int i;
  
  // get vel for my dof out of vector udot
  for (i=0; i<numDOF; i++) {
    int loc = myID(i);	    	
    if (loc >= 0) 
      vel(i) = udot(loc);  
  }

  myNode->setTrialVel(vel);
}


// void setNodeAccel(const Vector &udotdot);
//	Method to set the corresponding nodes accelerations to the
//	values in udotdot, components identified by myID;
void
DOF_Group::setNodeAccel(const Vector &udotdot)
{

  assert(myNode != nullptr);

  Vector &accel = *unbalance;;
  accel = myNode->getTrialAccel();
  int i;
  
  // get disp for the unconstrained dof
  for (i=0; i<numDOF; i++) {
      int loc = myID(i);	    		
      if (loc >= 0)
          accel(i) = udotdot(loc); 
  }

  myNode->setTrialAccel(accel);
}


// void setNodeIncrDisp(const Vector &u);
//	Method to set the corresponding nodes displacements to the
//	values in u, components identified by myID;

void
DOF_Group::incrNodeDisp(const Vector &u)
{
  assert(myNode != nullptr);

  Vector &disp = *unbalance;

  assert(disp.Size() != 0);

  // get disp for my dof out of vector u
  for (int i=0; i<numDOF; i++) {
      int loc = myID(i);	    			
      if (loc >= 0)
        disp(i) = u(loc);  
      else 
        disp(i) = 0.0;  
  }

  myNode->incrTrialDisp(disp);
}
	
	
// void setNodeincrVel(const Vector &udot);
//	Method to set the corresponding nodes velocities to the
//	values in udot, components identified by myID;

void
DOF_Group::incrNodeVel(const Vector &udot)
{
  assert(myNode != nullptr);
    
  Vector &vel = *unbalance;
  
  // get vel for my dof out of vector udot
  for (int i=0; i<numDOF; i++) {
      int loc = myID(i);
      if (loc >= 0)
        vel(i) = udot(loc);  // -1 for dof labelled 1 through ndof
      else
        vel(i) = 0.0;
  }
  myNode->incrTrialVel(vel);
}



// void setNodeIncrAccel(const Vector &udotdot);
//	Method to set the corresponding nodes accelerations to the
//	values in udotdot, components identified by myID;

void
DOF_Group::incrNodeAccel(const Vector &udotdot)
{

  assert(myNode != nullptr);

  Vector &accel = *unbalance;
  
  // get disp for the unconstrained dof
  for (int i=0; i<numDOF; i++) {
      int loc = myID(i);
      if (loc >= 0)
        accel(i) = udotdot(loc); 
      else 
        accel(i) = 0.0;
  }    
  myNode->incrTrialAccel(accel);
}

const Vector & 
DOF_Group::getTrialDisp()
{
  assert(myNode != nullptr);
  return myNode->getTrialDisp();
}


const Vector & 
DOF_Group::getTrialVel()
{
  assert(myNode != nullptr);
  return myNode->getTrialVel();
}


const Vector & 
DOF_Group::getTrialAccel()
{
  assert(myNode != nullptr);
  return myNode->getTrialAccel();
}


void
DOF_Group::setEigenvector(int mode, const Vector &theVector)
{
  assert(myNode != nullptr);
  Vector &eigenvector = *unbalance;
  
  // get disp for the unconstrained dof
  for (int i=0; i<numDOF; i++) {
      int loc = myID(i);
      if (loc >= 0)
          eigenvector(i) = theVector(loc); 
      else eigenvector(i) = 0.0;
  }    
  myNode->setEigenvector(mode, eigenvector);
}



const Matrix &
DOF_Group::getEigenvectors()
{
  assert(myNode != nullptr);
  return myNode->getEigenvectors();
}



Matrix *
DOF_Group::getT()
{
  return 0;
}



void  
DOF_Group::addLocalM_Force(const Vector &accel, double fact)
{
  assert(myNode != nullptr);
  unbalance->addMatrixVector(1.0, myNode->getMass(), accel, fact);
}


// AddingSensitivity:BEGIN ////////////////////////////////////////
const Vector &
DOF_Group::getDispSensitivity(int gradNumber)
{
  Vector &result = *unbalance;
  for (int i=0; i<numDOF; i++) {
    result(i) = myNode->getDispSensitivity(i+1,gradNumber);
  }
  return result;
}
	
const Vector &
DOF_Group::getVelSensitivity(int gradNumber)
{
    Vector &result = *unbalance;
    for (int i=0; i<numDOF; i++)
      result(i) = myNode->getVelSensitivity(i+1,gradNumber);

    return result;
}


const Vector &
DOF_Group::getAccSensitivity(int gradNumber)
{
  Vector &result = *unbalance;
  for (int i=0; i<numDOF; i++)
    result(i) = myNode->getAccSensitivity(i+1,gradNumber);

  return result;
}
	
	

int 
DOF_Group::saveDispSensitivity(const Vector &v, int gradNum, int numGrads)
{
  Vector &dudh = *unbalance;

  for (int i = 0; i < numDOF; i++) {
    int loc = myID(i);
    if (loc >= 0)
      dudh(i) = v(loc);
    else
      dudh(i) = 0.0;
  }

  return myNode->saveDispSensitivity(dudh, gradNum, numGrads);
}


int 
DOF_Group::saveVelSensitivity(const Vector &v, int gradNum, int numGrads)
{
  Vector &dudh = *unbalance;

  for (int i = 0; i < numDOF; i++) {
    int loc = myID(i);
    if (loc >= 0)
      dudh(i) = v(loc);
    else
      dudh(i) = 0.0;
  }

  return myNode->saveVelSensitivity(dudh, gradNum, numGrads);
}

int 
DOF_Group::saveAccSensitivity(const Vector &v, int gradNum, int numGrads)
{
  Vector &dudh = *unbalance;

  for (int i = 0; i < numDOF; i++) {
    int loc = myID(i);
    if (loc >= 0)
      dudh(i) = v(loc);
    else
      dudh(i) = 0.0;
  }

  return myNode->saveAccelSensitivity(dudh, gradNum, numGrads);
}


int 
DOF_Group::saveSensitivity(const Vector &v, const Vector &vdot,
			   const Vector &vdotdot, int gradNum, int numGrads)
{
  int ok = 0;

  ok += this->saveDispSensitivity(v, gradNum, numGrads);
  ok += this->saveVelSensitivity(vdot, gradNum, numGrads);
  ok += this->saveAccSensitivity(vdotdot, gradNum, numGrads);

  return ok;
}

void  
DOF_Group::addM_ForceSensitivity(const Vector &Udotdot, double fact)
{
  assert(myNode != nullptr);

  Vector accel(numDOF);
  // get accel for the unconstrained dof
  for (int i=0; i<numDOF; i++) {
      int loc = myID(i);
      if (loc >= 0)
          accel(i) = Udotdot(loc); 
      else accel(i) = 0.0;
  }
      
  unbalance->addMatrixVector(1.0, myNode->getMassSensitivity(), accel, fact);
}

void
DOF_Group::addD_Force(const Vector &Udot, double fact)
{
  assert(myNode != nullptr);

  Vector vel(numDOF);
  // get velocity for the unconstrained dof
  for (int i=0; i<numDOF; i++) {
      int loc = myID(i);
      if (loc >= 0)
          vel(i) = Udot(loc);
      else vel(i) = 0.0;
  }

  unbalance->addMatrixVector(1.0, myNode->getDamp(), vel, fact);
}

void
DOF_Group::addD_ForceSensitivity(const Vector &Udot, double fact)
{
  assert(myNode != nullptr);

  Vector vel(numDOF);
  // get velocity for the unconstrained dof
  for (int i=0; i<numDOF; i++) {
      int loc = myID(i);
      if (loc >= 0)
          vel(i) = Udot(loc);
      else vel(i) = 0.0;
  }

  unbalance->addMatrixVector(1.0, myNode->getDampSensitivity(), vel, fact);
}

// AddingSensitivity:END //////////////////////////////////////////
void
DOF_Group::resetNodePtr()
{
  myNode = 0;
}


double 
DOF_Group::getDampingBetaFactor(int mode, double ratio, double wn)
{
  // to return 2.0 * ratio * wn * phi(mode)' * M * v
  double beta = 0;
  const Matrix & mass = myNode->getMass();
  const Matrix & eigenVectors = myNode->getEigenvectors();
  const Vector & vel = myNode->getTrialVel();
  int numDOF = eigenVectors.noRows();
  int numMode = eigenVectors.noCols();
  Vector Mv = mass * vel;
  if (mode < numMode) {
    for (int i = 0; i < numDOF; i++) 
      beta += 2.0 * ratio * wn * eigenVectors(i, mode) * Mv(i);
  }
  return beta;
}

const Vector &
DOF_Group::getDampingBetaForce(int mode, double beta)
{
  // to return beta * M * phi(mode)
  const Matrix & mass = myNode->getMass();
  const Matrix & eigenVectors = myNode->getEigenvectors();
  int numDOF = eigenVectors.noRows();

  Vector eigenvector(numDOF);

  for (int i=0; i<numDOF; i++)
    eigenvector(i) = eigenVectors(i,mode);

  unbalance->addMatrixVector(0.0, mass, eigenvector, -beta);
  return *unbalance;
}

