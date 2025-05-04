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
// Purpose: This file contains the implementation of class SP_Constraint.
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
#include <string.h>
#include <SP_Constraint.h>
#include <classTags.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Domain.h>
#include <Node.h>
#include <ID.h>
#include <TaggedObject.h> // for PRINT flags

static int numSPs = 0;
static int nextTag = 0;

// 2 little procedures needed for parallel processing all due to fact that SP's need 
// to keep unique tags among processes in parallel

int SP_Constraint_SetNextTag(int next) {
  nextTag = next;
  return nextTag;
}

int SP_Constraint_GetNextTag() {
  return nextTag;
}

// constructor for FEM_ObjectBroker
SP_Constraint::SP_Constraint(int clasTag)
:DomainComponent(0,clasTag),
 nodeTag(0), dofNumber(0), valueR(0.0), valueC(0.0), initialValue(0.0), 
 initialized(false), isConstant(true), retZeroInitValue(true),
 loadPatternTag(-1)
{
  numSPs++;
}

// constructor for a subclass to use
SP_Constraint::SP_Constraint(int node, int ndof, int clasTag)
:DomainComponent(nextTag++, clasTag),
 nodeTag(node), dofNumber(ndof), valueR(0.0), valueC(0.0), initialValue(0.0), 
 initialized(false), isConstant(true), retZeroInitValue(true),
 loadPatternTag(-1)
 // valueC is set to 1.0 so that homo will be false when recvSelf() invoked
 // should be ok as valueC cannot be used by subclasses and subclasses should
 // not be used if it is a homogeneous constraint.
{
  numSPs++;
}

// constructor for object of type SP_Constraint
SP_Constraint::SP_Constraint(int node, int ndof, double value, bool ISconstant, bool retZeroInitValue)
:DomainComponent(nextTag++, CNSTRNT_TAG_SP_Constraint),
 nodeTag(node), dofNumber(ndof), valueR(value), valueC(value), initialValue(0.0), 
 initialized(false), isConstant(ISconstant), retZeroInitValue(retZeroInitValue),
 loadPatternTag(-1)
{
  numSPs++;
}

SP_Constraint::~SP_Constraint()
{
  numSPs--;
  if (numSPs == 0)
    nextTag = 0;
}

int
SP_Constraint::getNodeTag() const
{
    // return id of constrained node
    return nodeTag;
}

int
SP_Constraint::getDOF_Number() const
{
    //  return the number of the constrained DOF    
    return dofNumber;
}


double
SP_Constraint::getValue()
{
    // return the value of the constraint
    return valueC;
}

double
SP_Constraint::getInitialValue()
{
  if (retZeroInitValue == false)    
    // return the initial value of the constraint if retZeroInitValue is false,
    //   - constraint handlers will subtract this off the current getValue()    
    return initialValue;
  
  else
    return 0;
}

int
SP_Constraint::applyConstraint(double loadFactor)
{
    // as SP_Constraint objects are time invariant nothing is done
    if (isConstant == false)
      valueC = loadFactor*valueR;
    return 0;
}


bool
SP_Constraint::isHomogeneous() const
{
    if (valueR == 0.0)
	return true;
    else
	return false;
}

void
SP_Constraint::setLoadPatternTag(int tag)
{
  loadPatternTag = tag;
}

int
SP_Constraint::getLoadPatternTag() const
{
  return loadPatternTag;
}

void
SP_Constraint::setDomain(Domain* theDomain)
{
    // store initial state
    if (theDomain) {
        if (!initialized) { // don't do it if setDomain called after recvSelf when already initialized!
            Node* theNode = theDomain->getNode(nodeTag);
            if (theNode == 0) {
                opserr << "FATAL SP_Constraint::setDomain() - Constrained";
                opserr << " Node does not exist in Domain\n";
                opserr << nodeTag << endln;
                return; // -1;
            }
            const Vector& U = theNode->getTrialDisp();
            if (dofNumber < 0 || dofNumber >= U.Size()) {
                opserr << "SP_Constraint::setDomain FATAL Error: Constrained DOF " << dofNumber << " out of bounds [0-" << U.Size() << "]\n";
                return; // -1;
            }
            initialValue = U(dofNumber);
            initialized = true;
        }
    }

    // call base class implementation
    return DomainComponent::setDomain(theDomain);
}

int 
SP_Constraint::sendSelf(int cTag, Channel &theChannel)
{
    static Vector data(10);  // we send as double to avoid having 
                     // to send two messages.
    data(0) = this->getTag(); 
    data(1) = nodeTag;
    data(2) = dofNumber;
    data(3) = valueC;
    if (isConstant == true)
	data(4) = 1.0;
    else
	data(4) = 0.0;
    data(5) = valueR;
    data(6) = this->getLoadPatternTag();

    data(7) = nextTag;
    data(8) = initialValue;
    data(9) = static_cast<double>(initialized);

    int result = theChannel.sendVector(this->getDbTag(), cTag, data);
    if (result != 0) {
      // opserr << "WARNING SP_Constraint::sendSelf - error sending Vector data\n";
      return result;
    }

    return 0;
}

int 
SP_Constraint::recvSelf(int cTag, Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
    static Vector data(10);  // we sent the data as double to avoid having to send
                     // two messages
    int result = theChannel.recvVector(this->getDbTag(), cTag, data);
    if (result < 0) {
	// opserr << "WARNING SP_Constraint::recvSelf - error receiving Vector data\n";
	return result;
    }
    
    // if o.k. set the data
    this->setTag((int)data(0));
    nodeTag = (int)data(1);
    dofNumber = (int)data(2);
    valueC = data(3);

    if (data(4) == 1.0)
	isConstant = true;
    else
	isConstant = false;
    valueR = data(5);
    valueC = valueR;
    this->setLoadPatternTag((int)data(6));

    nextTag = (int)data(7);
    initialValue = data(8);
    initialized = static_cast<bool>(data(9));

    return 0;
}


void
SP_Constraint::Print(OPS_Stream &s, int flag) 
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_ELEM_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"node\": " << nodeTag << ", ";
    s << "\"dof\": " << dofNumber+1 << ", ";
    s << "\"ref_value\": " << valueR << "}";

  } else {
    s << "SP_Constraint: " << this->getTag();
    s << "\t Node: " << nodeTag << " DOF: " << dofNumber+1;
    s << " ref value: " << valueR << " current value: " << valueC << " initial value: " << initialValue << endln;
  }
}








