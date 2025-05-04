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
// Purpose: This file contains the implementation of class MP_Constraint.
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
#include <MP_Constraint.h>
#include <string.h>
#include <Matrix.h>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Domain.h>
#include <Node.h>

static int numMPs = 0;
static int nextTag = 0;

 
// constructor for FEM_ObjectBroker			// Arash
MP_Constraint::MP_Constraint(int clasTag )		
:DomainComponent(nextTag++, clasTag),
 nodeRetained(0),nodeConstrained(0),constraint(0),constrDOF(0),retainDOF(0), initialized(false),
 dbTag1(0), dbTag2(0)
{
  numMPs++;
}

// constructor for Subclass
MP_Constraint::MP_Constraint(int nodeRetain, int nodeConstr, 
			     ID &constrainedDOF, 
			     ID &retainedDOF, int clasTag)
:DomainComponent(nextTag++, clasTag),
 nodeRetained(nodeRetain), nodeConstrained(nodeConstr), 
 constraint(0), constrDOF(0), retainDOF(0), initialized(false), dbTag1(0), dbTag2(0)
{
  numMPs++;
  
  constrDOF = new ID(constrainedDOF);
  retainDOF = new ID(retainedDOF);    
  if (constrDOF == 0 || constrainedDOF.Size() != constrDOF->Size() ||
      retainDOF == 0 || retainedDOF.Size() != retainDOF->Size()) { 
    opserr << "MP_Constraint::MP_Constraint - ran out of memory 1\n";
    exit(-1);
  }    

  // resize initial state
  Uc0.resize(constrDOF->Size());
  Uc0.Zero();
  Ur0.resize(retainDOF->Size());
  Ur0.Zero();
}


// general constructor for ModelBuilder
MP_Constraint::MP_Constraint(int nodeRetain, int nodeConstr, Matrix &constr,
			     ID &constrainedDOF, ID &retainedDOF)
:DomainComponent(nextTag++, CNSTRNT_TAG_MP_Constraint), 
 nodeRetained(nodeRetain), nodeConstrained(nodeConstr), 
 constraint(0), constrDOF(0), retainDOF(0), initialized(false), dbTag1(0), dbTag2(0)
{
  numMPs++;    
  constrDOF = new ID(constrainedDOF);
  retainDOF = new ID(retainedDOF);    
  
  constraint = new Matrix(constr);     

  // resize initial state
  Uc0.resize(constrDOF->Size());
  Uc0.Zero();
  Ur0.resize(retainDOF->Size());
  Ur0.Zero();
}



MP_Constraint::~MP_Constraint()
{
    // invoke the destructor on the matrix and the two ID objects
    if (constraint != 0)
	delete constraint;
    if (constrDOF != 0)
	delete constrDOF;
    if (retainDOF != 0)
	delete retainDOF;    
    
    numMPs--;
    if (numMPs == 0)
      nextTag = 0;
}

void MP_Constraint::setDomain(Domain* theDomain)
{
    // store initial state
    if (theDomain) {
        if (!initialized) { // don't do it if setDomain called after recvSelf when already initialized!
            Node* theRetainedNode = theDomain->getNode(nodeRetained);
            Node* theConstrainedNode = theDomain->getNode(nodeConstrained);
            if (theRetainedNode == 0 || theConstrainedNode == 0) {
                opserr << "FATAL MP_Constraint::setDomain() - Constrained or Retained";
                opserr << " Node does not exist in Domain\n";
                opserr << nodeRetained << " " << nodeConstrained << endln;
                exit(-1);
            }
            const Vector& Uc = theConstrainedNode->getTrialDisp();
            const Vector& Ur = theRetainedNode->getTrialDisp();
            const ID& idc = getConstrainedDOFs();
            const ID& idr = getRetainedDOFs();
            for (int i = 0; i < idc.Size(); ++i) {
                int cdof = idc(i);
                if (cdof < 0 || cdof >= Uc.Size()) {
                    opserr << "MP_Constraint::setDomain FATAL Error: Constrained DOF " << cdof << " out of bounds [0-" << Uc.Size() << "]\n";
                    exit(-1);
                }
                Uc0(i) = Uc(cdof);
            }
            for (int i = 0; i < idr.Size(); ++i) {
                int rdof = idr(i);
                if (rdof < 0 || rdof >= Ur.Size()) {
                    opserr << "MP_Constraint::setDomain FATAL Error: Retained DOF " << rdof << " out of bounds [0-" << Ur.Size() << "]\n";
                    exit(-1);
                }
                Ur0(i) = Ur(rdof);
            }
            initialized = true;
        }
    }

    // call base class implementation
    DomainComponent::setDomain(theDomain);
}

int
MP_Constraint::getNodeRetained() const
{
    // return id of retained node
    return nodeRetained;
}

int
MP_Constraint::getNodeConstrained() const
{
    // return id of constrained node    
    return nodeConstrained;
}


const ID &
MP_Constraint::getConstrainedDOFs() const
{
  assert(constrDOF != nullptr);
  // return the ID corresponding to constrained DOF of Ccr
  return *constrDOF;    
}


// return the ID corresponding to retained DOF of Ccr
const ID &
MP_Constraint::getRetainedDOFs() const
{
    assert(retainDOF != nullptr);
    return *retainDOF;    
}

// does nothing MP_Constraint objects are time invariant
int 
MP_Constraint::applyConstraint(double timeStamp)
{
    return 0;
}

bool
MP_Constraint::isTimeVarying() const
{
    return false;
}


// return the constraint matrix Ccr
const Matrix &
MP_Constraint::getConstraint()
{
    assert(constraint != nullptr);
    return *constraint;    
}

const Vector& MP_Constraint::getConstrainedDOFsInitialDisplacement() const
{
    return Uc0;
}

const Vector& MP_Constraint::getRetainedDOFsInitialDisplacement() const
{
    return Ur0;
}

int 
MP_Constraint::sendSelf(int cTag, Channel &theChannel)
{
    static ID data(11);
    int dataTag = this->getDbTag();

    data(0) = this->getTag(); 
    data(1) = nodeRetained;
    data(2) = nodeConstrained;
    if (constraint == 0) data(3) = 0; else data(3) = constraint->noRows();
    if (constraint == 0) data(4) = 0; else data(4) = constraint->noCols();    
    if (constrDOF == 0) data(5) = 0; else data(5) = constrDOF->Size();    
    if (retainDOF == 0) data(6) = 0; else data(6) = retainDOF->Size();        
    
    // need two database tags for ID objects
    if (constrDOF != 0 && dbTag1 == 0) 
      dbTag1 = theChannel.getDbTag();
    if (retainDOF != 0 && dbTag2 == 0) 
      dbTag2 = theChannel.getDbTag();

    data(7) = dbTag1;
    data(8) = dbTag2;
    data(9) = nextTag;
    data(10) = static_cast<int>(initialized);

    int result = theChannel.sendID(dataTag, cTag, data);
    if (result < 0) {
	// opserr << "WARNING MP_Constraint::sendSelf - error sending ID data\n";
	return result;  
    }    
    
    if (constraint != 0 && constraint->noRows() != 0) {
	int result = theChannel.sendMatrix(dataTag, cTag, *constraint);
	if (result < 0) {
	    // opserr << "WARNING MP_Constraint::sendSelf ";
	    // opserr << "- error sending Matrix data\n"; 
	    return result;  
	}
    }

    if (constrDOF != 0 && constrDOF->Size() != 0) {
	int result = theChannel.sendID(dbTag1, cTag, *constrDOF);
	if (result < 0) {
	    // opserr << "WARNING MP_Constraint::sendSelf ";
	    // opserr << "- error sending constrained data\n"; 
	    return result;  
	}
    }

    if (retainDOF != 0 && retainDOF->Size() != 0) {
	int result = theChannel.sendID(dbTag2, cTag, *retainDOF);
	if (result < 0) {
	    // opserr << "WARNING MP_Constraint::sendSelf ";
	    // opserr << "- error sending retained data\n"; 
	    return result;  
	}
    }
    
    // send initial displacement vectors.
    // we need 2 database tags because they have the same size,
    // but we can reuse the tags used for ID objects, since they go into different files
    if (Uc0.Size() > 0) {
        int result = theChannel.sendVector(dbTag1, cTag, Uc0);
        if (result < 0) {
            opserr << "WARNING MP_Constraint::sendSelf ";
            opserr << "- error sending constrained initial displacement\n";
            return result;
        }
    }
    if (Ur0.Size() > 0) {
        int result = theChannel.sendVector(dbTag2, cTag, Ur0);
        if (result < 0) {
            opserr << "WARNING MP_Constraint::sendSelf ";
            opserr << "- error sending retained initial displacement\n";
            return result;
        }
    }

    return 0;
}


int 
MP_Constraint::recvSelf(int cTag, Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
    int dataTag = this->getDbTag();
    static ID data(11);
    int result = theChannel.recvID(dataTag, cTag, data);
    if (result < 0) {
	// opserr << "WARNING MP_Constraint::recvSelf - error receiving ID data\n";
	return result;  
    }    

    this->setTag(data(0));
    nodeRetained = data(1);
    nodeConstrained = data(2);
    int numRows = data(3); 
    int numCols = data(4);
    dbTag1 = data(7);
    dbTag2 = data(8);
    nextTag = data(9);
    initialized = static_cast<bool>(data(10));

    if (numRows != 0 && numCols != 0) {
	constraint = new Matrix(numRows,numCols);
	
	int result = theChannel.recvMatrix(dataTag, cTag, *constraint);
	if (result < 0) {
	    // opserr << "WARNING MP_Constraint::recvSelf ";
	    // opserr << "- error receiving Matrix data\n"; 
	    return result;  
	}
    }    
    int size = data(5);
    if (size != 0) {
	constrDOF = new ID(size);
	int result = theChannel.recvID(dbTag1, cTag, *constrDOF);
	if (result < 0) {
	    // opserr << "WARNING MP_Constraint::recvSelf ";
	    // opserr << "- error receiving constrained data\n"; 
	    return result;  
	}	
    }
    
    size = data(6);
    if (size != 0) {
	retainDOF = new ID(size);
	int result = theChannel.recvID(dbTag2, cTag, *retainDOF);
	if (result < 0) {
	    opserr << "WARNING MP_Constraint::recvSelf ";
	    opserr << "- error receiving retained data\n"; 
	    return result;  
	}	
    }    
    
    // recv initial displacement vectors.
    // we need 2 database tags because they have the same size,
    // but we can reuse the tags used for ID objects, since they go into different files
    if (constrDOF && constrDOF->Size() > 0)
        Uc0.resize(constrDOF->Size());
    else
        Uc0 = Vector();
    if (retainDOF && retainDOF->Size() > 0)
        Ur0.resize(retainDOF->Size());
    else
        Ur0 = Vector();
    if (Uc0.Size() > 0) {
        int result = theChannel.recvVector(dbTag1, cTag, Uc0);
        if (result < 0) {
            opserr << "WARNING MP_Constraint::recvSelf ";
            opserr << "- error receiving constrained initial displacement\n";
            return result;
        }
    }
    if (Ur0.Size() > 0) {
        int result = theChannel.recvVector(dbTag2, cTag, Ur0);
        if (result < 0) {
            opserr << "WARNING MP_Constraint::recvSelf ";
            opserr << "- error receiving retained initial displacement\n";
            return result;
        }
    }

    return 0;
}



void
MP_Constraint::Print(OPS_Stream &s, int flag)
{     

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
      const char *indent = " ";
      const char *newln = " ";
      s << OPS_PRINT_JSON_ELEM_INDENT << "{";
      s << "\"name\": \"" << this->getTag() << "\"," << newln;
      s << indent << "\"node\": " << nodeConstrained << "," << newln;
      s << indent << "\"ref_node\": " << nodeRetained << "," << newln;
      if (constrDOF != 0 && retainDOF != 0) {
        s << indent << "\"dofs\": [";
        const int nc = (*constrDOF).Size();
        for (int i=0; i < nc; i++)
          s << (*constrDOF)(i)+1 << (i < nc-1? ", " : ""); 
        s << "]," << newln;

        const int nr = (*retainDOF).Size();
        s << indent << "\"ref_dof\": [";
        for (int i=0; i < nr; i++)
          s << (*retainDOF)(i)+1 << (i < nr-1? ", " : "");
        s << "]," << newln;

        if (constraint != 0) {
          s << indent << "\"matrix\": [";
          for (int i=0; i<constraint->noRows(); i++) {
            s << "[";
            for (int j=0; j<constraint->noCols(); j++) {
              s << (*constraint)(i,j);
              if (j < constraint->noCols()-1)
                s << ", ";
            }
            s << "]";
            if (i < constraint->noRows()-1)
              s << ", ";
          }
          s << "]";
        }
      }
      s << "}";
      return;

    } else {
      s << "MP_Constraint: " << this->getTag() << "\n";
      s << "\tNode Constrained: " << nodeConstrained;
      s << " node Retained: " << nodeRetained << "\n";
      if (constrDOF != 0 && retainDOF != 0) {
        s << " constrained dof: ";
        for (int i=0; i<(*constrDOF).Size(); i++)
          s << (*constrDOF)(i)+1 << " ";
        s << "\n";
          s << " retained dof: ";        
        for (int i=0; i<(*retainDOF).Size(); i++)
          s << (*retainDOF)(i)+1 << " ";
        s << "\n";
        if (constraint != 0)
          s << " constraint matrix: " << *constraint << "\n";
      }
    }
}


