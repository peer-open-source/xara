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
// Description: This file contains the implementation of PlainHandler.
//
// Written: fmk 
// Revision: A
//
#include <PlainHandler.h>
#include <stdlib.h>
#include <Matrix.h>
#include <AnalysisModel.h>
#include <Domain.h>
#include <FE_Element.h>
#include <DOF_Group.h>
#include <Node.h>
#include <Element.h>
#include <NodeIter.h>
#include <ElementIter.h>
#include <SP_ConstraintIter.h>
#include <SP_Constraint.h>
#include <MP_ConstraintIter.h>
#include <MP_Constraint.h>
#include <ID.h>
#include <Subdomain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <map>

PlainHandler::PlainHandler()
 : ConstraintHandler(HANDLER_TAG_PlainHandler)
{

}

PlainHandler::~PlainHandler()
{

}

int
PlainHandler::handle(const ID *nodesLast)
{
  // first check links exist to a Domain and an AnalysisModel
  AnalysisModel *theModel = this->getAnalysisModelPtr();
  Domain *theDomain = theModel->getDomainPtr();

  // get all SPs
  std::multimap<int,SP_Constraint*> allSPs;
  SP_ConstraintIter &theSPs = theDomain->getDomainAndLoadPatternSPs();
  SP_Constraint *theSP; 
  while ((theSP = theSPs()) != nullptr) {
    if (theSP->isHomogeneous() == false) {
      opserr << "WARNING PlainHandler::handle() - ";
      opserr << " non-homogeneos constraint";
      opserr << " for node " << theSP->getNodeTag();
      opserr << " homogeneous constraint assumed\n";
    }
    allSPs.insert(std::make_pair(theSP->getNodeTag(),theSP));
  }

  // initialise the DOF_Groups and add them to the AnalysisModel.
  //    : must of course set the initial IDs
  NodeIter &theNod = theDomain->getNodes();
  Node *nodPtr;
  SP_Constraint *spPtr;

  int numDOF = 0;
  int count3 = 0;
  int countDOF =0;
  while ((nodPtr = theNod()) != nullptr) {
    DOF_Group *dofPtr;
    if ((dofPtr = new DOF_Group(numDOF++, nodPtr)) == 0) {
      opserr << "WARNING PlainHandler::handle() - ran out of memory";
      opserr << " creating DOF_Group " << numDOF << endln;	
      return -4;    		
    }
    // initially set all the ID value to -2
    const ID &id = dofPtr->getID();
    for (int j=0; j < id.Size(); j++) {
      dofPtr->setID(j,-2);
      countDOF++;
    }
    // loop through the SP_Constraints to see if any of the
    // DOFs are constrained, if so set initial ID value to -1
    int nodeID = nodPtr->getTag();
    std::multimap<int,SP_Constraint*>::iterator first = allSPs.lower_bound(nodeID);
    std::multimap<int,SP_Constraint*>::iterator last = allSPs.upper_bound(nodeID);
    for (auto it=first; it!=last; it++) {
	    spPtr = it->second;
	    const ID &id = dofPtr->getID();
	    int dof = spPtr->getDOF_Number();		
	    if (id(dof) == -2) {
        dofPtr->setID(spPtr->getDOF_Number(),-1);
        countDOF--;	
      } else {
        opserr << "WARNING PlainHandler::handle() - ";
        opserr << " multiple single pointconstraints at DOF " << dof;
        opserr << " for node " << spPtr->getNodeTag() << endln;
      }
    }

    // loop through the MP_Constraints to see if any of the
    // DOFs are constrained, note constraint matrix must be diagonal
    // with 1's on the diagonal
    MP_ConstraintIter &theMPs = theDomain->getMPs();
    MP_Constraint *mpPtr;
    while ((mpPtr = theMPs()) != nullptr) {
	    if (mpPtr->getNodeConstrained() == nodeID) {
        if (mpPtr->isTimeVarying() == true) {
          opserr << "WARNING PlainHandler::handle() - ";
          opserr << " time-varying constraint";
          opserr << " for node " << nodeID;
          opserr << " non-varyng assumed\n";
        }
        const Matrix &C = mpPtr->getConstraint();
        int numRows = C.noRows();
        int numCols = C.noCols();
        if (numRows != numCols) {
          opserr << "WARNING PlainHandler::handle() - ";
          opserr << " constraint matrix not diagonal, ignoring constraint";
          opserr << " for node " << nodeID << endln;
          opserr << " non-varyng assumed\n";
        } else {
          int ok = 0;
          for (int i=0; i<numRows; i++) {
            if (C(i,i) != 1.0)
              ok = 1;
            for (int j=0; j<numRows; j++)
              if (i != j)
                if (C(i,j) != 0.0)
                  ok = 1;
          }
          if (ok != 0) {
            opserr << "WARNING PlainHandler::handle() - ";
            opserr << " constraint matrix not identity, ignoring constraint";
            opserr << " for node " << nodeID << endln;
            opserr << " non-varyng assumed\n";
          } else {
            const ID &dofs = mpPtr->getConstrainedDOFs();
            const ID &id = dofPtr->getID();				
            for (int i=0; i<dofs.Size(); i++) {
              int dof = dofs(i);	
              if (id(dof) == -2) {
                dofPtr->setID(dof,-4);
                countDOF--;	
              } else {
                opserr << "WARNING PlainHandler::handle() - ";
                opserr << " constraint at dof " << dof << " already specified for constrained node";
                opserr << " in MP_Constraint at node " << nodeID << endln;
              }		      
            }
          }
        }
      }
    }

    nodPtr->setDOF_GroupPtr(dofPtr);
    theModel->addDOF_Group(dofPtr);
  }

  // set the number of eqn in the model
  theModel->setNumEqn(countDOF);

  // now see if we have to set any of the dof's to -3
  //    int numLast = 0;
  if (nodesLast != 0) 
    for (int i=0; i<nodesLast->Size(); i++) {
        int nodeID = (*nodesLast)(i);
        Node *nodPtr = theDomain->getNode(nodeID);
        if (nodPtr != nullptr) {
          DOF_Group* dofPtr = nodPtr->getDOF_GroupPtr();
          
          const ID &id = dofPtr->getID();
          // set all the dof values to -3
          for (int j=0; j < id.Size(); j++) 
            if (id(j) == -2) {
              dofPtr->setID(j,-3);
              count3++;
            } else {
              opserr << "WARNING PlainHandler::handle() ";
              opserr << " - boundary sp constraint in subdomain";
              opserr << " this should not be - results suspect \n";
            }
        }
    }

  // initialise the FE_Elements and add to the AnalysisModel.
  ElementIter &theEle = theDomain->getElements();
  Element *elePtr;

  int numFe = 0;
  while ((elePtr = theEle()) != nullptr) {

    // only create an FE_Element for a subdomain element if it does not
    // do independent analysis .. then subdomain part of this analysis so create
    // an FE_element & set subdomain to point to it.
    if (elePtr->isSubdomain() == true) {

      Subdomain *theSub = (Subdomain *)elePtr;
      if (theSub->doesIndependentAnalysis() == false) {

        FE_Element *fePtr = new FE_Element(numFe++, elePtr);
        theModel->addFE_Element(fePtr);
        theSub->setFE_ElementPtr(fePtr);
      }

    } else {
      // just a regular element .. create an FE_Element for it & add to AnalysisModel
      theModel->addFE_Element(new FE_Element(numFe++, elePtr));
    }
  }

  return count3;
}
 

int
PlainHandler::sendSelf(int cTag, Channel &theChannel)
{
  return 0;
}

int
PlainHandler::recvSelf(int cTag, 
		       Channel &theChannel, 
		       FEM_ObjectBroker &theBroker)
{
  return 0;
}
