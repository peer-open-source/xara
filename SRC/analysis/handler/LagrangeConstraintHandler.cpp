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

// $Revision: 1.5 $
// $Date: 2005-11-29 22:04:40 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/handler/LagrangeConstraintHandler.cpp,v $

// Written: fmk 
// Created: May 1998
//
#include <LagrangeConstraintHandler.h>
#include <stdlib.h>

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
#include <LagrangeDOF_Group.h>
#include <LagrangeSP_FE.h>
#include <LagrangeMP_FE.h>


LagrangeConstraintHandler::LagrangeConstraintHandler(double sp, double mp)
:ConstraintHandler(HANDLER_TAG_LagrangeConstraintHandler),
 alphaSP(sp), alphaMP(mp)
{

}

LagrangeConstraintHandler::~LagrangeConstraintHandler()
{

}

int
LagrangeConstraintHandler::handle(const ID *nodesLast)
{
  // first check links exist to a Domain and an AnalysisModel object
  AnalysisModel *theModel = this->getAnalysisModelPtr();
  Domain *theDomain = theModel->getDomainPtr();

  // get number of elements and nodes in the domain 
  // and init the theFEs and theDOFs arrays

  int numConstraints = 0;
  SP_ConstraintIter &theSPss = theDomain->getDomainAndLoadPatternSPs();
  SP_Constraint *spPtr;
  while ((spPtr = theSPss()) != nullptr)
    numConstraints++;

  numConstraints += theDomain->getNumMPs();

  // create a DOF_Group for each Node
  //    : must of course set the initial IDs
  NodeIter &theNod = theDomain->getNodes();
  Node *nodPtr;
  MP_Constraint *mpPtr;    
  DOF_Group *dofPtr;
  
  int numDofGrp = 0;
  int count3 = 0;
  int countDOF =0;
  while ((nodPtr = theNod()) != nullptr) {
    dofPtr = new DOF_Group(numDofGrp++, nodPtr);

    // initially set all the ID value to -2
    
    const ID &id = dofPtr->getID();
    for (int j=0; j < id.Size(); j++) {
      dofPtr->setID(j,-2);
      countDOF++;
    }

    nodPtr->setDOF_GroupPtr(dofPtr);
    theModel->addDOF_Group(dofPtr);
  }

  // create the FE_Elements for the Elements and add to the AnalysisModel
  ElementIter &theEle = theDomain->getElements();
  Element *elePtr;

  int numFeEle = 0;
  while ((elePtr = theEle()) != nullptr) {

    // only create an FE_Element for a subdomain element if it does not
    // do independent analysis .. then subdomain part of this analysis so create
    // an FE_element & set subdomain to point to it.
    if (elePtr->isSubdomain() == true) {
      Subdomain *theSub = (Subdomain *)elePtr;
      if (theSub->doesIndependentAnalysis() == false) {
        FE_Element *fePtr = new FE_Element(numFeEle++, elePtr);
        theModel->addFE_Element(fePtr);
        theSub->setFE_ElementPtr(fePtr);
      }
    }
    else {
      // just a regular element
      theModel->addFE_Element(new FE_Element(numFeEle++, elePtr));
    }
  }

  // create the LagrangeSP_FE for the SP_Constraints and 
  // add to the AnalysisModel
  SP_ConstraintIter &theSPs = theDomain->getDomainAndLoadPatternSPs();
  while ((spPtr = theSPs()) != nullptr) {
    dofPtr = new LagrangeDOF_Group(numDofGrp++, *spPtr);
    const ID &id = dofPtr->getID();
    for (int j=0; j < id.Size(); j++) {
      dofPtr->setID(j,-2);
      countDOF++;
    }
    // Add the DOF_Group to the model
    theModel->addDOF_Group(dofPtr);    		

    // Create the FE
    theModel->addFE_Element(new LagrangeSP_FE(
        numFeEle++, *theDomain, *spPtr, *dofPtr, alphaSP
    ));
  }	    

  // create the LagrangeMP_FE for the MP_Constraints and 
  // add to the AnalysisModel    

  MP_ConstraintIter &theMPs = theDomain->getMPs();
  while ((mpPtr = theMPs()) != nullptr) {
      dofPtr = new LagrangeDOF_Group(numDofGrp++, *mpPtr);

      const ID &id = dofPtr->getID();
      for (int j=0; j < id.Size(); j++) {
        dofPtr->setID(j,-2);
        countDOF++;
      }

      theModel->addDOF_Group(dofPtr);

      theModel->addFE_Element(new LagrangeMP_FE(
        numFeEle++, *theDomain, *mpPtr, *dofPtr, alphaMP
      ));
  }
  
  theModel->setNumEqn(countDOF);

  // set the number of eqn in the model
  // now see if we have to set any of the dof's to -3
  //    int numLast = 0;
  if (nodesLast != nullptr) 
    for (int i=0; i<nodesLast->Size(); i++) {
      int nodeID = (*nodesLast)(i);
      Node *nodPtr = theDomain->getNode(nodeID);
      if (nodPtr != nullptr) {
        dofPtr = nodPtr->getDOF_GroupPtr();
        
        const ID &id = dofPtr->getID();
        // set all the dof values to -3
        for (int j=0; j < id.Size(); j++) {
          if (id(j) == -2) {
            dofPtr->setID(j,-3);
            count3++;
          } else {
            opserr << "WARNING LagrangeConstraintHandler::handle() ";
            opserr << " - boundary sp constraint in subdomain";
            opserr << " this should not be - results suspect \n";
          }
        }
      }
    }

  return count3;
}


int
LagrangeConstraintHandler::sendSelf(int cTag, Channel &theChannel)
{
  Vector data(2);
  int result = 0;
  data(0) = alphaSP;
  data(1) = alphaMP;
  result = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (result != 0) 
    opserr << "LagrangeConstraintHandler::sendSelf() - error sending Vector\n";
  return result;
}

int
LagrangeConstraintHandler::recvSelf(int cTag, 
				   Channel &theChannel, 
				   FEM_ObjectBroker &theBroker)  
{
  Vector data(2);
  int result = 0;
  result = theChannel.recvVector(this->getDbTag(), cTag, data);
  alphaSP = data(0);
  alphaMP = data(1);
  if (result != 0) 
    opserr << "LagrangeConstraintHandler::recvSelf() - error receiving Vector\n";
  return result;
}

