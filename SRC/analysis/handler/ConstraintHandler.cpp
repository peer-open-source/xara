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
                                                                        
// $Revision: 1.4 $
// $Date: 2005-11-28 21:34:03 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/handler/ConstraintHandler.cpp,v $
                                                                        
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Description: This file contains the implementation of ConstraintHandler.
//
// What: "@(#) ConstraintHandler.h, revA"

#include <ConstraintHandler.h>
#include <Domain.h>
#include <AnalysisModel.h>
#include <FE_EleIter.h>
#include <FE_Element.h>
#include <DOF_GrpIter.h>
#include <DOF_Group.h>

ConstraintHandler::ConstraintHandler(int clasTag)
 : MovableObject(clasTag)
 , theAnalysisModelPtr(nullptr)
{

}


ConstraintHandler::~ConstraintHandler()
{
    
}


int
ConstraintHandler::doneNumberingDOF()
{

  // iterate through the DOF_Groups telling them that their ID has now been set
  AnalysisModel *theModel1=this->getAnalysisModelPtr();
  DOF_GrpIter &theDOFS = theModel1->getDOFs();
  DOF_Group *dofPtr;
  while ((dofPtr = theDOFS()) != nullptr) {
    dofPtr->doneID();
  }


  // iterate through the FE_Element getting them to set their IDs
  AnalysisModel *theModel=this->getAnalysisModelPtr();
  FE_EleIter &theEle = theModel->getFEs();
  FE_Element *elePtr;
  while ((elePtr = theEle()) != nullptr) {
    elePtr->setID();
  }

  return 0;
}


void 
ConstraintHandler::setLinks(AnalysisModel &theModel)
{
  theAnalysisModelPtr = &theModel;
}
	

int
ConstraintHandler::update()
{
  return 0;
}

int
ConstraintHandler::applyLoad()
{
  return 0;
}


AnalysisModel *
ConstraintHandler::getAnalysisModelPtr() const
{
  return theAnalysisModelPtr;
}
