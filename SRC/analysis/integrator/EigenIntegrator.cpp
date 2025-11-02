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
                                                                        
// $Revision: 1.3 $
// $Date: 2005-12-19 22:43:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/EigenIntegrator.cpp,v $
                                                                        
// Written: Jun Peng
// Created: Wed Jan 27, 1999
// Revision: A
//
// Description: This file contains the class definition of EigenIntegrator.
// EigenIntegrator is an algorithmic class for setting up the finite element 
// equations for a eigen problem analysis.
//
// This class is inheritanted from the base class of Integrator which was
// created by fmk (Frank).


#include "EigenIntegrator.h"
#include <FE_Element.h>
#include <AnalysisModel.h>
#include <EigenSOE.h>
#include <Vector.h>
#include <DOF_Group.h>
#include <FE_EleIter.h>
#include <DOF_GrpIter.h>

EigenIntegrator::EigenIntegrator(AnalysisModel &theModel, EigenSOE &theSysOE)
  : theSOE(&theSysOE), theAnalysisModel(&theModel)
{

}

EigenIntegrator::~EigenIntegrator()
{

}


int 
EigenIntegrator::formEleTangent(FE_Element *theEle)
{
  if (flagK == 0)
    return this->formEleTangK(theEle);
  else
    return this->formEleTangM(theEle);
}

int 
EigenIntegrator::formNodTangent(DOF_Group *theDof)
{
  return this->formNodTangM(theDof);
}

int 
EigenIntegrator::formEleResidual(FE_Element *)
{
  return 0;
}

int 
EigenIntegrator::formNodUnbalance(DOF_Group *theDof)
{
  return 0;
}


int 
EigenIntegrator::getLastResponse(Vector &result, const ID &id)
{
  return 0;
}

int
EigenIntegrator::formK()
{    
  // the loops to form and add the tangents are broken into two for 
  // efficiency when performing parallel computations

  // loop through the FE_Elements getting them to form the tangent
  // FE_EleIter &theEles1 = theAnalysisModel->getFEs();
  FE_Element *elePtr;

  flagK = 0;

  theSOE->zeroA();

  //while((elePtr = theEles1()) != 0) 
  //  elePtr->formTangent(this);
  
  // loop through the FE_Elements getting them to add the tangent    
  int result = 0;
  FE_EleIter &theEles2 = theAnalysisModel->getFEs();    
  while((elePtr = theEles2()) != 0) {
    if (theSOE->addA(elePtr->getTangent(this), elePtr->getID()) < 0) {
        opserr << "WARNING EigenIntegrator::formK -";
        opserr << " failed in addA for ID " << elePtr->getID();	    
        result = -2;
    }
  }

  return result;
}


int
EigenIntegrator::formM()
{  
  // the loops to form and add the tangents are broken into two for 
  // efficiency when performing parallel computations

  // loop through the FE_Elements getting them to form the tangent
  // FE_EleIter &theEles1 = theAnalysisModel->getFEs();

  flagK = 1;
  theSOE->zeroM();

  // while((elePtr = theEles1()) != 0) 
  //     elePtr->formTangent(this);
  
  // loop through the FE_Elements getting them to add the tangent    
  int result = 0;
  FE_Element *elePtr;
  FE_EleIter &theEles2 = theAnalysisModel->getFEs();    
  while((elePtr = theEles2()) != nullptr) {     
    if (theSOE->addM(elePtr->getTangent(this), elePtr->getID()) < 0) {
      opserr << "WARNING EigenIntegrator::formK -";
      opserr << " failed in addA for ID " << elePtr->getID();	    
      result = -2;
    }
  }

  DOF_Group *dofPtr;
  DOF_GrpIter &theDofs = theAnalysisModel->getDOFs();    
  while((dofPtr = theDofs()) != nullptr) {
    //   	dofPtr->formTangent(this);
    if (theSOE->addM(dofPtr->getTangent(this),dofPtr->getID()) < 0) {
      opserr << "WARNING EigenIntegrator::formM -";
      opserr << " failed in addM for ID " << dofPtr->getID();	    
      result = -3;
    }
  }

  return result;
}

int 
EigenIntegrator::formEleTangK(FE_Element *theEle)
{
  theEle->zeroTangent();
  theEle->addKtToTang(1.0);
  return 0;
}

int 
EigenIntegrator::formEleTangM(FE_Element *theEle)
{
  theEle->zeroTangent();
  theEle->addMtoTang(1.0);
  return 0;
}

int 
EigenIntegrator::formNodTangM(DOF_Group *theDof)
{
  theDof->zeroTangent();
  theDof->addMtoTang(1.0);
  return 0;
}
