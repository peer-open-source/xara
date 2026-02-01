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

FE_Element::FE_Element(int tag, int numDOF_Group, int ndof)
  :TaggedObject(tag),
   myDOF_Groups(numDOF_Group), numDOF(ndof), theModel(nullptr),
   theResidual(nullptr), theTangent(nullptr), theIntegrator(nullptr)
{
  // this is for a subtype, the subtype must set the myDOF_Groups ID array
  // as subtypes have no access to the tangent or residual we don't set them
  // this way we can detect if subclass does not provide all methods it should
}



FE_Element::~FE_Element()
{
  // delete tangent and residual if created specially
  if (numDOF > MAX_NUM_DOF) {
    if (theTangent != nullptr)
      delete theTangent;
    if (theResidual != nullptr) 
      delete theResidual;
  }
}


const ID &
FE_Element::getDOFtags() const
{
  return myDOF_Groups;
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
  assert(theModel != nullptr);
  if (this->setID(*theModel) == 0) {
    return 0;
  }
  return -1;
}


void
FE_Element::zeroTangent()
{
}


void
FE_Element::addKtToTang(double fact)
{
}


void
FE_Element::addCtoTang(double fact)
{
}


void
FE_Element::addMtoTang(double fact)
{
}


void
FE_Element::addKiToTang(double fact)
{
}


void
FE_Element::addKpToTang(double fact, int numP)
{
}

int
FE_Element::storePreviousK(int numP)
{
}

//
// RESIDUAL
//

void
FE_Element::zeroResidual()
{
}


void
FE_Element::addRtoResidual(double fact)
{
}


void
FE_Element::addRIncInertiaToResidual(double fact)
{
}


const Vector &
FE_Element::getTangForce(const Vector &disp, double fact)
{
  return *theResidual;
}



const Vector &
FE_Element::getK_Force(const Vector &disp, double fact)
{
  return *theResidual;
}


const Vector &
FE_Element::getKi_Force(const Vector &disp, double fact)
{
  return *theResidual;
}

const Vector &
FE_Element::getM_Force(const Vector &disp, double fact)
{
  return *theResidual;
}

const Vector &
FE_Element::getC_Force(const Vector &disp, double fact)
{
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
  return *theResidual;
}


void
FE_Element::addD_Force(const Vector &accel, double fact)
{
}

void
FE_Element::addK_Force(const Vector &disp, double fact)
{
}



// AddingSensitivity:BEGIN /////////////////////////////////
void
FE_Element::addResistingForceSensitivity(int gradNumber, double fact)
{
  // theResidual->addVector(1.0, myEle->getResistingForceSensitivity(gradNumber), -fact);
}

void
FE_Element::addM_ForceSensitivity(int gradNumber, const Vector &vect, double fact)
{
}

void
FE_Element::addD_ForceSensitivity(int gradNumber, const Vector &vect, double fact)
{
}


int
FE_Element::commitSensitivity(int gradNum, int numGrads)
{
  return 0;
}

// AddingSensitivity:END ////////////////////////////////////

