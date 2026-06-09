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
// Written: fmk 11/98
// Revised:
//
// Purpose: This file contains the class definition for EarthquakePattern.
// EarthquakePattern is an abstract class.

#include <EarthquakePattern.h>
#include <GroundMotion.h>

#include <Domain.h>
#include <NodeIter.h>
#include <Node.h>
#include <ElementIter.h>
#include <Element.h>
#include <stdlib.h>
#include <Channel.h>
#include <ErrorHandler.h>

#include <string.h>
#include <stdlib.h>

EarthquakePattern::EarthquakePattern(int tag, int _classTag)
  :LoadPattern(tag, _classTag, 1.0)
  , theMotions(0), numMotions(0), uDotG(0), uDotDotG(0), currentTime(0.0), parameterID(0)
{

}


EarthquakePattern::~EarthquakePattern()
{
  // invoke the destructor on all ground motions supplied
  for (int i=0; i<numMotions; i++)
    delete theMotions[i];
  
  if (theMotions != 0)
    delete [] theMotions;

  if (uDotG != 0)
    delete uDotG;

  if (uDotDotG != 0)
    delete uDotDotG;
}



void 
EarthquakePattern::applyLoad(double time)
{
  // see if quick return, i.e. no Ground Motions or domain set
  if (numMotions == 0)
    return;

  // check if setLoadConstant() has been called
  if (isConstant != true)
    currentTime = time;

  Domain *theDomain = this->getDomain();
  if (theDomain == nullptr)
    return;

#if 0
  // set the vel and accel vector
  for (int i=0; i<numMotions; i++)
    (*uDotDotG)(i) = theMotions[i]->getAccel(currentTime);

  NodeIter &theNodes = theDomain->getNodes();
  Node *theNode;
  while ((theNode = theNodes()) != nullptr)
    theNode->addInertiaLoadToUnbalance(*uDotDotG, 1.0);


  ElementIter &theElements = theDomain->getElements();
  Element *theElement;
  while ((theElement = theElements()) != nullptr) 
    theElement->addInertiaLoadToUnbalance(*uDotDotG);
#endif
}


void 
EarthquakePattern::applyLoadSensitivity(double time)
{
  // see if quick return, i.e. no Ground Motions or domain set
  if (numMotions == 0)
    return;

  Domain *theDomain = this->getDomain();
  if (theDomain == nullptr)
    return;


  // set the vel and accel vector
  for (int i=0; i<numMotions; i++) {
    (*uDotG)(i) = theMotions[i]->getVel(time);
    if (parameterID != 0) { // Something is random in the motions
      (*uDotDotG)(i) = theMotions[i]->getAccelSensitivity(time);
    }
    else {
      (*uDotDotG)(i) = theMotions[i]->getAccel(time);
    }
  }

  bool somethingRandomInMotions = false;
  if (parameterID != 0) {
    somethingRandomInMotions = true;
  }


  NodeIter &theNodes = theDomain->getNodes();
  Node *theNode;
  while ((theNode = theNodes()) != 0) 
    theNode->addInertiaLoadSensitivityToUnbalance(*uDotDotG, 1.0,  somethingRandomInMotions);


  ElementIter &theElements = theDomain->getElements();
  Element *theElement;
  while ((theElement = theElements()) != 0) 
    theElement->addInertiaLoadSensitivityToUnbalance(*uDotDotG,  somethingRandomInMotions);
}


int
EarthquakePattern::addMotion(GroundMotion &theMotion)
{
  // make space for new
  GroundMotion **newMotions = new GroundMotion *[numMotions+1];
  //  GroundMotion **newMotions = (GroundMotion **)malloc(sizeof(GroundMotion *)*(numMotions+1));
  if (newMotions == 0) {
    opserr << "EarthquakePattern::addMotion - could not add new, out of mem\n";
    return -1;
  }
  
  // copy old
  for (int i=0; i<numMotions; i++)
    newMotions[i] = theMotions[i];

  // add the new motion to new
  newMotions[numMotions] = &theMotion;

  // delete the old
  if (theMotions != 0)
    delete [] theMotions;

  // reset
  theMotions = newMotions;
  numMotions++;

  // delete old vectors and create new ones
  if (uDotG != 0)
    delete uDotG;
  uDotG = new Vector(numMotions);

  if (uDotDotG != 0)
    delete uDotDotG;
  uDotDotG = new Vector(numMotions);

  return 0;
}


bool
EarthquakePattern::addSP_Constraint(SP_Constraint *)
{
  return false;
}

bool
EarthquakePattern::addNodalLoad(NodalLoad *)
{
  return false;
}

bool
EarthquakePattern::addElementalLoad(ElementalLoad *)
{
  return false;
}

// void 
// EarthquakePattern::Print(OPS_Stream &s, int flag)
// {   

// }


// AddingSensitivity:BEGIN ////////////////////////////////////
// int
// EarthquakePattern::setParameter(const char **argv, int argc, Parameter &param)
// {
//   if (argc < 2)
//     return -1;

//   // This is not general, yet

//   if (strstr(argv[0],"randomProcessDiscretizer") != 0) {
//     return theMotions[0]->setParameter(&argv[1], argc-1, param);
//   }
//   else
//     return 0;
// }

// int
// EarthquakePattern::updateParameter(int pparameterID, Information &info)
// {
//   return theMotions[0]->updateParameter(pparameterID, info);
// }

// int
// EarthquakePattern::activateParameter(int pparameterID)
// {
//   parameterID = pparameterID;

//   return theMotions[0]->activateParameter(pparameterID);
// }
// AddingSensitivity:END ////////////////////////////////////
