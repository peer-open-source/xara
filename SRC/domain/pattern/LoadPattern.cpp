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
// Purpose: This file contains the method definitions for class LoadPattern.
// LoadPattern is a container class.
//
// Written: fmk 07/99
// Revised:
//
#include <string.h>

#include <LoadPattern.h>
#include <ID.h>
#include <TimeSeries.h>
#include <NodalLoad.h>
#include <ElementalLoad.h>
#include <SP_Constraint.h>
#include <MapOfTaggedObjects.h>
#include <ElementalLoadIter.h>
#include <NodalLoadIter.h>
#include <SingleDomSP_Iter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <OPS_Globals.h>


LoadPattern::LoadPattern(int tag, int clasTag, double fact)
  : TaggedObject(tag), MovableObject(clasTag), 
    isConstant(false), 
    loadFactor(0.0), scaleFactor(fact), 
    theSeries(nullptr), 
    theDomain(nullptr),
    currentGeoTag(0), lastGeoSendTag(-1),
    theNodalLoads(nullptr), theElementalLoads(nullptr), theSPs(nullptr), theNodIter(nullptr),
    theEleIter(nullptr), theSpIter(nullptr), lastChannel(0)
{
  // constructor for subclass
  theNodalLoads     = new MapOfTaggedObjects();
  theElementalLoads = new MapOfTaggedObjects();
  theSPs            = new MapOfTaggedObjects();

  theEleIter = new ElementalLoadIter(theElementalLoads);
  theNodIter = new NodalLoadIter(theNodalLoads);
  theSpIter  = new SingleDomSP_Iter(theSPs);

  randomLoads = 0;
  dLambdadh   = 0;
}

LoadPattern::LoadPattern()
  : TaggedObject(0), MovableObject(PATTERN_TAG_LoadPattern),
  isConstant(false),
  loadFactor(0.0), scaleFactor(1.0), 
  theSeries(nullptr),
  theDomain(nullptr),
  currentGeoTag(0), lastGeoSendTag(-1), 
  dbSPs(0), dbNod(0), dbEle(0), theNodalLoads(0),
  theElementalLoads(0), 
  theSPs(0), 
  theNodIter(0), 
  theEleIter(0), 
  theSpIter(0),
  lastChannel(0)
{
  theNodalLoads     = new MapOfTaggedObjects();
  theElementalLoads = new MapOfTaggedObjects();
  theSPs            = new MapOfTaggedObjects();

  theEleIter = new ElementalLoadIter(theElementalLoads);
  theNodIter = new NodalLoadIter(theNodalLoads);
  theSpIter  = new SingleDomSP_Iter(theSPs);

  randomLoads = 0;
  dLambdadh   = 0;
}


LoadPattern::~LoadPattern()
{
  if (theSeries != 0)
    delete theSeries;

  if (theNodalLoads != 0)
    delete theNodalLoads;

  if (theElementalLoads != 0)
    delete theElementalLoads;

  if (theSPs != 0)
    delete theSPs;

  if (theEleIter != 0)
    delete theEleIter;

  if (theNodIter != 0)
    delete theNodIter;

  if (theSpIter != 0)
    delete theSpIter;

  if (randomLoads != 0)
    delete randomLoads;

  if (dLambdadh != 0)
    delete dLambdadh;
}


void
LoadPattern::setDomain(Domain *theDomain)
{
  // if subclass does not implement .. check for 0 pointer
  if (theNodalLoads != nullptr) {
    // TODO: why sps and eles are guarded??


#ifdef OLD_LOAD_PATTERN
    NodalLoad *nodLoad;
    NodalLoadIter &theNodalIter = this->getNodalLoads();
    while ((nodLoad = theNodalIter()) != nullptr)
      nodLoad->setDomain(theDomain);

    ElementalLoad *eleLoad;
    ElementalLoadIter &theElementalIter = this->getElementalLoads();
    while ((eleLoad = theElementalIter()) != nullptr)
      eleLoad->setDomain(theDomain);
#endif
    SP_Constraint *theSP;
    SP_ConstraintIter &theSpConstraints = this->getSPs();
    while ((theSP = theSpConstraints()) != nullptr)
      theSP->setDomain(theDomain);
  }

  this->theDomain = theDomain;
}


bool
LoadPattern::addSP_Constraint(SP_Constraint *theSp)
{
  Domain *theDomain = this->getDomain();

  bool result = theSPs->addComponent(theSp);
  if (result == true) {
    if (theDomain != 0)
      theSp->setDomain(theDomain);
    theSp->setLoadPatternTag(this->getTag());
    currentGeoTag++;
  } else
    opserr << "WARNING: LoadPattern::addSP_Constraint() - load could not be "
              "added\n";
  return result;
}

SP_ConstraintIter &
LoadPattern::getSPs()
{
  theSpIter->reset();
  return *theSpIter;
}


void
LoadPattern::clearAll()
{
  theElementalLoads->clearAll();
  theNodalLoads->clearAll();
  theSPs->clearAll();
  currentGeoTag++;
  lastChannel = 0;
  if (dLambdadh != 0) {
    dLambdadh->Zero();
  }
}


SP_Constraint *
LoadPattern::removeSP_Constraint(int tag)
{
  TaggedObject *obj = theSPs->removeComponent(tag);
  if (obj == nullptr)
    return nullptr;
  SP_Constraint *result = (SP_Constraint *)obj;
  result->setDomain(nullptr);
  currentGeoTag++;
  return result;
}



void
LoadPattern::applyLoad(double pseudoTime)
{

#ifdef OLD_LOAD_PATTERN
  // first determine the load factor
  if (theSeries != nullptr && isConstant != true) {
    loadFactor = theSeries->getFactor(pseudoTime);
    loadFactor *= scaleFactor;
  }
  {
    Load *nodLoad;
    NodalLoadIter &theNodalIter = this->getNodalLoads();
    while ((nodLoad = theNodalIter()) != nullptr)
      nodLoad->applyLoad(loadFactor);
  }

  {
    Load *eleLoad;
    ElementalLoadIter &theElementalIter = this->getElementalLoads();
    while ((eleLoad = theElementalIter()) != nullptr)
      eleLoad->applyLoad(loadFactor);
  }
#else
  double loadFactor = this->getLoadFactor();
#endif
  SP_Constraint *sp;
  SP_ConstraintIter &theIter = this->getSPs();
  while ((sp = theIter()) != nullptr)
    sp->applyConstraint(loadFactor);
}

void LoadPattern::setLoadConstant() { isConstant = true; }

void LoadPattern::unsetLoadConstant() { isConstant = false; }


double
LoadPattern::getLoadFactor()
{
  if (theSeries != nullptr)
    return loadFactor;
  else
    return 0.0;
}


int
LoadPattern::sendSelf(int cTag, Channel &theChannel)
{
  return -1;
}

int 
LoadPattern::recvSelf(int cTag, Channel &theChannel,
                          FEM_ObjectBroker &theBroker)
{
  return -1;
}

#if 0
void 
LoadPattern::applyLoadSensitivity(double pseudoTime)
{
  // P*dfactor/dh
  if (theSeries != 0 && isConstant != true) {
    loadFactor = theSeries->getFactorSensitivity(pseudoTime);
    loadFactor *= scaleFactor;
  }

  NodalLoad *nodLoad;
  NodalLoadIter &theNodalIter = this->getNodalLoads();
  while ((nodLoad = theNodalIter()) != 0)
    nodLoad->applyLoad(loadFactor);

  // factor*dP/dh
  if (theSeries != 0 && isConstant != true) {
    loadFactor = theSeries->getFactor(pseudoTime);
    loadFactor *= scaleFactor;
  }

  NodalLoadIter &theNodalIter2 = this->getNodalLoads();
  while ((nodLoad = theNodalIter2()) != 0)
    nodLoad->applyLoadSensitivity(loadFactor);

  // Don't include element loads and sp constraints for now
  /*
    ElementalLoad *eleLoad;
    ElementalLoadIter &theElementalIter = this->getElementalLoads();
    while ((eleLoad = theElementalIter()) != 0)
    eleLoad->applyLoad(loadFactor);
    
    SP_Constraint *sp;
    SP_ConstraintIter &theIter = this->getSPs();
    while ((sp = theIter()) != 0)
    sp->applyConstraint(loadFactor);
  */
}


int 
LoadPattern::setParameter(const char **argv, int argc, Parameter &param)
{
  if (theSeries == nullptr) {
    opserr << "set/update/activate parameter is illegaly called in LoadPattern "
           << "\n";
    return 0;
  }

  if (argc < 1)
    return -1;

  // Nodal load
  if (strstr(argv[0], "loadAtNode") != 0) {

    if (argc < 3)
      return -1;

    RVisRandomProcessDiscretizer = false;

    int nodeNumber = atoi(argv[1]);
    NodalLoad *thePossibleNodalLoad;
    NodalLoad *theNodalLoad     = 0;
    NodalLoadIter &theNodalIter = this->getNodalLoads();

    while ((thePossibleNodalLoad = theNodalIter()) != 0) {
      if (nodeNumber == thePossibleNodalLoad->getNodeTag()) {
        theNodalLoad = thePossibleNodalLoad;
      }
    }

    if (theNodalLoad != 0)
      return theNodalLoad->setParameter(&argv[2], argc - 2, param);
    else
      return -1;
  }

  else if (strstr(argv[0], "elementPointLoad") != 0 ||
           strstr(argv[0], "elementLoad") != 0) {

    if (argc < 3)
      return -1;

    RVisRandomProcessDiscretizer = false;

    int eleNumber                     = atoi(argv[1]);
    ElementalLoad *theEleLoad         = 0;
    ElementalLoadIter &theEleLoadIter = this->getElementalLoads();
    while ((theEleLoad = theEleLoadIter()) != 0) {
      int eleTag = theEleLoad->getElementTag();
      if (eleNumber == eleTag) {
        return theEleLoad->setParameter(&argv[2], argc - 2, param);
      }
    }

    return -1;
  }

  else if (strstr(argv[0], "randomProcessDiscretizer") != 0) {

    if (argc < 2)
      return -1;

    RVisRandomProcessDiscretizer = true;
    return theSeries->setParameter(&argv[1], argc - 1, param);
  }

  // Unknown parameter
  else
    return -1;
}

int
LoadPattern::updateParameter(int parameterID, Information &info)
{
  if (theSeries == 0) {
    opserr << "set/update/activate parameter is illegaly called in LoadPattern "
           << "\n";
  }

  opserr << "LoadPattern::updateParameter -- no parameters defined, this "
            "method should not be called"
         << "\n";

  return 0;

  /*
  if (RVisRandomProcessDiscretizer) {
    return theSeries->updateParameter(parameterID,info);
  }
  else {
    NodalLoad *thePossibleNodalLoad = 0;
    NodalLoad *theNodalLoad = 0;
    NodalLoadIter &theNodalIter = this->getNodalLoads();
    
    switch (parameterID) {
    case 1: case -1:  // Not implemented.
      return -1;
    default:
      if (parameterID > 1000  &&  parameterID < 2000)  {
	int nodeNumber = parameterID-1000;
	while ((thePossibleNodalLoad = theNodalIter()) != 0)  {
	  if ( nodeNumber == thePossibleNodalLoad->getNodeTag() )  {
	    theNodalLoad = thePossibleNodalLoad;
	  }
	}
	return theNodalLoad->updateParameter(1, info);
      }
      else if (parameterID > 2000  &&  parameterID < 3000)  {
	int nodeNumber = parameterID-2000;
	while ((thePossibleNodalLoad = theNodalIter()) != 0)  {
	  if ( nodeNumber == thePossibleNodalLoad->getNodeTag() )  {
	    theNodalLoad = thePossibleNodalLoad;
	  }
	}
	return theNodalLoad->updateParameter(2, info);
      }
      else if (parameterID > 3000  &&  parameterID < 4000)  {
	int nodeNumber = parameterID-3000;
	while ((thePossibleNodalLoad = theNodalIter()) != 0)  {
	  if ( nodeNumber == thePossibleNodalLoad->getNodeTag() )  {
	    theNodalLoad = thePossibleNodalLoad;
	  }
	}
	return theNodalLoad->updateParameter(3, info);
            }
      else
	return -1;
    }
  }
  */
}

int LoadPattern::activateParameter(int parameterID)
{
  if (theSeries == nullptr) {
    opserr << "set/update/activate parameter is illegaly called in LoadPattern "
           << "\n";
  }

  opserr << "LoadPattern::activateParameter -- no parameters defined, this "
            "method should not be called"
         << "\n";

  return 0;

  /*
  if (RVisRandomProcessDiscretizer) {
    return theSeries->activateParameter(parameterID);
  }
  else {
    
    
    // Don't set flag here in the load pattern itself.
    // (Assume there always may be random loads)
    
    NodalLoad *theNodalLoad = 0;
    NodalLoadIter &theNodalIter = this->getNodalLoads();
    
    if (parameterID == 0) {
      
      // Go through all nodal loads and zero out gradientIdentifier
      // (Remember: the identifier is only zero if we are in
      // the process of zeroing out all sensitivity flags).
      while ((theNodalLoad = theNodalIter()) != 0)  {
	theNodalLoad->activateParameter(parameterID);
      }
      
    }
    else {
      
      // Find the right nodal load and set the flag
      if (parameterID > 1000  &&  parameterID < 2000)  {
	int nodeNumber = parameterID-1000;
	while ((theNodalLoad = theNodalIter()) != 0)  {
	  if ( nodeNumber == theNodalLoad->getNodeTag() )  {
	    theNodalLoad->activateParameter(1);
	  }
	}
      }
      else if (parameterID > 2000  &&  parameterID < 3000)  {
	int nodeNumber = parameterID-2000;
	while ((theNodalLoad = theNodalIter()) != 0)  {
	  if ( nodeNumber == theNodalLoad->getNodeTag() )  {
	    theNodalLoad->activateParameter(2);
	  }
	}
      }
      else if (parameterID > 3000  &&  parameterID < 4000)  {
	int nodeNumber = parameterID-3000;
	while ((theNodalLoad = theNodalIter()) != 0)  {
	  if ( nodeNumber == theNodalLoad->getNodeTag() )  {
	    theNodalLoad->activateParameter(3);
	  }
	}
      }
      else {
        opserr << "LoadPattern::gradient() -- error in identifier. " << "\n";
      }
    }
  }
  return 0;
  */
}

const Vector &
LoadPattern::getExternalForceSensitivity(int gradNumber)
{

  // THIS METHOD IS CURRENTLY ONLY USED FOR THE STATIC CASE
  // IT SHOULD BE DELETED AND REPLACED BY THE DYNAMIC CASE

  // Initial declarations
  Vector tempRandomLoads(1);
  int sizeRandomLoads;

  // Start with a fresh return vector
  if (randomLoads == 0) {
    randomLoads = new Vector(1);
  } else {
    delete randomLoads;
    randomLoads = new Vector(1);
  }

  // Prepare the vector identifying which loads are random.
  NodalLoad *theNodalLoad     = nullptr;
  NodalLoadIter &theNodalIter = this->getNodalLoads();

  // Loop through the nodal loads to pick up possible contributions
  int nodeNumber;
  int dofNumber;
  while ((theNodalLoad = theNodalIter()) != 0) {
    const Vector &gradientVector =
        theNodalLoad->getExternalForceSensitivity(gradNumber);
    if (gradientVector(0) != 0.0) {

      // Found a random load! Get nodeNumber and dofNumber
      nodeNumber = theNodalLoad->getNodeTag();
      dofNumber  = (int)gradientVector(0);

      // Update the randomLoads vector
      sizeRandomLoads = randomLoads->Size();
      if (sizeRandomLoads == 1) {
        delete randomLoads;
        randomLoads       = new Vector(2);
        (*randomLoads)(0) = (double)nodeNumber;
        (*randomLoads)(1) = (double)dofNumber;
      } else {
        tempRandomLoads = (*randomLoads);
        delete randomLoads;
        randomLoads = new Vector(sizeRandomLoads + 2);
        for (int i = 0; i < sizeRandomLoads; i++) {
          (*randomLoads)(i) = tempRandomLoads(i);
        }
        (*randomLoads)(sizeRandomLoads)     = nodeNumber;
        (*randomLoads)(sizeRandomLoads + 1) = dofNumber;
      }
    }
  }

  return (*randomLoads);
}

int LoadPattern::saveLoadFactorSensitivity(double dlambdadh, int gradIndex,
                                           int numGrads)
{

  if (dLambdadh == 0) {
    dLambdadh = new Vector(numGrads);
  }

  if (dLambdadh == 0 || dLambdadh->Size() != numGrads) {
    if (dLambdadh != 0)
      delete dLambdadh;
    dLambdadh = new Vector(numGrads);
  }

  if (gradIndex >= 0 && gradIndex < numGrads) {
    (*dLambdadh)(gradIndex) = dlambdadh;
    return 0;
  } else {
    opserr
        << "LoadPattern::saveLoadFactorSensitivity -- gradIndex out of bounds"
        << "\n";
    return -1;
  }
}

double LoadPattern::getLoadFactorSensitivity(int gradIndex)
{
  if (dLambdadh != 0 && gradIndex >= 0 && gradIndex < dLambdadh->Size())
    return (*dLambdadh)(gradIndex);
  else
    return 0.0;
}

// AddingSensitivity:END //////////////////////////////////////

#endif

