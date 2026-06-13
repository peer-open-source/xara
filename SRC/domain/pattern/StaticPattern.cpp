//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
//===----------------------------------------------------------------------===//
//
#include "StaticPattern.h"
#include <ID.h>
#include <TimeSeries.h>
#include <NodalLoad.h>
#include <ElementalLoad.h>
#include <MapOfTaggedObjects.h>
#include <ElementalLoadIter.h>
#include <NodalLoadIter.h>
#include <SingleDomSP_Iter.h>
#include <SP_Constraint.h>


StaticPattern::StaticPattern(int tag, double fact, TimeSeries *theSeries)
  : LoadPattern(tag, PATTERN_TAG_StaticPattern, fact)
  , theDomain(nullptr)
  , theSeries(theSeries)
  , loadFactor(0.0)
  , scaleFactor(fact)
  , currentGeoTag(0)
  , lastGeoSendTag(-1)
  , theNodalLoads(new MapOfTaggedObjects())
  , theElementalLoads(new MapOfTaggedObjects())
  , theSPs(new MapOfTaggedObjects())
  , theNodIter(new NodalLoadIter(theNodalLoads))
  , theEleIter(new ElementalLoadIter(theElementalLoads))
  , theSpIter(new SingleDomSP_Iter(theSPs))
{
  randomLoads = 0;
  dLambdadh   = 0;
}


StaticPattern::~StaticPattern()
{
  delete theNodalLoads;
  delete theElementalLoads;
  delete theSPs;
  delete theNodIter;
  delete theEleIter;
  delete theSpIter;
  delete theSeries;
  if (randomLoads != nullptr)
    delete randomLoads;
  if (dLambdadh != nullptr)
    delete dLambdadh;
}


double
StaticPattern::getLoadFactor()
{
  if (theSeries != nullptr)
    return loadFactor;
  else
    return 0.0;
}

void
StaticPattern::applyLoad(double pseudoTime)
{
  // first determine the load factor
  if (theSeries != nullptr && isConstant != true) {
    loadFactor = theSeries->getFactor(pseudoTime);
    loadFactor *= scaleFactor;
  }

  {
    NodalLoad *nodLoad;
    NodalLoadIter &theNodalIter = this->getNodalLoads();
    while ((nodLoad = theNodalIter()) != nullptr)
      nodLoad->applyLoad(loadFactor);
  }

  {
    ElementalLoad *eleLoad;
    ElementalLoadIter &theElementalIter = this->getElementalLoads();
    while ((eleLoad = theElementalIter()) != nullptr)
      eleLoad->applyLoad(loadFactor);
  }

  // call at end so that loadFactor is set
  this->LoadPattern::applyLoad(pseudoTime);
  // SP_Constraint *sp;
  // SP_ConstraintIter &theIter = this->getSPs();
  // while ((sp = theIter()) != nullptr)
  //   sp->applyConstraint(loadFactor);
}



bool
StaticPattern::addNodalLoad(NodalLoad *load)
{
  Domain *theDomain = this->getDomain();

  bool result = theNodalLoads->addComponent(load);

  if (result != true) {
    opserr << "WARNING: Load could not be added\n";
    return result;
  }

  if (theDomain != nullptr)
    load->setDomain(theDomain);

  load->setLoadPatternTag(this->getTag());
  currentGeoTag++;

  return result;
}


bool
StaticPattern::addElementalLoad(Domain& domain, ElementalLoad *load)
{
  Domain* theDomain = &domain; // TODO: just use domain.(...) below
  // Domain *theDomain = this->getDomain();

  bool result = theElementalLoads->addComponent(load);
  if (result != true) {
    opserr << "WARNING: LoadPattern::addElementalLoad() - load could not be "
              "added\n";
    return false;
  }
  load->setDomain(theDomain);
  load->setLoadPatternTag(this->getTag());
  currentGeoTag++;
  return result;
}


NodalLoad *
StaticPattern::removeNodalLoad(int tag)
{
  TaggedObject *obj = theNodalLoads->removeComponent(tag);
  if (obj == nullptr)
    return nullptr;

  NodalLoad *result = (NodalLoad *)obj;
  result->setDomain(nullptr);
  currentGeoTag++;
  return result;
}


ElementalLoad *
StaticPattern::removeElementalLoad(int tag)
{
  TaggedObject *obj = theElementalLoads->removeComponent(tag);
  if (obj == nullptr)
    return nullptr;

  ElementalLoad *result = (ElementalLoad *)obj;
  result->setDomain(nullptr);
  currentGeoTag++;
  return result;
}


NodalLoadIter &
StaticPattern::getNodalLoads()
{
  theNodIter->reset();
  return *theNodIter;
}

ElementalLoadIter &
StaticPattern::getElementalLoads()
{
  theEleIter->reset();
  return *theEleIter;
}


void
StaticPattern::clearAll()
{
  theElementalLoads->clearAll();
  theNodalLoads->clearAll();
  theSPs->clearAll();
  currentGeoTag++;
  if (dLambdadh != nullptr) {
    dLambdadh->Zero();
  }

  this->LoadPattern::clearAll();
}

void
StaticPattern::setDomain(Domain *theDomain)
{
  this->theDomain = theDomain;

  // if subclass does not implement .. check for 0 pointer
  if (theNodalLoads != nullptr) {
    NodalLoad *nodLoad;
    NodalLoadIter &theNodalIter = this->getNodalLoads();
    while ((nodLoad = theNodalIter()) != nullptr)
      nodLoad->setDomain(theDomain);

    ElementalLoad *eleLoad;
    ElementalLoadIter &theElementalIter = this->getElementalLoads();
    while ((eleLoad = theElementalIter()) != nullptr)
      eleLoad->setDomain(theDomain);

    SP_Constraint *theSP;
    SP_ConstraintIter &theSpConstraints = this->getSPs();
    while ((theSP = theSpConstraints()) != nullptr)
      theSP->setDomain(theDomain);
  }


  this->LoadPattern::setDomain(theDomain);
}



void 
StaticPattern::applyLoadSensitivity(double pseudoTime)
{
  // P*dfactor/dh
  if (theSeries != 0 && isConstant != true) {
    loadFactor = theSeries->getFactorSensitivity(pseudoTime);
    loadFactor *= scaleFactor;
  }

  NodalLoad *nodLoad;
  NodalLoadIter &theNodalIter = this->getNodalLoads();
  while ((nodLoad = theNodalIter()) != nullptr)
    nodLoad->applyLoad(loadFactor);

  // factor*dP/dh
  if (theSeries != nullptr && isConstant != true) {
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
StaticPattern::setParameter(const char **argv, int argc, Parameter &param)
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
StaticPattern::updateParameter(int parameterID, Information &info)
{
  if (theSeries == nullptr) {
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

int
StaticPattern::activateParameter(int parameterID)
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
StaticPattern::getExternalForceSensitivity(int gradNumber)
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


int
StaticPattern::saveLoadFactorSensitivity(double dlambdadh, int gradIndex,
                                           int numGrads)
{

  if (dLambdadh == nullptr) {
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

double 
StaticPattern::getLoadFactorSensitivity(int gradIndex)
{
  if (dLambdadh != 0 && gradIndex >= 0 && gradIndex < dLambdadh->Size())
    return (*dLambdadh)(gradIndex);
  else
    return 0.0;
}


void
StaticPattern::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"type\": \"Static\"" << ", ";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"scale\": " << scaleFactor << ", ";
    if (theSeries) {
      s << "\"series\": ";
      theSeries->Print(s, flag); 
      s << ", ";
    }
    s << "\"nodes\": [\n";
    theNodalLoads->Print(s, flag);
    // if ((theSPs->getNumComponents() > 0) && (theNodalLoads->getNumComponents() > 0))
    //   s << ",\n";
    // theSPs->Print(s, flag);
    s << "\n" << OPS_PRINT_JSON_MATE_INDENT <<  "],\n";
    s << OPS_PRINT_JSON_MATE_INDENT << "\"elements\": [\n";
    //theElementalLoads->Print(s, flag);
    s << "\n" << OPS_PRINT_JSON_MATE_INDENT <<  "]\n";
    s << OPS_PRINT_JSON_MATE_INDENT << "}";
    return;
  }
  
  else {
    s << "Load Pattern: " << this->getTag() << "\n";
    // s << "  Scale Factor: " << scaleFactor
      // << "\n";
    if (theSeries != 0)
      theSeries->Print(s, flag);
    s << "  Nodal Loads: \n";
    theNodalLoads->Print(s, flag);
    s << "\n  Elemental Loads: \n";
    theElementalLoads->Print(s, flag);
    s << "\n  Single Point Constraints: \n";
    // theSPs->Print(s, flag);
  }
}