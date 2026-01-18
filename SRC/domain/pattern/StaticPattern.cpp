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

StaticPattern::StaticPattern(int tag, double fact)
  : LoadPattern(tag, PATTERN_TAG_StaticPattern, fact)
  , currentGeoTag(0)
  , lastGeoSendTag(-1)
  , theNodalLoads(new MapOfTaggedObjects())
  , theElementalLoads(new MapOfTaggedObjects())
  , theSPs(new MapOfTaggedObjects())
  , theNodIter(new NodalLoadIter(theNodalLoads))
  , theEleIter(new ElementalLoadIter(theElementalLoads))
  , theSpIter(new SingleDomSP_Iter(theSPs))
{

}

StaticPattern::~StaticPattern()
{
    delete theNodalLoads;
    delete theElementalLoads;
    delete theSPs;
    delete theNodIter;
    delete theEleIter;
    delete theSpIter;
}

#if 0
bool
StaticPattern::addNodalLoad(NodalLoad *load)
{
  Domain *theDomain = this->getDomain();

  bool result = theNodalLoads->addComponent(load);

  if (result != true) {
    opserr << "WARNING: LoadPattern::addNodalLoad() - load could not be added\n";
    return result;
  }

  if (theDomain != nullptr)
    load->setDomain(theDomain);

  load->setLoadPatternTag(this->getTag());
  currentGeoTag++;

  return result;
}


bool
StaticPattern::addElementalLoad(ElementalLoad *load)
{
  Domain *theDomain = this->getDomain();

  bool result = theElementalLoads->addComponent(load);
  if (result == true) {
    if (theDomain != 0)
      load->setDomain(theDomain);
    load->setLoadPatternTag(this->getTag());
    currentGeoTag++;
  } else
    opserr << "WARNING: LoadPattern::addElementalLoad() - load could not be "
              "added\n";

  return result;
}


NodalLoad *
StaticPattern::removeNodalLoad(int tag)
{
  TaggedObject *obj = theNodalLoads->removeComponent(tag);
  if (obj == 0)
    return 0;
  NodalLoad *result = (NodalLoad *)obj;
  result->setDomain(nullptr);
  currentGeoTag++;
  return result;
}


ElementalLoad *
StaticPattern::removeElementalLoad(int tag)
{
  TaggedObject *obj = theElementalLoads->removeComponent(tag);
  if (obj == 0)
    return 0;

  ElementalLoad *result = (ElementalLoad *)obj;
  result->setDomain(nullptr);
  currentGeoTag++;
  return result;
}


NodalLoadIter &
LoadPattern::getNodalLoads()
{
  theNodIter->reset();
  return *theNodIter;
}

ElementalLoadIter &
LoadPattern::getElementalLoads()
{
  theEleIter->reset();
  return *theEleIter;
}

#endif