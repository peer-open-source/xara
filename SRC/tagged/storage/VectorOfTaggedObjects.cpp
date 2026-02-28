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
// Purpose: This file contains the implementation of the VectorOfTaggedObjects
// class.
//
// Written: cmp 
// Created: Jan 2026
// Revision: A
//
#include <vector>
#include <TaggedObject.h>
#include "VectorOfTaggedObjects.h"
#include <OPS_Stream.h>

VectorOfTaggedObjects::VectorOfTaggedObjects()
:myIter(*this)
{
  // creates the iter with this as the argument
}

VectorOfTaggedObjects::~VectorOfTaggedObjects()
{
  this->clearAll();
}


int
VectorOfTaggedObjects::setSize(int newSize)
{
  // no setSize for map template .. can only check enough space available
  int maxSize = int(container.max_size());
  if (newSize > maxSize)
    return -1;

  container.reserve(newSize);
  return 0;
}


bool 
VectorOfTaggedObjects::addComponent(TaggedObject *newComponent)
{
  int tag = newComponent->getTag();

  // check if the ele already in map, if not we add
  for (auto const& item : container) {
    if (item->getTag() == tag) {
      return false;
    }
  }

  container.push_back(newComponent);
  return true;  // o.k.
}


TaggedObject *
VectorOfTaggedObjects::removeComponent(int tag)
{
  for (auto it = container.begin(); it != container.end(); ++it) {
    if ((*it)->getTag() == tag) {
      TaggedObject* item = *it;
      container.erase(it);
      return item;
    }
  }
  return nullptr;
}


int
VectorOfTaggedObjects::getNumComponents() const
{
  return int(container.size());
}


TaggedObject *
VectorOfTaggedObjects::getComponentPtr(int tag)
{
  for (TaggedObject* item : container) {
    if (item->getTag() == tag) {
      return item;
    }
  }
  return nullptr;
}


VectorOfTaggedObjects::Iterator &
VectorOfTaggedObjects::getIterRef()
{
  myIter.reset();
  return myIter;
}

VectorOfTaggedObjects::Iterator
VectorOfTaggedObjects::getIter()
{
  return Iterator(*this);
}


TaggedObjectStorage *
VectorOfTaggedObjects::getEmptyCopy()
{
  VectorOfTaggedObjects *theCopy = new VectorOfTaggedObjects();
  return theCopy;
}


void
VectorOfTaggedObjects::clearAll(bool invokeDestructor)
{
  // invoke the destructor on all the tagged objects stored
  if (invokeDestructor == true) {
    for (TaggedObject* item : container) {
      delete item;
    }
  }

  // now clear the map of all entries
  container.clear();
}


void
VectorOfTaggedObjects::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    bool first = true;
    for (TaggedObject* p : container) {
      if (!first)
        s << ", ";
      first = false;
      p->Print(s, flag);
    }
    return;
  }

  for (TaggedObject* p : container)
    p->Print(s, flag);
}



VectorOfTaggedObjects::Iterator::Iterator(VectorOfTaggedObjects &theComponents)
{
  container = &(theComponents.container);
  this->reset();
}


VectorOfTaggedObjects::Iterator::~Iterator()
{

}    

void
VectorOfTaggedObjects::Iterator::reset()
{
  currentComponent = container->begin();
}


TaggedObject *
VectorOfTaggedObjects::Iterator::operator()()
{
  if (currentComponent != container->end()) {
      TaggedObject *result = *currentComponent;
      currentComponent++;
      return result;
  } else
      return 0;
}

