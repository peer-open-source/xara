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
// Purpose: This file contains the implementation of the MapOfTaggedObjects
// class.
//
// Written: fmk 
// Created: 02/00
// Revision: A
//
#include <map>
#include <TaggedObject.h>
#include <MapOfTaggedObjects.h>
#include <OPS_Globals.h>

// some typedefs that will be useful
typedef Map<int, TaggedObject*> MAP_TAGGED;
typedef MAP_TAGGED::value_type        MAP_TAGGED_TYPE;
typedef MAP_TAGGED::iterator          MAP_TAGGED_ITERATOR;

MapOfTaggedObjects::MapOfTaggedObjects()
:myIter(*this)
{
    // creates the iter with this as the argument
}

MapOfTaggedObjects::~MapOfTaggedObjects()
{
    this->clearAll();
}


int
MapOfTaggedObjects::setSize(int newSize)
{
    // no setSize for map template .. can only check enough space available
    int maxSize = int(theMap.max_size());
    if (newSize > maxSize) {
      opserr << "MapOfTaggedObjects::setSize - failed as map STL has a max size of " << maxSize << "\n";
      return -1;
    } 
   
    return 0;
}


bool 
MapOfTaggedObjects::addComponent(TaggedObject *newComponent)
{
    MAP_TAGGED_ITERATOR theEle;
    int tag = newComponent->getTag();

    // check if the ele already in map, if not we add
    std::pair<MAP_TAGGED_ITERATOR,bool> res = theMap.insert(MAP_TAGGED_TYPE(tag,newComponent));    
    if (res.second == false) {
      opserr << "MapOfTaggedObjects::addComponent - not adding as one with similar tag exists, tag: " <<
	tag << "\n";
      return false;
    }

    /*
    theEle = theMap.find(tag);
    if (theEle == theMap.end()) {
      theMap.insert(MAP_TAGGED_TYPE(tag,newComponent));
      
      // check if sucessfully added 
      theEle = theMap.find(tag);
      if (theEle == theMap.end()) {
	opserr << "MapOfTaggedObjects::addComponent - map STL failed to add object with tag : " << 
	  newComponent->getTag() << "\n";
	return false;
      }
    }  else {	
      opserr << "MapOfTaggedObjects::addComponent - not adding as one with similar tag exists, tag: " <<
	tag << "\n";
      return false;
    }
    */

    return true;  // o.k.
}


TaggedObject *
MapOfTaggedObjects::removeComponent(int tag)
{
    TaggedObject *removed =0;
    MAP_TAGGED_ITERATOR theEle;

    // return 0 if component does not exist, otherwise remove it
    theEle = theMap.find(tag);
    if (theEle == theMap.end()) // the object has not been added
	return nullptr;
    else { // the object exists so we remove it
	removed = (*theEle).second;
	int ok = int(theMap.erase(tag));
	if (ok != 1) { // ensure the map did remove the object
	  opserr << "MapOfTaggedObjects::removeComponent - map STL failed to remove object with tag " << 
	    tag << "\n";
	  return nullptr;
	}
    }

    return removed;
}


int
MapOfTaggedObjects::getNumComponents(void) const
{
    return int(theMap.size());
}


TaggedObject *
MapOfTaggedObjects::getComponentPtr(int tag)
{
    TaggedObject *removed = nullptr;
    MAP_TAGGED_ITERATOR theEle;
    
    // return 0 if component does not exist, otherwise remove it
    theEle = theMap.find(tag);
    if (theEle == theMap.end()) 
	return nullptr;
    else 
	removed = (*theEle).second;
    
    return removed;
}


TaggedObjectIter &
MapOfTaggedObjects::getComponents()
{
    myIter.reset();
    return myIter;
}


MapOfTaggedObjectsIter 
MapOfTaggedObjects::getIter()
{
    return MapOfTaggedObjectsIter(*this);
}


TaggedObjectStorage *
MapOfTaggedObjects::getEmptyCopy(void)
{
    MapOfTaggedObjects *theCopy = new MapOfTaggedObjects();
    
    if (theCopy == nullptr) {
      opserr << "MapOfTaggedObjects::getEmptyCopy-out of memory\n";
    }

    return theCopy;
}

void
MapOfTaggedObjects::clearAll(bool invokeDestructor)
{

    // invoke the destructor on all the tagged objects stored
    if (invokeDestructor == true) {
	MAP_TAGGED_ITERATOR p = theMap.begin();
	while (p != theMap.end()) {
	    delete (*p).second;
	    p++;
	}
    }

    // now clear the map of all entries
    theMap.clear();
}

void
MapOfTaggedObjects::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
      int i = 0;
      MAP_TAGGED_ITERATOR p = theMap.begin();
      while (p != theMap.end()) {
        if (i++)
          s << ",\n";
        ((*p).second)->Print(s, flag);
        p++;
      }
      return;
    }

    // s << "\nnumComponents: " << this->getNumComponents();
    // go through the array invoking Print on non-zero entries
    MAP_TAGGED_ITERATOR p = theMap.begin();
    while (p != theMap.end()) {
	((*p).second)->Print(s, flag);
	p++;
    }
}


