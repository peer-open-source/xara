//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// File: ~/OOP/domain/loadcase/ElementalLoadIter.C
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
// Description: This file contains the method definitions for class 
// ElementalLoadIter. ElementalLoadIter is a class for iterating through the 
// elemental loads of a load case.

#include "ElementalLoadIter.h"

#include <ElementalLoad.h>
#include <TaggedObjectIter.h>
#include <TaggedObjectStorage.h>


ElementalLoadIter::ElementalLoadIter(TaggedObjectStorage *theStorage)
  :myIter(theStorage->getComponents())
{
}


ElementalLoadIter::~ElementalLoadIter()
{
}    

void
ElementalLoadIter::reset()
{
    myIter.reset();
}    


ElementalLoad *
ElementalLoadIter::operator()()
{
    // check if we still have elements in the model
    // if not return 0, indicating we are done
    TaggedObject *theComponent = myIter();
    if (theComponent == 0)
        return 0;
    else {
        ElementalLoad *result = (ElementalLoad *)theComponent;
        return result;
    }
}

