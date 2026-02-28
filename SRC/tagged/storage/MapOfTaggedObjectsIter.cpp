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
// Description: This file contains the implementation of MapOfTaggedObjects.
//
// $Date: 2000-09-15 08:23:30 $
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
#include <map>
#include <MapOfTaggedObjectsIter.h>
#include <MapOfTaggedObjects.h>

// some typedefs that will be useful
typedef std::map<int, TaggedObject *> MAP_TAGGED;
typedef MAP_TAGGED::value_type   MAP_TAGGED_TYPE;
typedef MAP_TAGGED::iterator     MAP_TAGGED_ITERATOR;

//	constructor that takes the model, just the basic iter
MapOfTaggedObjectsIter::MapOfTaggedObjectsIter(MapOfTaggedObjects &theComponents)
{
    theMap = &(theComponents.theMap);
}


MapOfTaggedObjectsIter::~MapOfTaggedObjectsIter()
{

}    

void
MapOfTaggedObjectsIter::reset()
{
    currentComponent = theMap->begin();
}

TaggedObject *
MapOfTaggedObjectsIter::operator()(void)
{
    if (currentComponent != theMap->end()) {
	TaggedObject *result = (*currentComponent).second;
	currentComponent++;
	return result;
    } else
	return 0;
}

    
    
