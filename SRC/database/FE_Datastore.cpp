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
                                                                        
// $Revision: 1.8 $
// $Date: 2005-11-07 23:53:45 $
// $Source: /usr/local/cvs/OpenSees/SRC/database/FE_Datastore.cpp,v $
                                                                        
// Written: fmk 
// Created: 10/98
// Revision: A
//
// Description: This file contains the class implementation for FE_Datastore.
// FE_Datastore is an abstract base class. An FE_Datastore object is used
// in the program to store/restore the geometry and state information in
// a domain at a particular instance in the analysis.
//
// What: "@(#) FE_Datastore.C, revA"

#include "FE_Datastore.h"
#include <FEM_ObjectBroker.h>
#include <MovableObject.h>
#include <Domain.h>
#include <OPS_Globals.h>
#include <ID.h>

int FE_Datastore::lastDbTag(0);

// FE_Datastore(int tag, int noExtNodes);
// 	constructor that takes the FE_Datastore's unique tag and the number
//	of external nodes for the FE_Datastore.

FE_Datastore::FE_Datastore(Domain &thDomain, FEM_ObjectBroker &theBroker) 
  :theObjectBroker(&theBroker), theDomain(&thDomain)
{

}


FE_Datastore::~FE_Datastore() 
{
    // does nothing
}

int
FE_Datastore::isDatastore()
{
  return 1;
}

/********************************************************************
 *                   CHANNEL METHODS  THAT DO NOTHING               *
 ********************************************************************/

char *
FE_Datastore::addToProgram()
{
  return 0;
}

int 
FE_Datastore::setUpConnection()
{
  return 0;
}

int 
FE_Datastore::setNextAddress(const ChannelAddress &otherChannelAddress)
{
  return 0;
}


ChannelAddress *
FE_Datastore::getLastSendersAddress(void)
{
  return 0;
}


int 
FE_Datastore::sendObj(int commitTag,
		      MovableObject &theObject, 
		      ChannelAddress *theAddress)
{
  return theObject.sendSelf(commitTag, *this);
}

int 
FE_Datastore::recvObj(int commitTag,
		      MovableObject &theObject, 
		      FEM_ObjectBroker &theNewBroker,
		      ChannelAddress *theAddress)
{
  return theObject.recvSelf(commitTag, *this, theNewBroker);
}

		


int
FE_Datastore::commitState(int commitTag)
{
  // 09/2026 - Removed functionality 
  return -1;
}



int
FE_Datastore::restoreState(int commitTag)
{
  // 09/2026 - Removed functionality 
  return -1;
}




int 
FE_Datastore::createTable(const char*table, int numColumns, char *columns[])
{
  opserr << "FE_Datastore::createTable - not yet implemented\n";
  return -1;
}


int 
FE_Datastore::insertData(const char *tableName, char *columns[], 
			int commitTag, const Vector &data)
{
  opserr << "FE_Datastore::insertData - not yet implemented\n";
  return -1;
}


int 
FE_Datastore::getData(const char *table, char *column[], int commitTag, Vector &data)
{
  opserr << "FE_Datastore::getData - not yet implemented\n";
  return -1;
}

int
FE_Datastore::getDbTag(void)
{
  lastDbTag++;
  return lastDbTag;
}
