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
// Purpose: This file contains the implementation of MovableObject.
//
// Written: fmk
// Created: 11/96
//
#include <MovableObject.h>
#include <OPS_Globals.h>

MovableObject::MovableObject(ClassTag cTag, int dTag)
:classTag(cTag), dbTag(dTag)
{

}


MovableObject::MovableObject(ClassTag theTag)
:classTag(theTag), dbTag(0)
{

}

MovableObject::~MovableObject()
{

}




ClassTag 
MovableObject::getClassTag() const
{
  return classTag;
}

static char unknownClassType[] = {"UnknownMovableObject"};

const char *
MovableObject::getClassType() const
{
  return unknownClassType;
}


int 
MovableObject::getDbTag() const
{
    return dbTag;
}

void
MovableObject::setDbTag(int newTag)
{
  dbTag = newTag;
}

int
MovableObject::setParameter(const char **argv, int argc, Parameter &param)
{
  return -1;
}

int
MovableObject::updateParameter(int parameterID, Information &info)
{
  return 0;
}

int
MovableObject::activateParameter(int parameterID)
{
  return 0;
}

int 
MovableObject::setVariable(const char *variable, Information &theInfo)
{
  return -1;
}

int 
MovableObject::getVariable(const char *variable, Information &theInfo)
{
  return -1;
}
