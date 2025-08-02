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
// Purpose: This file contains the implementation for the Load class.
//                                                                        
// File: ~/domain/load/Load.C
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
#include <Load.h>


Load::Load(int tag, int clasTag)
:TaggedObject(tag), MovableObject(clasTag), loadPatternTag(-1)
{
    // does nothing
}


Load::~Load()
{
    // does nothing
}

void
Load::setLoadPatternTag(int tag)
{
  loadPatternTag = tag;
}

int
Load::getLoadPatternTag(void) const
{
  return loadPatternTag;
}

void
Load::setDomain(Domain *model)
{
    // sets the pointer 
    theDomain = model;
}


Domain *
Load::getDomain(void) const
{
  // returns the current pointer
  return theDomain;
}
