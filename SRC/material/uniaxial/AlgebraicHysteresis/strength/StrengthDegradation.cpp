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
                                                                        
// $Revision$
// $Date$
// $Source$

// Written: MHS
// Created: August 2000
//
// Description: This file contains the interface for StrengthDegradation,
// which models hysteretic strength degradation.

#include <StrengthDegradation.h>

#include <TaggedObject.h>
#include <MapOfTaggedObjects.h>

StrengthDegradation::StrengthDegradation(int tag, int classTag)
  :MaterialState(tag,classTag)
{
  
}

StrengthDegradation::~StrengthDegradation()
{
  
}

StrengthDegradation*
StrengthDegradation::getCopy(UniaxialMaterial *u)
{
  return this->getCopy();
}
