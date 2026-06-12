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
// Description: This file contains the methods for class ElementalLoad.
//
// Written: fmk 11/95
//          modified 11/01 for new design
//

#include <ElementalLoad.h>
#include <Logging.h>
#include <Element.h>
#include <Domain.h>

ElementalLoad::ElementalLoad(int tag, int cTag, int theEleTag)
  : TaggedObject(tag)
  , MovableObject(cTag)
  , eleTag(theEleTag)
  , theElement(nullptr)
  , theDomain(nullptr)
  , loadPatternTag(-1)
{

}

ElementalLoad::ElementalLoad(int tag, int cTag)
  : 
  TaggedObject(tag)
  , MovableObject(cTag)
  , eleTag(0)
  , theElement(nullptr)
{

}


// provided for the FEM_Object broker; the tag and elementTag need
// to be supplied in recvSelf();
ElementalLoad::ElementalLoad(int cTag)
 : TaggedObject(0)
 , MovableObject(cTag)
 , eleTag(0)
 , theElement(nullptr)
{

}


ElementalLoad::~ElementalLoad()
{

}


int
ElementalLoad::setDomain(Domain *domain)
{
  this->theDomain = domain;

  if (theDomain == nullptr) {
    theElement = nullptr;
    return 0;
  }

  theElement = theDomain->getElement(eleTag);
  if (theElement == nullptr) {
    opserr << "WARNING - ElementalLoad::setDomain - no ele with tag ";
    opserr << eleTag << " exists in the domain\n";
    return -1;
  }
  return 0;
}

void 
ElementalLoad::applyLoad(double loadFactor) 
{
  if (theElement != nullptr)
    theElement->addLoad(this, loadFactor);
}

void 
ElementalLoad::applyLoad(const Vector &loadFactors) 
{
  if (theElement != nullptr)
    theElement->addLoad(this, loadFactors);
}

const Vector&
ElementalLoad::getSensitivityData(int gradIndex)
{
  static Vector dummy(10);
  return dummy;
}

int
ElementalLoad::getElementTag() const 
{
  return eleTag;
}

