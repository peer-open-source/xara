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
// Purpose: This file contains the class definition for ElementalLoad.
// ElementalLoad is an abstract class.
//
// Written: fmk 
//
#ifndef ElementalLoad_h
#define ElementalLoad_h


#include <Load.h>
#include <Vector.h>

class Element;

class ElementalLoad : public Load
{
  public:
    ElementalLoad(int tag, int classTag, int eleTag);
    ElementalLoad(int tag, int classTag);
    ElementalLoad(int classTag);
    ~ElementalLoad();

    virtual void setDomain(Domain *theDomain);
    virtual void applyLoad(double loadfactor);
    virtual void applyLoad(const Vector &loadfactors);
    virtual const Vector &getData(int &type, double loadFactor) = 0;
    virtual const Vector &getSensitivityData(int gradIndex);

    virtual int getElementTag();

  protected:
    int eleTag;
    Element *theElement;

  private:

};

#endif

