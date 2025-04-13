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
#ifndef Material_h
#define Material_h
//
// Description: This file contains the class definition for Material.
// Material is an abstract base class and thus no objects of it's type
// can be instantiated. It has pure virtual functions which must be
// implemented in it's derived classes. 
//
// Written: fmk 
// Created: 05/98
// Revision: A
//
#include <DomainComponent.h>
#include <MovableObject.h>

class OPS_Stream;
class Information;
class Response;

class Material : public TaggedObject, public MovableObject
{
  public:
    Material(int tag, int classTag);    
    virtual ~Material();


    virtual Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    virtual int getResponse(int responseID, Information &info);
    virtual int getResponseSensitivity(int responseID, int gradIndex,
				       Information &info);

    // method for this material to update itself according to its new parameters
    virtual void update() {
      return;
    }

  protected:
    
  private:
};


#endif

