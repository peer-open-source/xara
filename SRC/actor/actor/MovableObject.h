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
#ifndef MovableObject_h
#define MovableObject_h

// Written: fmk
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the class definition for MovableObject.
// MovableObject is meant to be an abstract base class and thus no objects 
// of it's type can be instantiated. A movable object is an object which
// can send/receive itself to/from a Channel object.
//
// What: "@(#) MovableObject.h, revA"

#include <classTags.h>
#include <cstddef>

typedef int ClassTag;

class Channel;
class FEM_ObjectBroker;
class Information;
class Parameter;

class MovableObject
{
  public:
    MovableObject(ClassTag classTag, int dbTag);        
    MovableObject(ClassTag classTag);    
    virtual ~MovableObject();

    ClassTag getClassTag() const;
    virtual const char *getClassType() const;

    int getDbTag() const;
    void setDbTag(int dbTag);

    virtual int sendSelf(int commitTag, Channel &) {
      return -1;
    }
    virtual int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &) {
      return -1;
    }
    
    // methods for sensitivity studies
    virtual int setParameter(const char **argv, int argc, Parameter &param);
    virtual int updateParameter(int parameterID, Information &info);
    virtual int activateParameter(int parameterID);

    virtual int setVariable(const char *variable, Information &);
    virtual int getVariable(const char *variable, Information &);


  private:
    ClassTag classTag;
    int dbTag;
};

#endif
