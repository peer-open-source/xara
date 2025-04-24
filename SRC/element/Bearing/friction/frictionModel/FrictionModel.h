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
#ifndef FrictionModel_h
#define FrictionModel_h

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 02/06
// Revision: A
//
// Description: This file contains the class definition for FrictionModel.
// FrictionModel is an abstract base class and thus no objects of it's type
// can be instantiated. It has pure virtual functions which must be
// implemented in it's derived classes. 

#include <DomainComponent.h>
#include <MovableObject.h>
#include <TaggedObject.h>
#include <Vector.h>
#include <api/runtimeAPI.h>

class Response;

class FrictionModel : public TaggedObject, public MovableObject
{
public:
    FrictionModel(int tag, int classTag);
    virtual ~FrictionModel();
    
    // public methods to set and obtain response
    virtual int setTrial(double normalForce, double velocity = 0.0) = 0;
    virtual double getNormalForce();
    virtual double getVelocity();
    virtual double getFrictionForce() = 0;
    virtual double getFrictionCoeff() = 0;
    virtual double getDFFrcDNFrc() = 0;
    virtual double getDFFrcDVel() = 0;
    
    virtual int commitState() = 0;
    virtual int revertToLastCommit() = 0;
    virtual int revertToStart() = 0;
    
    virtual FrictionModel *getCopy() = 0;
    
    virtual Response *setResponse(const char **argv, int argc, OPS_Stream &);
    virtual int getResponse(int responseID, Information &);
    
    virtual int sendSelf(int commitTag, Channel &) = 0;
    virtual int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &) = 0;
    
    virtual void Print(OPS_Stream &s, int flag = 0) = 0;
    
protected:
    double trialN;      // trial normal contact force
    double trialVel;    // trial velocity
    
private:

};

#define OPS_getFrictionModel(tag) G3_getSafeBuilder(rt)->getTypedObject<FrictionModel>(tag)

#endif
