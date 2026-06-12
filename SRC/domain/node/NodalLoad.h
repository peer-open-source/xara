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
#ifndef NodalLoad_h
#define NodalLoad_h

// Written: fmk 
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the class interface for NodalLoad.
// NodalLoad is a class for applying nodal loads to the model.

#include <MovableObject.h>
#include <TaggedObject.h>
#include <Node.h>
#include <Vector.h>
class Domain;

class NodalLoad : public TaggedObject, public MovableObject
{
  public:
    NodalLoad(int classTag);
    NodalLoad(int tag, int node, int classTag);
    NodalLoad(int tag, int node, const Vector &load, bool isLoadConstant = false);
    virtual ~NodalLoad();

    int setDomain(Domain *);
    Domain *getDomain() const {return theDomain;}

    int getNodeTag() const;
    virtual void applyLoad(double loadFactor);
    virtual void applyLoadSensitivity(double loadFactor);


    void setLoadPatternTag(int tag) {loadPatternTag = tag;}
    int  getLoadPatternTag() const {return loadPatternTag;}
    
    // Sensitivity
    int setParameter(const char **argv, int argc, Parameter &);
    int            updateParameter(int parameterID, Information &);
    int            activateParameter(int parameterID);
    const Vector & getExternalForceSensitivity(int gradNumber);

    //Change made by Liming for NodalThermalAction [SIF]
    virtual void applyLoad(Vector& loadFactors);
    virtual const Vector &getData(int& type);
    // Change made by Liming for NodalThermalAction [SIF]


    int sendSelf(int commitTag, Channel &) override;
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &) override;
    
    void Print(OPS_Stream &s, int flag) override;

  private:
    int  myNode;        // tag indicating associated Node objects tag
    Node *myNodePtr;    // pointer to Node object on which load acts
    Vector *load;       // the reference load - pointer to new copy or 0
    bool  konstant;     // true if load is load factor independent
    // Sensitivity
    int parameterID;
    static Vector gradientVector;
    Domain* theDomain;

    int loadPatternTag;
};

#endif

