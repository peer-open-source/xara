//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
// Written: fmk 
// Created: 07/99
// Revision: A
//
// Purpose: This file contains the class definition for LoadPattern.
// LoadPattern is an *abstract* class in Xara. 
//
#pragma once
#include <TaggedObject.h>
#include <MovableObject.h>
#include <Vector.h>

class Domain;
class NodalLoad;
class TimeSeries;
class ElementalLoad;
class SP_Constraint;
class NodalLoadIter;
class ElementalLoadIter;
class SingleDomSP_Iter;
class SP_ConstraintIter;
class TaggedObjectStorage;
class GroundMotion;
class MapOfTaggedObjects;

// for applyResidual
class AnalysisModel;
class LinearSOE;

class LoadPattern : public TaggedObject, public MovableObject
{
  public:
    LoadPattern(int tag, int classTag, double fact); // for subclasses

    virtual ~LoadPattern();

    // method to set the associated TimeSeries and Domain
    virtual void setDomain(Domain *);


    // methods to apply loads

    // apply load at start of a step.
    virtual int applyResidual(AnalysisModel&, Vector&, double) {return 0;}
    void setLoadConstant();
    void unsetLoadConstant();


    // methods to add loads
    virtual bool addSP_Constraint(SP_Constraint *);
    virtual SP_ConstraintIter &getSPs();
    virtual SP_Constraint *removeSP_Constraint(int tag);

    virtual double getLoadFactor();


protected:
    virtual Domain* getDomain() {return theDomain;}

public:
    virtual void applyLoad(double pseudoTime = 0.0)=0;
    // Sensitivity
    virtual void applyLoadSensitivity(double pseudoTime = 0.0) {}
    virtual int  setParameter(const char **argv, int argc, Parameter &) {return -1;};
    virtual int  updateParameter(int parameterID, Information &) {return -1;};
    virtual int  activateParameter(int parameterID) {return -1;};
    virtual const Vector & getExternalForceSensitivity(int gradNumber) {
      static Vector dummy(0);
      return dummy;
    }
    virtual int saveLoadFactorSensitivity(double dlambdadh, int gradIndex, int numGrads) {return 0;}
    virtual double getLoadFactorSensitivity(int gradIndex) {return 0.0;}


    virtual void clearAll();


    // TaggedObject
    virtual void Print(OPS_Stream &s, int flag)=0;

    enum : int {
      PATTERN_TAG_StaticPattern = 1000
    };

  protected:
    bool   isConstant;     // to indicate whether setConstant has been called


  private:
    double loadFactor;     // current load factor
    double scaleFactor;    // factor to scale load factor from time series

    TimeSeries *theSeries; // pointer to associated TimeSeries

    int	   currentGeoTag;
    int    lastGeoSendTag;
    int    dbSPs, dbNod, dbEle; // database tags for storing info about components
    
    // storage objects for the loads and constraints
    MapOfTaggedObjects  *theNodalLoads;
    TaggedObjectStorage  *theElementalLoads;
    TaggedObjectStorage  *theSPs;

    // iterator objects for the objects added to the storage objects
    NodalLoadIter       *theNodIter;
    ElementalLoadIter   *theEleIter;
    SingleDomSP_Iter    *theSpIter;
//
    Vector *randomLoads;
    bool RVisRandomProcessDiscretizer;
    Vector *dLambdadh;

  //
    int lastChannel; 
    Domain* theDomain;
};
