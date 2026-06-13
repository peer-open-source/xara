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
//===----------------------------------------------------------------------===//
//
// Written: cmp
//          Stanford University
//
#pragma once 
#include <LoadPattern.h>

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
class MapOfTaggedObjects;


class StaticPattern : public LoadPattern
{
  public:
    StaticPattern(int tag, double fact, TimeSeries *);
    ~StaticPattern() override;

    bool addNodalLoad(NodalLoad *);
    bool addElementalLoad(Domain&, ElementalLoad *);
    NodalLoadIter     &getNodalLoads();
    ElementalLoadIter &getElementalLoads();
    NodalLoad *removeNodalLoad(int tag);
    ElementalLoad *removeElementalLoad(int tag);

    //
    // LoadPattern interface
    //
#if 0
    bool addSP_Constraint(SP_Constraint *) final;
    SP_ConstraintIter &getSPs() final;
    SP_Constraint *removeSP_Constraint(int tag) final;
#endif

    void setDomain(Domain *) override;
    Domain* getDomain() {return theDomain;}
    void applyLoad(double pseudoTime) override;
    double getLoadFactor();
    void clearAll() override;


    // Sensitivity
    void applyLoadSensitivity(double pseudoTime = 0.0) override;
    int  setParameter(const char **argv, int argc, Parameter &) override;
    int  updateParameter(int parameterID, Information &) override;
    int  activateParameter(int parameterID) override;
    const Vector & getExternalForceSensitivity(int gradNumber) override;
    int saveLoadFactorSensitivity(double dlambdadh, int gradIndex, int numGrads) override;
    double getLoadFactorSensitivity(int gradIndex) override;

    //
    // TaggedObject interface
    //
    void Print(OPS_Stream &s, int flag) override;

private:
    Domain* theDomain;

    double loadFactor;     // current load factor
    double scaleFactor;    // factor to scale load factor from time series
    TimeSeries *theSeries; // pointer to associated TimeSeries

    // storage objects for the loads and constraints
    MapOfTaggedObjects  *theNodalLoads;
    TaggedObjectStorage  *theElementalLoads;
    TaggedObjectStorage  *theSPs;

    // iterator objects for the objects added to the storage objects
    NodalLoadIter       *theNodIter;
    ElementalLoadIter   *theEleIter;
    SingleDomSP_Iter    *theSpIter;

    int	   currentGeoTag;
    int    lastGeoSendTag;


    Vector *randomLoads;
    bool RVisRandomProcessDiscretizer;
    Vector *dLambdadh;
};