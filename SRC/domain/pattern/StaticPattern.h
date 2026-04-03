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
class GroundMotion;


class StaticPattern : public LoadPattern
{
  public:
    StaticPattern(int tag, double fact = 1.0);
    ~StaticPattern() override;

#if 0
    LoadPattern *getCopy() override;
    bool addSP_Constraint(SP_Constraint *) final;
    SP_ConstraintIter &getSPs() final;
    SP_Constraint *removeSP_Constraint(int tag) final;
private:
    TaggedObjectStorage  *theSPs; 	  

#endif
#if 0
    bool addNodalLoad(NodalLoad *) final;
    bool addElementalLoad(ElementalLoad *) final;
    NodalLoadIter     &getNodalLoads() final;
    ElementalLoadIter &getElementalLoads() final;

    // methods to remove loads
    NodalLoad *removeNodalLoad(int tag) final;
    ElementalLoad *removeElementalLoad(int tag) final;
    void clearAll();
#endif
private:
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
};