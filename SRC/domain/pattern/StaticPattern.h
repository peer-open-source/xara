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
class GroundMotion;


class StaticPattern : public LoadPattern
{
  public:
    StaticPattern(int tag, double fact = 1.0);
    ~StaticPattern() override;

#if 0
    LoadPattern *getCopy() override;
    bool addSP_Constraint(SP_Constraint *) final;
    bool addNodalLoad(NodalLoad *) final;
    bool addElementalLoad(ElementalLoad *) final;
    NodalLoadIter     &getNodalLoads() final;
    ElementalLoadIter &getElementalLoads() final;    
    SP_ConstraintIter &getSPs() final;
    // methods to remove loads
    NodalLoad *removeNodalLoad(int tag) final;
    ElementalLoad *removeElementalLoad(int tag) final;
    SP_Constraint *removeSP_Constraint(int tag) final;
    void clearAll();
#endif
private:
#if 0
    // storage objects for the loads and constraints
    TaggedObjectStorage  *theNodalLoads;
    TaggedObjectStorage  *theElementalLoads;
    TaggedObjectStorage  *theSPs; 	  

    // iterator objects for the objects added to the storage objects
    NodalLoadIter       *theNodIter;
    ElementalLoadIter   *theEleIter;
    SingleDomSP_Iter    *theSpIter;
#endif
};