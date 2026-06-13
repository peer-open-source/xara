# pragma once
#include "LoadPattern.h"
#include <TaggedIterator.hpp>
#include <MapOfTaggedObjects.h>
#include <SP_ConstraintIter.h>
#include <ranges>
#include <vector>

// case <tag> {
//   pattern <type> <tag> ...
//.  
//.  analyze 
// }


class Domain;
class LoadPattern;

class ElementalLoad;
class NodalLoad;
class SP_Constraint;
class SP_ConstraintIter;

namespace OpenSees {

class LoadCase {
public:
    LoadCase(Domain& domain);

    void clearAll();

    // int domainChange(); // TODO

    void applyLoad(double pseudoTime);
    void setLoadConstant();
    void unsetLoadConstant();

    virtual int getNumLoadPatterns(void) const;
    bool addLoadPattern(LoadPattern *);
    LoadPattern* getLoadPattern(int tag);
    LoadPattern* removeLoadPattern(int tag);
    bool addSP_Constraint(SP_Constraint *, int loadPatternTag);

    using PatternStorage = MapOfTaggedObjects;
    using PatternIterator = TaggedIterator<LoadPattern, PatternStorage>;
    PatternIterator   &getLoadPatterns();

    auto getConstraints() {
        allSPs.clear();
        PatternIterator &thePatterns = this->getLoadPatterns();
        LoadPattern *lp;
        while ((lp = thePatterns()) != nullptr) {
            SP_ConstraintIter &spIter = lp->getSPs();
            SP_Constraint *sp;
            while ((sp = spIter()) != nullptr) {
                allSPs.push_back(sp);
            }
        }
        return std::views::all(allSPs);
    }

private:
    Domain& domain;
    PatternStorage  *theLoadPatterns;
    PatternIterator *theLoadPatternIter;
    std::vector<SP_Constraint*> allSPs;
};

} // namespace OpenSees