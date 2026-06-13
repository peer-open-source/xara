

#include "LoadCase.h"
#include <MapOfTaggedObjects.h>
#include <LoadPattern.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <MP_Constraint.h>
#include <MP_ConstraintIter.h>
#include <NodalLoad.h>
#include <NodalLoadIter.h>
#include <ElementalLoad.h>
#include <ElementalLoadIter.h>

#include <Element.h>
#include <ElementIter.h>
#include <Node.h>
#include <NodeIter.h>
#include <Domain.h>
#include <Matrix.h>

using namespace OpenSees;

LoadCase::LoadCase(Domain& domain_)
: domain(domain_)
, theLoadPatterns(new PatternStorage())
, theLoadPatternIter(new PatternIterator(theLoadPatterns))
{

}

void 
LoadCase::clearAll()
{
  // clear the loads and constraints from any load pattern
  PatternIterator &thePatterns = this->getLoadPatterns();
  LoadPattern *thePattern;
  while ((thePattern = thePatterns()) != nullptr)
    thePattern->clearAll();
}

bool 
LoadCase::addLoadPattern(LoadPattern *load)
{
  // first check if a load pattern with a similar tag exists
  int tag = load->getTag();
  TaggedObject *other = theLoadPatterns->getComponentPtr(tag);
  if (other != nullptr) {
    opserr << "LoadCase::addLoadPattern - cannot add as LoadPattern with tag "
            << tag << " already exists in model\n";
    return false;
  }    

  int numSPs = 0;
  SP_Constraint *theSP_Constraint;
  SP_ConstraintIter &theSP_Constraints = load->getSPs();
  while ((theSP_Constraint = theSP_Constraints()) != nullptr)
    numSPs++;
  
  // now we add the load pattern to the container for load patterns
  bool result = theLoadPatterns->addComponent(load);
  if (result == true) {
    load->setDomain(&domain);
    if (numSPs > 0)
      domain.domainChange();
  }
  else 
    opserr << "cannot add LoadPattern with tag " 
            << tag << " to the container\n";              	
  return result;
}    


bool
LoadCase::addSP_Constraint(SP_Constraint *spConstraint, int pattern)
{

  // now add it to the pattern
  TaggedObject *thePattern = theLoadPatterns->getComponentPtr(pattern);
  if (thePattern == nullptr) {
      opserr << "LoadCase::addSP_Constraint - cannot add as pattern with tag "
             << pattern << " does not exist in domain\n"; 

      return false;
  }

  LoadPattern *theLoadPattern = (LoadPattern *)thePattern;
  bool result = theLoadPattern->addSP_Constraint(spConstraint);
  if (result == false) {
    opserr << "LoadCase::addSP_Constraint - " 
           << pattern << " pattern could not add the SP_Constraint\n"; 
			      
    return false;
  }

  spConstraint->setDomain(&domain);
  domain.domainChange();  

  return true;
}


LoadPattern *
LoadCase::removeLoadPattern(int tag)
{
  // remove the object from the container            
  TaggedObject *obj = theLoadPatterns->removeComponent(tag);
  
  // if not there return nullptr    
  if (obj == nullptr)
    return nullptr;
  
  // perform a downward cast, set the objects domain pointer to 0
  // and return the result of the cast            
  LoadPattern *result = (LoadPattern *)obj;
  // result->setDomain(0);

  //
  // now set the Domain pointer for all loads and SP constraints 
  // in the loadPattern to be 0
  //
#if 0
  NodalLoad *theNodalLoad;
  NodalLoadIter &theNodalLoads = result->getNodalLoads();
  while ((theNodalLoad = theNodalLoads()) != nullptr) {
    // theNodalLoad->setDomain(0);
  }

  ElementalLoad *theElementalLoad;
  ElementalLoadIter &theElementalLoads = result->getElementalLoads();
  while ((theElementalLoad = theElementalLoads()) != nullptr) {
    // theElementalLoad->setDomain(0);
  }
#endif

  int numSPs = 0;
  SP_Constraint *theSP_Constraint;
  SP_ConstraintIter &theSP_Constraints = result->getSPs();
  while ((theSP_Constraint = theSP_Constraints()) != nullptr) {
    numSPs++;
  // theSP_Constraint->setDomain(0);
  }

  // mark the domain has having changed if numSPs > 0
  // as the constraint handlers have to be redone
  if (numSPs > 0)
    domain.domainChange();

  // finally return the load pattern
  return result;    
}    


LoadCase::PatternIterator &
LoadCase::getLoadPatterns()
{
  theLoadPatternIter->reset();
  return *theLoadPatternIter;
}


LoadPattern *
LoadCase::getLoadPattern(int tag) 
{
  TaggedObject *mc = theLoadPatterns->getComponentPtr(tag);
  // if not there return 0 otherwise perform a cast and return that  
  if (mc == 0) 
    return 0;
  LoadPattern *result = (LoadPattern *)mc;
  return result;
}


int 
LoadCase::getNumLoadPatterns() const
{
  return theLoadPatterns->getNumComponents();
}



void
LoadCase::applyLoad(double scale)
{
#if 0 // CMP: Keep this part in Domain
  // set the current pseudo time in the domain to be newTime
  currentTime = scale;
  dT = currentTime - committedTime;
  ops_Dt = dT;

  //
  // first zero all loads
  //
  Node *nodePtr;
  NodeIter &theNodeIter = domain.getNodes();
  while ((nodePtr = theNodeIter()) != nullptr)
    nodePtr->zeroUnbalancedLoad();

  Element *elePtr;
  ElementIter &theElemIter = domain.getElements();    
  while ((elePtr = theElemIter()) != nullptr)
    if (elePtr->isSubdomain() == false)
        elePtr->zeroLoad();    
#endif 

  //
  // now loop over load patterns, invoking applyLoad on them
  //
  LoadPattern *thePattern;
  PatternIterator &thePatterns = this->getLoadPatterns();
  while((thePattern = thePatterns()) != nullptr)
    thePattern->applyLoad(scale);

#if 0 // CMP: Keep this part in Domain
  //
  // finally loop over the MP_Constraints and SP_Constraints
  //

  MP_ConstraintIter &theMPs = this->getMPs();
  MP_Constraint *theMP;
  while ((theMP = theMPs()) != nullptr)
    theMP->applyConstraint(scale);

  SP_ConstraintIter &theSPs = this->getSPs();
  SP_Constraint *theSP;
  while ((theSP = theSPs()) != nullptr)
    theSP->applyConstraint(scale);
#endif
}


void
LoadCase::setLoadConstant()
{
  // loop over all the load patterns that are currently added to the domain
  // getting them to set their loads as now constant
  LoadPattern *thePattern;
  PatternIterator &thePatterns = this->getLoadPatterns();
  while((thePattern = thePatterns()) != nullptr)
    thePattern->setLoadConstant();
}


void
LoadCase::unsetLoadConstant()
{
  // loop over all the load patterns that are currently added to the domain
  // getting them to set their loads as now constant
  LoadPattern *thePattern;
  PatternIterator &thePatterns = this->getLoadPatterns();
  while((thePattern = thePatterns()) != nullptr)
    thePattern->unsetLoadConstant();
}

