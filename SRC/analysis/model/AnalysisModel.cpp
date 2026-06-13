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
// Purpose: This file contains the class definition for AnalysisModel
// AnalysisModel is a container class. The class is responsible for holding
// and providing access to the Elements, Nodes, LoadCases, SP_Constraints 
// and MP_Constraints. These objects are all added to the AnalysisModel by a 
// ModelBuilder.
//
// Written: fmk 
// Created: Fri Sep 20 15:27:47: 1996
// Revision: A
//
#include <assert.h>
#include <cmath>

#include <AnalysisModel.h>
#include <Domain.h>
#include <FE_Element.h>
#include <DOF_Group.h>
#include <Graph.h>
#include <Vertex.h>
#include <Node.h>
#include <NodeIter.h>
#include <ConstraintHandler.h>

#include <Vector.h>
#include <Matrix.h>
#include <LinearSOE.h>
#include <Node.h>
#include <NodeIter.h>
#include <LoadPattern.h>
#include <Parameter.h>
#include <ParameterIter.h>

#include <MapOfTaggedObjects.h>
#include <VectorOfTaggedObjects.h>

#include <analysis/damping/ModalDamping.h>

#define START_EQN_NUM 0
#define START_VERTEX_NUM 0


static int 
ApproxDOF(Domain& domain)
{
  int nf = 0;
  Node *nodePtr;
  NodeIter &theNodes = domain.getNodes();
  while ((nodePtr = theNodes()) != nullptr)
    nf += nodePtr->getNumberDOF();

  return nf;
}


AnalysisModel::AnalysisModel(Domain& domain)
: MovableObject(AnaMODEL_TAGS_AnalysisModel),
  myDomain(&domain)
 , myHandler(nullptr)
 , myDOFGraph(0)
 , myGroupGraph(0)
 , numFE_Ele(0), numDOF_Grp(0), numEqn(0)
 , eigenVectors(0), eigenValues(0)
 , modalDamping(nullptr)
{
#if 0
  theFEs     = new VectorOfTaggedObjects(); // 256);
  theDOFs    = new VectorOfTaggedObjects(); // 256);
#else 

  theFEs     = new ArrayOfTaggedObjects(1024);
  theDOFs    =  new ArrayOfTaggedObjects(1024);
#endif
  theFEiter  = new FE_EleIter(theFEs);
  theDOFiter = new DOF_GrpIter(theDOFs);

  theDOFs->setSize(ApproxDOF(domain));
  theFEs->setSize(domain.getNumElements());
} 



AnalysisModel::~AnalysisModel()
{
  if (theFEs != nullptr) {
    theFEs->clearAll();
    delete theFEs;
  }

  if (theDOFs != nullptr) {
    theDOFs->clearAll();
    delete theDOFs;
  }

  if (theFEiter != nullptr)
    delete theFEiter;

  if (theDOFiter != nullptr)
    delete theDOFiter;

  if (myGroupGraph != nullptr) {
    delete myGroupGraph;    
  }        
  
  if (myDOFGraph != nullptr) {
    delete myDOFGraph;
  }

  if (modalDamping != nullptr)
    delete modalDamping;
}    

void
AnalysisModel::setLinks(Domain &theDomain, ConstraintHandler &theHandler)
{
  myDomain  = &theDomain;
  myHandler = &theHandler;
}


bool
AnalysisModel::addFE_Element(FE_Element *theElement)
{
  // check we don't add a null pointer or this is a subclass
  // trying to use this method when it shouldn't
  if (theElement == 0 || theFEs == 0)
    return false;

  // check if an Element with a similar tag already exists in the Domain
  int tag = theElement->getTag();
  TaggedObject *other = theFEs->getComponentPtr(tag);
  if (other != 0) {
    opserr << "AnalysisModel::addFE_Element - element with tag " 
           << tag << " already exists in model\n"; 
    return false;
  }

  // add the element to the container object for the elements
  bool result = theFEs->addComponent(theElement);
  if (result == true) {
    theElement->setAnalysisModel(*this);
    numFE_Ele++;
    return true;  // o.k.
  }

  return result;
}



bool
AnalysisModel::addDOF_Group(DOF_Group *theGroup)
{
  // check we don't add a null pointer or this is a subclass trying
  // to use a method it shouldn't be using
  if (theGroup == 0 || theDOFs == 0)
      return false;

  // check if an Element with a similar tag already exists in the Domain
  int tag = theGroup->getTag();
  TaggedObject *other = theDOFs->getComponentPtr(tag);
  if (other != nullptr) {
    opserr << "AnalysisModel::addDOF_Group - group with tag "
           << tag << " already exists in model\n"; 
    return false;
  }

  // add the element to the container object for the elements
  bool result = theDOFs->addComponent(theGroup);
  if (result == true) {
    numDOF_Grp++;
    return true; // o.k.
  } else
    return false;
}


void
AnalysisModel::clearAll() 
{
  // if the graphs have been constructed delete them
  this->clearDOFGroupGraph();

  this->clearDOFGraph(); 

  theFEs->clearAll();
  theDOFs->clearAll();
  
  numFE_Ele =0;
  numDOF_Grp = 0;
  numEqn = 0;

  if (myHandler != nullptr)
    myHandler->clearAll();
  
  if (modalDamping != nullptr)
    delete modalDamping;
  modalDamping = nullptr;

  
  // for the nodes reset the DOF_Group pointers to 0
  Domain *theDomain = this->getDomainPtr();
  if (theDomain == nullptr)
    return;
  
  NodeIter &theNod = theDomain->getNodes();
  Node *nodPtr;
  while ((nodPtr = theNod()) != nullptr)
    nodPtr->setDOF_GroupPtr(nullptr);
}

void
AnalysisModel::clearDOFGraph() 
{
  if (myDOFGraph != nullptr)
    delete myDOFGraph;

  myDOFGraph = nullptr;
}

void
AnalysisModel::clearDOFGroupGraph() 
{
  if (myGroupGraph != nullptr)
    delete myGroupGraph;    
  
  myGroupGraph = nullptr;
}




int
AnalysisModel::getNumDOF_Groups() const
{
  return numDOF_Grp;
}


DOF_Group *
AnalysisModel::getDOF_GroupPtr(int tag)
{
  TaggedObject *other = theDOFs->getComponentPtr(tag);
  if (other == nullptr) {
    return nullptr;
  }
  DOF_Group *result = (DOF_Group *)other;
  return result;
}


FE_EleIter &
AnalysisModel::getFEs()
{
  theFEiter->reset();
  return *theFEiter;
}

DOF_GrpIter &
AnalysisModel::getDOFs()
{
  theDOFiter->reset();
  return *theDOFiter;
}

void 
AnalysisModel::setNumEqn(int theNumEqn)
{
  if (modalDamping != nullptr && (numEqn != theNumEqn)) {
    delete modalDamping;
    modalDamping = nullptr;
  }

  numEqn = theNumEqn;
}

int 
AnalysisModel::getNumEqn() const
{
  return numEqn;
}


Graph &
AnalysisModel::getDOFGraph()
{
  if (myDOFGraph == nullptr) {
    // int numVertex = this->getNumDOF_Groups();
    //    myDOFGraph = new Graph(numVertex);
    MapOfTaggedObjects *graphStorage = new MapOfTaggedObjects();
    myDOFGraph = new Graph(*graphStorage);

    //
    // create a vertex for each dof
    //
    
    DOF_Group *dofPtr =0;
    DOF_GrpIter &theDOFs = this->getDOFs();
    while ((dofPtr = theDOFs()) != 0) {
      const ID &id = dofPtr->getID();
      int size = id.Size();
      for (int i=0; i<size; i++) {
        int dofTag = id(i);
        if (dofTag >= START_EQN_NUM) {
          Vertex *vertexPtr = myDOFGraph->getVertexPtr(dofTag);
          if (vertexPtr == 0) {
            Vertex *vertexPtr = new Vertex(dofTag, dofTag);      
            if (vertexPtr == 0) {
              opserr << "WARNING AnalysisModel::getDOFGraph";
              opserr << " - Not Enough Memory to create " << i+1 << "th Vertex\n";
              return *myDOFGraph;
            }
            if (myDOFGraph->addVertex(vertexPtr, false) == false) {
              opserr << "WARNING AnalysisModel::getDOFGraph - error adding vertex\n";
              return *myDOFGraph;
            }
          }
        }
      }
    }
    
    // now add the edges, by looping over the FE_elements, getting their
    // IDs and adding edges between DOFs for equation numbers >= START_EQN_NUM
    
    FE_Element *elePtr =0;
    FE_EleIter &eleIter = this->getFEs();
    
    myDOFGraph->startAddEdge();
    while((elePtr = eleIter()) != 0) {
      const ID &id = elePtr->getID();
      const int size = id.Size();
      for (int i=0; i<size; i++) {
        int eqn1 = id(i);
        
        // if eqnNum of DOF is a valid eqn number add an edge
        // to all other DOFs with valid eqn numbers.
        if (eqn1 >=START_EQN_NUM) {
          for (int j=i+1; j<size; j++) {
            int eqn2 = id(j);
            if (eqn2 >=START_EQN_NUM)
              myDOFGraph->addEdgeFast(eqn1-START_EQN_NUM+START_VERTEX_NUM,
                                  eqn2-START_EQN_NUM+START_VERTEX_NUM);
          }
        }
      }
    }
  }    

  return *myDOFGraph;
}


Graph &
AnalysisModel::getDOFGroupGraph()
{
  if (myGroupGraph == nullptr) {
    // int numVertex = this->getNumDOF_Groups();
    // assert(numVertex != 0);
    // myGroupGraph = new Graph(numVertex);
    MapOfTaggedObjects *graphStorage = new MapOfTaggedObjects();
    myGroupGraph = new Graph(*graphStorage);
        

    // now create the vertices with a reference equal to the DOF_Group number.
    // and a tag which ranges from 0 through numVertex-1
    DOF_Group   *dofPtr;
    DOF_GrpIter &dofIter2 = this->getDOFs();
    // int count = START_VERTEX_NUM;
    while ((dofPtr = dofIter2()) != 0) {
        int DOF_GroupTag = dofPtr->getTag();
        int DOF_GroupNodeTag = dofPtr->getNodeTag();
        int numDOF = dofPtr->getNumFreeDOF();
        Vertex *vertexPtr = new Vertex(DOF_GroupTag, DOF_GroupNodeTag, 0, numDOF);

        myGroupGraph->addVertex(vertexPtr);
    }

    // now add the edges, by looping over the Elements, getting their
    // IDs and adding edges between DOFs for equation numbers >= START_EQN_NUM
    
    FE_Element *elePtr;
    FE_EleIter &eleIter = this->getFEs();

    while((elePtr = eleIter()) != 0) {
      const ID &id = elePtr->getDOFtags();
      int size = id.Size();
      for (int i=0; i<size; i++) {
        int dof1 = id(i);
        for (int j=0; j<size; j++) 
          if (i != j) {
            int dof2 = id(j);
            myGroupGraph->addEdge(dof1,dof2);
          }
      }
    }
  }

  return *myGroupGraph;
}




void 
AnalysisModel::setResponse(const Vector &disp,
                           const Vector &vel, 
                           const Vector &accel)
{
  DOF_GrpIter &theDOFGrps = this->getDOFs();
  DOF_Group   *dofPtr;

  while ((dofPtr = theDOFGrps()) != nullptr) {
    dofPtr->setNodeDisp(disp);
    dofPtr->setNodeVel(vel);
    dofPtr->setNodeAccel(accel);
  }
}


void
AnalysisModel::setStateGradient(
                    const Vector &u,
                    const Vector &v, 
                    const Vector &a,
                    int grad, int ngrad)
{
  DOF_GrpIter &theDOFGrps = this->getDOFs();
  DOF_Group         *group;

  while ((group = theDOFGrps()) != nullptr) {
    group->saveSensitivity(u, v, a, grad, ngrad);
  }
}

void 
AnalysisModel::getStateGradient(
                    Vector &du, 
                    Vector &dv, 
                    Vector &da, 
                    int grad)
{

  DOF_GrpIter &theDOFs = this->getDOFs();
  DOF_Group *group;
  while ((group = theDOFs()) != nullptr) {

    const ID &id = group->getID();
    const int idSize = id.Size();
    const Vector &dispSens = group->getDispSensitivity(grad);
    for (int i = 0; i < idSize; i++)
      if (int loc = id(i); loc >= 0)
        du(loc) = dispSens(i);

    const Vector &velSens = group->getVelSensitivity(grad);
    for (int i = 0; i < idSize; i++)
      if (int loc = id(i); loc >= 0)
        dv(loc) = velSens(i);

    const Vector &accelSens = group->getAccSensitivity(grad);
    for (int i = 0; i < idSize; i++)
      if (int loc = id(i); loc >= 0)
        da(loc) = accelSens(i);
  }
}

int 
AnalysisModel::commitGradient(int gradNum, int numGrads)
{
  // Loop through the FE_Elements and set unconditional sensitivities
  FE_Element *elePtr;
  FE_EleIter &theEles = this->getFEs();    
  while ((elePtr = theEles()) != nullptr)
    elePtr->commitSensitivity(gradNum, numGrads);

  return 0;
}

void 
AnalysisModel::setDisp(const Vector &disp)
{
  DOF_GrpIter &theDOFGrps = this->getDOFs();
  DOF_Group   *dofPtr;

  while ((dofPtr = theDOFGrps()) != nullptr)
    dofPtr->setNodeDisp(disp);
}        
        
void 
AnalysisModel::setVel(const Vector &vel)
{
  DOF_GrpIter &theDOFGrps = this->getDOFs();
  DOF_Group   *dofPtr;
  
  while ((dofPtr = theDOFGrps()) != nullptr)
    dofPtr->setNodeVel(vel);
}
        

void 
AnalysisModel::setAccel(const Vector &accel)
{
  DOF_GrpIter &theDOFGrps = this->getDOFs();
  DOF_Group   *dofPtr;
  
  while ((dofPtr = theDOFGrps()) != 0) 
    dofPtr->setNodeAccel(accel);        
}

void 
AnalysisModel::incrDisp(const Vector &disp)
{
  DOF_GrpIter &theDOFGrps = this->getDOFs();
  DOF_Group         *dofPtr;

  while ((dofPtr = theDOFGrps()) != nullptr)
    dofPtr->incrNodeDisp(disp);
}
        
void 
AnalysisModel::incrVel(const Vector &vel)
{
  DOF_GrpIter &theDOFGrps = this->getDOFs();
  DOF_Group         *dofPtr;    
  while ((dofPtr = theDOFGrps()) != nullptr)
    dofPtr->incrNodeVel(vel);
}


#if 0
void 
AnalysisModel::incrAccel(const Vector &accel)
{
  DOF_GrpIter &theDOFGrps = this->getDOFs();
    DOF_Group         *dofPtr;
    
    while ((dofPtr = theDOFGrps()) != 0) 
      dofPtr->incrNodeAccel(accel);        
}        
#endif



int
AnalysisModel::getState(Vector &U, Vector &Udot, Vector &Udotdot, int flag)
{

  if (U.Size() != numEqn || Udot.Size() != numEqn || Udotdot.Size() != numEqn) {
    return -1;
  }

  DOF_GrpIter &theDOFs = this->getDOFs();
  DOF_Group *dofPtr;
  while ((dofPtr = theDOFs()) != 0)  {
    const ID &id = dofPtr->getID();
    int idSize = id.Size();

    const Vector &disp = dofPtr->getCommittedDisp();
    for (int i=0; i < idSize; i++)  {
        int loc = id(i);
        if (loc >= 0)  {
          U(loc) = disp(i);
        }
    }
    
    const Vector &vel = dofPtr->getCommittedVel();
    for (int i=0; i < idSize; i++)  {
      int loc = id(i);
      if (loc >= 0)  {
        Udot(loc) = vel(i);
      }
    }
    
    const Vector &accel = dofPtr->getCommittedAccel();
    for (int i=0; i < idSize; i++)  {
      int loc = id(i);
      if (loc >= 0)  {
        Udotdot(loc) = accel(i);
      }
    }
  }

  return 0;
}


int 
AnalysisModel::applyResidual(Integrator& assm, LinearSOE& soe)
{

  // loop through the DOF_Groups and add the unbalance
  // same as IncrementalIntegrator::formNodalUnbalance
  DOF_GrpIter &theDOFs = this->getDOFs();
  DOF_Group *dofPtr;
  int res = 0;

  while ((dofPtr = theDOFs()) != nullptr) {
    if (soe.addB( dofPtr->getUnbalance(&assm), dofPtr->getID()) <0) [[unlikely]] {
      res = -1;
    }
  }

  // loop through the FE_Elements and add the residual
  // same as IncrementalIntegrator::formElementResidual
  FE_Element *elePtr;
  FE_EleIter &theEles2 = this->getFEs();    
  while ((elePtr = theEles2()) != nullptr) {
    if (soe.addB(elePtr->getResidual(&assm), elePtr->getID()) < 0) [[unlikely]] {
      res = -2;
    }
  }

  return res;
}



int 
AnalysisModel::applyInertia(const Vector &v, Vector &res)
{
  // int n = v.Size();
  assert(v.Size() == numEqn);
  assert(res.Size() == numEqn);

  res.Zero();

  // loop over the FE_Elements
  FE_Element *elePtr;
  FE_EleIter &theEles = this->getFEs();    
  while((elePtr = theEles()) != nullptr) {
    const Vector &b = elePtr->getM_Force(v, 1.0);
    res.Assemble(b, elePtr->getID(), 1.0);
  }

  // loop over the DOF_Groups
  DOF_Group *dofPtr;
  DOF_GrpIter &theDofs = this->getDOFs();
  while ((dofPtr = theDofs()) != nullptr) {
    const Vector &a = dofPtr->getM_Force(v, 1.0);      
    res.Assemble(a, dofPtr->getID(), 1.0);
  }
  return 0;
}


int 
AnalysisModel::applyInertia(const Vector &v, LinearSOE& soe, double fact)
{
  assert(v.Size() == numEqn);
  assert(soe.getNumEqn() == numEqn);

  // loop over the FE_Elements
  FE_Element *elePtr;
  FE_EleIter &theEles = this->getFEs();    
  while((elePtr = theEles()) != nullptr) {
    const Vector &b = elePtr->getM_Force(v, 1.0);
    soe.addB(b, elePtr->getID(), fact);
  }

  // loop over the DOF_Groups
  DOF_Group *dofPtr;
  DOF_GrpIter &theDofs = this->getDOFs();
  while ((dofPtr = theDofs()) != 0) {
    const Vector &a = dofPtr->getM_Force(v, 1.0);      
    soe.addB(a, dofPtr->getID(), fact);
  }
  return 0;
}


void 
AnalysisModel::setNumEigenvectors(int numEigenvectors)
{
  Node *theNode;
  NodeIter &theNodes = myDomain->getNodes();
  while ((theNode = theNodes()) != 0)
    theNode->setNumEigenvectors(numEigenvectors);
}


void 
AnalysisModel::setEigenvalues(const Vector &eigenvalues)
{
  myDomain->setEigenvalues(eigenvalues);
}        


const Vector &
AnalysisModel::getEigenvalues()
{
  return myDomain->getEigenvalues();
}        



const Vector *
AnalysisModel::getModalDampingFactors()
{
  return myDomain->getModalDampingFactors();
}


bool 
AnalysisModel::inclModalDampingMatrix()
{
  return myDomain->inclModalDampingMatrix();
}

int 
AnalysisModel::setModalDamping(const Vector &modalDampingFactors)
{
  if (modalDamping != nullptr)
    delete modalDamping;
  modalDamping = new ModalDamping(*this, modalDampingFactors, this->getNumEqn());

  return 0;
}

void 
AnalysisModel::setEigenvector(int mode, const Vector &eigenvalue)
{
  DOF_GrpIter &theDOFGrps = this->getDOFs();
  DOF_Group   *dofPtr;
  
  while ((dofPtr = theDOFGrps()) != nullptr) 
    dofPtr->setEigenvector(mode, eigenvalue);        
}        


void 
AnalysisModel::applyLoadDomain(double time)
{
  assert(myDomain != nullptr);
  myDomain->applyLoad(time);
  myHandler->applyLoad();
}


int
AnalysisModel::applyLoadGradient() 
{
  // 1) Zero the unbalaced load
  Node *node;
  NodeIter &theNodeIter = myDomain->getNodes();
  while ((node = theNodeIter()) != nullptr)
    node->zeroUnbalancedLoad();

  // 2) Add external load sensitivity
  LoadPattern *pattern;
  auto &thePatterns = myDomain->getLoadPatterns();
  while ((pattern = thePatterns()) != nullptr)
    pattern->applyLoadSensitivity(myDomain->getCurrentTime());

  return 0;
}


int
AnalysisModel::updateDomain()
{
  assert(myDomain != nullptr);

  int res = myDomain->update();

  if (res == 0)
    return myHandler->update();

  return res;
}


int
AnalysisModel::updateDomain(double newTime, double dT)
{
  assert(myDomain != nullptr);

  // invoke the method
  int res = 0;
  myDomain->applyLoad(newTime);
  if (res == 0)
    res = myHandler->applyLoad();
  if (res == 0)
    res = myDomain->update();
  if (res == 0)
    res = myHandler->update();

  return res;
}


int
AnalysisModel::analysisStep(double dT)
{
  assert(myDomain != nullptr);
  return myDomain->analysisStep(dT);
}


int
AnalysisModel::commitDomain()
{
  assert(myDomain != nullptr);

  // commit the domain state
  if (myDomain->commit() < 0) {
    opserr << "WARNING: AnalysisModel::commitDomain - Domain::commit() failed\n";
    return -2;
  }

  return 0;
}

# if 0
int
AnalysisModel::revertDomainToLastCommit()
{
  assert(myDomain != nullptr);

  // invoke the method
  if (myDomain->revertToLastCommit() < 0) {
    opserr << "WARNING: AnalysisModel::revertDomainToLastCommit.";
    opserr << " Domain::revertToLastCommit() failed.\n";
    return -2;
  }
  return 0;
}
#endif


double
AnalysisModel::getCurrentDomainTime()
{
  assert(myDomain != nullptr);
  return myDomain->getCurrentTime();
}


void
AnalysisModel::setCurrentDomainTime(double newTime)
{
  assert(myDomain != nullptr);
  myDomain->setCurrentTime(newTime);
}



void
AnalysisModel::setRayleighDampingFactors(double alphaM, double betaK, double betaK0, double betaKc)
{
  assert(myDomain != nullptr);
  myDomain->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);
}




Domain *
AnalysisModel::getDomainPtr() const
{
  return myDomain;
}

void 
AnalysisModel::Print(OPS_Stream &s, int flag)
{
  opserr << "{\n";
  opserr << "  \"u\": [\n";
  Node* theNode;
  NodeIter &theNodes = myDomain->getNodes();
  bool firstNode = true;
  while ((theNode = theNodes()) != nullptr) {
    if (firstNode)
      firstNode = false;
    else
      opserr << ",\n";
    opserr << "    [";
    const Vector &disp = theNode->getDisp();
    for (int i = 0; i < disp.Size(); i++) {
      if (i != 0)
        opserr << ", ";
      opserr << disp(i);
    }
    opserr << "]";
  }
  opserr << "\n  ]\n";
  opserr << "}\n";
}

int
AnalysisModel::sendSelf(int cTag, Channel &theChannel)
{
  return 0;
}


int
AnalysisModel::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker) 
{
  return 0;
}

