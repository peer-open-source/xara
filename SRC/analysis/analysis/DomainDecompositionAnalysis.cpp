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
// Description: This file contains the class definition for 
// DomainDecompositionAnalysis. DomainDecompositionAnalysis is a subclass 
// of AnalysisAnalysis, it is used to perform the static condensation process
// on a subdomain.
//
// File: ~/analysis/Analysis/DomainDecompositionAnalysis.C
// 
// Written: fmk 
// Created: Tue Sept 17 16:34:47: 1996
// Revision: A
//
#include <DomainDecompositionAnalysis.h>
#include <ConstraintHandler.h>
#include <AnalysisModel.h>
#include <numberer/DOF_Numberer.h>
#include <LinearSOE.h>
#include <DomainDecompAlgo.h>
#include <DomainSolver.h>
#include <ConvergenceTest.h>
#include <IncrementalIntegrator.h>
#include <Subdomain.h>

#include <FE_Element.h>
#include <DOF_Group.h>
#include <Matrix.h>
#include <ID.h>
#include <Node.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>

DomainDecompositionAnalysis::DomainDecompositionAnalysis(Subdomain &the_Domain)
 : Analysis(the_Domain),
   MovableObject(DomDecompANALYSIS_TAGS_DomainDecompositionAnalysis),
   theSubdomain(&the_Domain),
   theHandler(nullptr),
   theNumberer(nullptr),
   theModel(nullptr),
   theAlgorithm(nullptr),
   theIntegrator(nullptr),
   theSOE(nullptr),
   theSolver(nullptr),
   theResidual(0),numEqn(0),numExtEqn(0),tangFormed(false),tangFormedCount(0),
   domainStamp(0),
   myChannel(0)
{
    theSubdomain->setDomainDecompAnalysis(*this);
}


DomainDecompositionAnalysis::DomainDecompositionAnalysis(int clsTag,
                                                         Subdomain &the_Domain)
 : Analysis(the_Domain),
   MovableObject(clsTag),
   theSubdomain(&the_Domain),
   theHandler(0),
   theNumberer(0),
   theModel(0),
   theAlgorithm(0),
   theIntegrator(0),
   theSOE(0),
   theSolver(0),
   theResidual(0),numEqn(0),numExtEqn(0),tangFormed(false),tangFormedCount(0),
   domainStamp(0),
   myChannel(0)
{

}

DomainDecompositionAnalysis::DomainDecompositionAnalysis(Subdomain &the_Domain,
                                                         ConstraintHandler &handler,
                                                         DOF_Numberer &numberer,
                                                         AnalysisModel &model,
                                                         DomainDecompAlgo &theSolnAlgo,
                                                         IncrementalIntegrator &integrator,
                                                         LinearSOE &theLinSOE,
                                                         DomainSolver &theDDSolver,
                                                         ConvergenceTest *theTest)


:Analysis(the_Domain),
 MovableObject(DomDecompANALYSIS_TAGS_DomainDecompositionAnalysis),
 theSubdomain( &the_Domain),
 theHandler( &handler),
 theNumberer( &numberer),
 theModel( &model),
 theAlgorithm( &theSolnAlgo),
 theIntegrator( &integrator),
 theSOE( &theLinSOE),
 theSolver( &theDDSolver),
 theResidual(0),numEqn(0),numExtEqn(0),tangFormed(false),tangFormedCount(0)
{
    theModel->setLinks(the_Domain, handler);
    theHandler->setLinks(*theModel);
    theNumberer->setLinks(*theModel);
    theIntegrator->setLinks(*theModel,*theSOE, theTest);
    theAlgorithm->setLinks(*theIntegrator,*theSOE,*theSolver,*theSubdomain);
    theSubdomain->setDomainDecompAnalysis(*this);
}    


DomainDecompositionAnalysis::~DomainDecompositionAnalysis()
{
  if (theResidual != 0)
    delete theResidual;
}    

void
DomainDecompositionAnalysis::clearAll()
{
    // invoke the destructor on all the objects in the aggregation
  if (theModel != 0)
    delete theModel;
  if (theHandler != 0)
    delete theHandler;
  if (theNumberer != 0)
    delete theNumberer;
  if (theIntegrator != 0)
    delete theIntegrator;
  if (theAlgorithm != 0)
    delete theAlgorithm;
  if (theSOE != 0)
    delete theSOE;
  
  // now set the pointers to NULL
  theModel =0;
  theHandler =0;
  theNumberer =0;
  theIntegrator =0;
  theAlgorithm =0;
  theSOE =0;
}    

int 
DomainDecompositionAnalysis::analyze(double dT)
{
    return 0;
}

int 
DomainDecompositionAnalysis::initialize()
{
    return 0;
}


bool
DomainDecompositionAnalysis::doesIndependentAnalysis()
{
    return false;
}

#if 0
int
DomainDecompositionAnalysis::domainChanged()
{
    // remove existing FE_elements and DOF_Groups from the Analysis
    theModel->clearAll();
    theHandler->clearAll();

    // now we invoke handle() on the constraint handler which
    // causes the creation of FE_Element and DOF_Group objects
    // and their addition to the AnalysisModel.

    numExtEqn = theHandler->handle(&(theSubdomain->getExternalNodes()));

    // we now get a node to number last

    const ID &theExtNodes = theSubdomain->getExternalNodes();
    int idSize = theExtNodes.Size();
    //    int theLastDOF = -1;

    ID theLastDOFs(1);
    int cnt = 0;

    // create an ID containing the tags of the DOF_Groups that are to
    // be numbered last
    for (int i=0; i<idSize; i++) {
        int nodeTag = theExtNodes(i);
        Node *nodePtr = theSubdomain->getNode(nodeTag);
        DOF_Group *dofGrpPtr = nodePtr->getDOF_GroupPtr();
        if (dofGrpPtr != nullptr) {
            const ID theID = dofGrpPtr->getID();
            int size = theID.Size();
            for (int j=0; j<size; j++)
                if (theID(j) == -3) {
                    theLastDOFs[cnt]  = dofGrpPtr->getTag();
                    cnt++;
                    j = size;
                }
        }
    }

    // we now invoke number() on the numberer which causes
    // equation numbers to be assigned to all the DOFs in the
    // AnalysisModel.    

    theNumberer->numberDOF(theLastDOFs);

    /*************************
    for (int i=0; i<idSize; i++) {
        int nodeTag = theExtNodes(i);
        Node *nodePtr = theSubdomain->getNode(nodeTag);
        DOF_Group *dofPtr = nodePtr->getDOF_GroupPtr();
        if (dofPtr != 0) {
            const ID theID = dofPtr->getID();
            int size = theID.Size();
            for (int j=0; j<size; j++)
                if (theID(j) == -3) {
                    theLastDOF = dofPtr->getTag();
                    i = idSize;
                    j=size;
                }
        }
    }
    theNumberer->numberDOF(theLastDOF);
    **********************/


    // we invoke setSize() on the LinearSOE which
    // causes that object to determine its size    
    
    theSOE->setSize(theModel->getDOFGraph());    
    numEqn = theSOE->getNumEqn();

    // we invoke domainChange() on the integrator and algorithm

    theIntegrator->domainChanged();
    //theAlgorithm->domainChanged();        

    // now set the variables to indicate that tangent has not been formed

    tangFormed = false;
    tangFormedCount = 0;
    
    return 0;
}
#endif

int
DomainDecompositionAnalysis::getNumExternalEqn()
{
  return numExtEqn;
}

int
DomainDecompositionAnalysis::getNumInternalEqn()
{
  return numEqn-numExtEqn;
}

#if 0
int  
DomainDecompositionAnalysis::newStep(double dT)
{
  return theIntegrator->newStep(dT);
}

int  
DomainDecompositionAnalysis::analysisStep(double dT)
{
  return theIntegrator->newStep(dT);
}
#endif

int  
DomainDecompositionAnalysis::eigenAnalysis(int numMode, bool generalized, bool findSmallest)
{
  opserr << "DomainDecompositionAnalysis::eigenAnalysis() - should not be called\n";
  return -1;;
}

int  
DomainDecompositionAnalysis::setEigenSOE(EigenSOE &theSOE)
{
  opserr << "DomainDecompositionAnalysis::setEigenSOE() - should not be called\n";
  return -1;;
}

int  
DomainDecompositionAnalysis::computeInternalResponse()
{
  return theAlgorithm->solveCurrentStep();
}


int  
DomainDecompositionAnalysis::formTangent()
{
    int result =0;

    Domain *the_Domain = this->getDomainPtr();

    // we check to see if the domain has changed 
    int stamp = the_Domain->hasDomainChanged();
    if (stamp != domainStamp) {
        domainStamp = stamp;
        this->domainChanged();
    }
    
    // if tangFormed == -1 then formTangent has already been
    // called for this state by formResidual() or formTangVectProduct()
    // so we won't be doing it again.

    if (tangFormedCount != -1) {
        result = theIntegrator->formTangent();
        if (result < 0)
            return result;
        result = theSolver->condenseA(numEqn-numExtEqn);
        if (result < 0)
            return result;
    }
        
    tangFormed = true;
    tangFormedCount++;
    
    return result;
}



int  
DomainDecompositionAnalysis::formResidual(Vector& G)
{
    int result =0;
    Domain *the_Domain = this->getDomainPtr();    
    
    // we check to see if the domain has changed 
    int stamp = the_Domain->hasDomainChanged();
    if (stamp != domainStamp) {
        domainStamp = stamp;
        this->domainChanged();
    }
    
    if (tangFormed == false) {
        result = this->formTangent();
        if (result < 0)
            return result;
        tangFormedCount = -1; // set to minus number so tangent 
                              // is not formed twice at same state
    }

    result = theIntegrator->formUnbalance(G);

    if (result < 0)
        return result;

    return theSolver->condenseRHS(numEqn-numExtEqn);
}



int  
DomainDecompositionAnalysis::formTangVectProduct(Vector &u)
{
    int result = 0;

    Domain *the_Domain = this->getDomainPtr();
    
    // we check to see if the domain has changed 
    int stamp = the_Domain->hasDomainChanged();
    if (stamp != domainStamp) {
        domainStamp = stamp;
        this->domainChanged();
    }
    
    if (tangFormed == false) {
        result = this->formTangent();
        if (result < 0)
            return result;
        tangFormedCount = -1; // set to minus number so tangent 
                              // is not formed twice at same state
    }    

    return theSolver->computeCondensedMatVect(numEqn-numExtEqn,u);
}



const Matrix &
DomainDecompositionAnalysis::getTangent()
{
    Domain *the_Domain = this->getDomainPtr();
    
    // we check to see if the domain has changed 
    int stamp = the_Domain->hasDomainChanged();
    if (stamp != domainStamp) {
        domainStamp = stamp;
        this->domainChanged();
    }

    if (tangFormed == false) {
        this->formTangent();
    }
    
    return (theSolver->getCondensedA());
}



const Vector &
DomainDecompositionAnalysis::getResidual()
{

    Domain *the_Domain = this->getDomainPtr();
    
    // we check to see if the domain has changed 
    int stamp = the_Domain->hasDomainChanged();
    if (stamp != domainStamp) {
        domainStamp = stamp;
        this->domainChanged();
        // this->formResidual();
    }
    
    if (theResidual == 0) {
        theResidual = new Vector(theSolver->getCondensedRHS());
        return *theResidual;        
    }
    else if (theResidual->Size() != numExtEqn) {
        delete theResidual;
        theResidual = new Vector(theSolver->getCondensedRHS());
        return *theResidual;                    
    }
    else {
        (*theResidual) = theSolver->getCondensedRHS();
    }

    return *theResidual;
}




const Vector &
DomainDecompositionAnalysis::getTangVectProduct()
{
    Domain *the_Domain = this->getDomainPtr();
    
    // we check to see if the domain has changed 
    int stamp = the_Domain->hasDomainChanged();
    if (stamp != domainStamp) {
        domainStamp = stamp;
        this->domainChanged();
    }
    
    return theSolver->getCondensedMatVect();
}





Subdomain  *
DomainDecompositionAnalysis::getSubdomainPtr() const
{
  return theSubdomain;
}




ConstraintHandler *
DomainDecompositionAnalysis::getConstraintHandlerPtr() const
{
  return theHandler;
}



DOF_Numberer *
DomainDecompositionAnalysis::getDOF_NumbererPtr() const
{
    return theNumberer;
}



AnalysisModel  *
DomainDecompositionAnalysis::getAnalysisModelPtr() const
{
    return theModel;
}



DomainDecompAlgo  *
DomainDecompositionAnalysis::getDomainDecompAlgoPtr() const
{
    return theAlgorithm;
}



IncrementalIntegrator *
DomainDecompositionAnalysis::getIncrementalIntegratorPtr() const
{
    return theIntegrator;
}



LinearSOE *
DomainDecompositionAnalysis::getLinSOEPtr() const
{
    return theSOE;    
}



DomainSolver *
DomainDecompositionAnalysis::getDomainSolverPtr() const
{
    return theSolver;
}



int 
DomainDecompositionAnalysis::setAlgorithm(EquiSolnAlgo &theAlgorithm)
{
  opserr << "DomainDecompositionAnalysis::setAlgorithm() - not implemented\n";
  return -1;
}

int 
DomainDecompositionAnalysis::setIntegrator(IncrementalIntegrator &theIntegrator) 
{
  opserr << "DomainDecompositionAnalysis::setIntegrator() - not implemented\n";
  return -1;
}


int 
DomainDecompositionAnalysis::setLinearSOE(LinearSOE &theSOE)
{
  opserr << "DomainDecompositionAnalysis::setLinearSOE() - not implemented\n";
  return -1;
}

int 
DomainDecompositionAnalysis::setConvergenceTest(ConvergenceTest &theTest)
{
  opserr << "DomainDecompositionAnalysis::setConvergenceTest() - not implemented\n";
  return -1;
}


int 
DomainDecompositionAnalysis::checkAllResult(int mine)
{
  static ID data(1);
  data(0) = mine;
  if (myChannel != 0) {
    myChannel->sendID(0,0,data);
    myChannel->recvID(0,0,data);
  }
  return data(0);
}
