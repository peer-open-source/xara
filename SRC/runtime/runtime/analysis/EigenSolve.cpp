#include "EigenSolve.h"
#include <Logging.h>
#include <Domain.h>
#include <AnalysisModel.h>
#include <ConstraintHandler.h>
#include <Graph.h>
#include <FE_Element.h>
#include <DOF_Group.h>
#include <LinearSOE.h>
#include <DOF_Numberer.h>

// For eigen
#include <FE_Element.h>
#include <DOF_Group.h>
#include <TransformationConstraintHandler.h>

// Default concrete analysis classes
#include <Newmark.h>
#include <EigenSOE.h>
#include <SymBandEigenSolver.h>
#include <SymBandEigenSOE.h>
#include <FullGenEigenSolver.h>
#include <FullGenEigenSOE.h>
#include <ArpackSOE.h>
#include <SymmGeneralizedEigenSOE.h>
#include <SymmGeneralizedEigenSolver.h>




void
EigenSolve::setup(AnalysisModel& theAnalysisModel,
                  LinearSOE* theSOE,
                  ConstraintHandler* theHandler,
                  DOF_Numberer* theNumberer,
                  const Options& options)
{
  EigenSOE *theEigenSOE = nullptr;

  if (theHandler == nullptr)
    theHandler = new TransformationConstraintHandler();

  if (theNumberer == nullptr)
    theNumberer = new DOF_Numberer(*(new RCM(false)));

  if (theSOE == nullptr)
    // TODO: CHANGE TO MORE GENERAL SOE
    theSOE = new ProfileSPDLinSOE(*(new ProfileSPDLinDirectSolver()));


  this->fillDefaults(EIGEN_ANALYSIS);
  this->setLinks(EIGEN_ANALYSIS);

  if (theEigenSOE == nullptr) {
    domainStamp = 0;
    if (typeSolver == EigenSOE_TAGS_SymBandEigenSOE) {
      SymBandEigenSolver *theEigenSolver = new SymBandEigenSolver();
      theEigenSOE = new SymBandEigenSOE(*theEigenSolver, theAnalysisModel);
    }
    else if (typeSolver == EigenSOE_TAGS_FullGenEigenSOE) {
      FullGenEigenSolver *theEigenSolver = new FullGenEigenSolver();
      theEigenSOE = new FullGenEigenSOE(*theEigenSolver, theAnalysisModel);
    }
    else if (typeSolver == EigenSOE_TAGS_SymmGeneralizedEigenSOE) {
      SymmGeneralizedEigenSolver *theEigenSolver = new SymmGeneralizedEigenSolver();
      theEigenSOE = new SymmGeneralizedEigenSOE(*theEigenSolver, theAnalysisModel);
    }
    else {
      theEigenSOE = new ArpackSOE(theAnalysisModel, shift);
    }

    //
    // set the eigen soe in the system
    //
    theEigenSOE->setLinks(theAnalysisModel);
    theEigenSOE->setLinearSOE(*theSOE);
  }
}


int
EigenSolve::solve(AnalysisModel& model,
                  LinearSOE& soe,
                  ConstraintHandler& handler,
                  const Options& options)
{
  return 0;
}


int
EigenSolve::eigen(AnalysisModel& theAnalysisModel,
                  LinearSOE& theSOE,
                  EigenSOE& theEigenSOE,
                  ConstraintHandler& theHandler)
{

  int result = 0;
  Domain *the_Domain = theAnalysisModel.getDomainPtr();

  // for parallel processing, want all analysis doing an eigenvalue analysis
  result = the_Domain->eigenAnalysis(numMode, generalized, findSmallest);

  int stamp = the_Domain->hasDomainChanged();

  if (stamp != domainStamp) {
    //domainStamp = stamp; // commented out so domainChanged() gets called with integrator,
                         //  which isnt updated here
//    result = this->domainChanged();

    this->number();

    Graph &theGraph = theAnalysisModel.getDOFGraph();

    result = theSOE.setSize(theGraph);

    result = theEigenSOE.setSize(theGraph);

    theAnalysisModel.clearDOFGraph();

    if (result < 0) {
      opserr << "EigenSolve::eigen() - domainChanged failed\n";
      return -1;
    }
  }
#if 0
  theEigenSOE->formSystem(*theAnalysisModel, generalized);
#else
  //
  // zero A and M
  //
  theEigenSOE.zeroA();
  theEigenSOE.zeroM();

  //
  // form K
  //
  FE_EleIter &theEles = theAnalysisModel.getFEs();
  FE_Element *elePtr;
  while ((elePtr = theEles()) != nullptr) {
    elePtr->zeroTangent();
    elePtr->addKtToTang(1.0);
    if (theEigenSOE.addA(elePtr->getTangent(nullptr), elePtr->getID()) < 0) {
      opserr << G3_WARN_PROMPT << "eigen -";
      opserr << " failed in addA for ID " << elePtr->getID();
      result = -2;
    }
  }

  //
  // If generalized is true, form M
  //
  if (generalized == true) {
    FE_EleIter &theEles2 = theAnalysisModel.getFEs();
    while ((elePtr = theEles2()) != nullptr) {
      elePtr->zeroTangent();
      elePtr->addMtoTang(1.0);
      if (theEigenSOE.addM(elePtr->getTangent(0), elePtr->getID()) < 0) {
        opserr << G3_WARN_PROMPT << "eigen -";
        opserr << " failed in addA for ID " << elePtr->getID() << "\n";
        result = -2;
      }
    }

    DOF_Group *dofPtr;
    DOF_GrpIter &theDofs = theAnalysisModel.getDOFs();
    while ((dofPtr = theDofs()) != nullptr) {
      dofPtr->zeroTangent();
      dofPtr->addMtoTang(1.0);
      if (theEigenSOE.addM(dofPtr->getTangent(0), dofPtr->getID()) < 0) {
        opserr << G3_WARN_PROMPT << "theEigenSOE failed in addM for ID " << dofPtr->getID() << "\n";
        result = -3;
      }
    }
  }
#endif

  //
  // Solve for the eigen values & vectors
  //
  if (theEigenSOE.solve(numMode, generalized, findSmallest) < 0) {
    return -4;
  }

  //
  // Store the eigenvalues and eigenvectors in the model
  //
  theAnalysisModel.setNumEigenvectors(numMode);
  Vector theEigenvalues(numMode);
  for (int i = 1; i <= numMode; i++) {
    theEigenvalues[i-1] = theEigenSOE.getEigenvalue(i);
    theAnalysisModel.setEigenvector(i, theEigenSOE.getEigenvector(i));
  }
  theAnalysisModel.setEigenvalues(theEigenvalues);
  this->numEigen = numMode;


  //
  delete theEigenSOE;
  theEigenSOE = nullptr;
  return 0;
}

