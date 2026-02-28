
// $Revision: 1.5 $
// $Date: 2009-05-14 22:45:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/eigenSOE/EigenSOE.cpp,v $

// Written: Jun Peng
// Created: Sat Feb. 6, 1999
// Revision: A
//
// Description: This file contains the class definition of EigenSOE.
// EigenSOE is a subclass of SystemOfEqn.
// It has pure virtual functions which must be implemented in it's derived
// subclasses. To solve the genreal eigen value equations means that
// by the given K and M, find the corresponding eigen value and eigen
// vectors.


#include <EigenSOE.h>
#include <EigenSolver.h>
#include <AnalysisModel.h>

// for formSystem
#include <FE_Element.h>
#include <DOF_Group.h>
#include <Vector.h>
#include <Matrix.h>
#include <Logging.h>

EigenSOE::EigenSOE(EigenSolver &theEigenSolver, int classTag)
  : MovableObject(classTag)
  , theSolver(&theEigenSolver)
{

}


EigenSOE::EigenSOE(int classTag)
  : MovableObject(classTag)
  , theSolver(nullptr)
{

}


EigenSOE::~EigenSOE()
{
  if (theSolver != nullptr)
    delete theSolver;
}

int 
EigenSOE::solve(int numModes, bool generalized, bool findSmallest)
{
  if (theSolver == nullptr)
    return -1;
  else
    return theSolver->solve(numModes, generalized, findSmallest);
}

int
EigenSOE::setSolver(EigenSolver &newSolver)
{
  theSolver = &newSolver;
  return 0;
}

EigenSolver *
EigenSOE::getSolver()
{	
  return theSolver;
}

const Vector &
EigenSOE::getEigenvector(int mode)
{
  return theSolver->getEigenvector(mode);
}

double 
EigenSOE::getEigenvalue(int mode)
{
  return theSolver->getEigenvalue(mode);
}


int 
EigenSOE::setLinks(AnalysisModel &)
{
  return 0;
}


int 
EigenSOE::formSystem(AnalysisModel &model, bool generalized)
{
  int result = 0;

  //
  // zero A and M
  //
  this->zeroA();
  this->zeroM();

  //
  // form K
  //
  FE_EleIter &theEles = model.getFEs();
  FE_Element *elePtr;
  while ((elePtr = theEles()) != nullptr) {
    elePtr->zeroTangent();
    elePtr->addKtToTang(1.0);
    if (this->addA(elePtr->getTangent(nullptr), elePtr->getID()) < 0) {
      opserr << G3_WARN_PROMPT << "eigen -";
      opserr << " failed in addA for ID " << elePtr->getID();
      result = -2;
    }
  }

  //
  // If generalized is true, form M
  //
  if (generalized == true) {
    FE_EleIter &theEles2 = model.getFEs();
    while ((elePtr = theEles2()) != nullptr) {
      elePtr->zeroTangent();
      elePtr->addMtoTang(1.0);
      if (this->addM(elePtr->getTangent(nullptr), elePtr->getID()) < 0) {
        opserr << "WARNING BasicAnalysisBuilder::eigen -";
        opserr << " failed in addA for ID " << elePtr->getID() << "\n";
        result = -2;
      }
    }

    DOF_Group *dofPtr;
    DOF_GrpIter &theDofs = model.getDOFs();
    while ((dofPtr = theDofs()) != nullptr) {
      dofPtr->zeroTangent();
      dofPtr->addMtoTang(1.0);
      if (this->addM(dofPtr->getTangent(0), dofPtr->getID()) < 0) {
        opserr << G3_WARN_PROMPT << "failed in addM for ID " << dofPtr->getID() << "\n";
        result = -3;
      }
    }
  }
  return 0;
}