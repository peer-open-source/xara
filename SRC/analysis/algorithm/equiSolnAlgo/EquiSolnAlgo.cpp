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
// File: ~/OOP/analysis/algorithm/EquiSolnAlgo.C 
// 
// Written: fmk 
// Created: 11/96 
// Revision: A 
//
//
// Description: This file contains the class implementation for 
// EquiSolnAlgo. EquiSolnAlgo is an abstract base class, 
// i.e. no objects of it's type can be created.  Its subclasses deifine
// the sequence of operations to be performed in the analysis by static
// equilibrium of a finite element model.  
//
#include <EquiSolnAlgo.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>
#include <ConvergenceTest.h>

EquiSolnAlgo::EquiSolnAlgo(int clasTag)
:SolutionAlgorithm(clasTag),
 theIntegrator(nullptr), 
 theSysOfEqn(nullptr), 
 theTest(nullptr)
{

}

EquiSolnAlgo::~EquiSolnAlgo()
{

}

void 
EquiSolnAlgo::setLinks(IncrementalIntegrator &theNewIntegrator,
                       LinearSOE &theSOE,
                       ConvergenceTest *theConvergenceTest)
{
  theIntegrator = &theNewIntegrator;
  theSysOfEqn = &theSOE;
  theTest = theConvergenceTest;

  this->setConvergenceTest(theConvergenceTest);
}


int 
EquiSolnAlgo::setConvergenceTest(ConvergenceTest *theConvergenceTest)
{
  theTest = theConvergenceTest;
  return 0;
}

ConvergenceTest *
EquiSolnAlgo::getConvergenceTest()
{
  return theTest;
}


IncrementalIntegrator *
EquiSolnAlgo::getIncrementalIntegratorPtr() const
{
  return theIntegrator;
}



LinearSOE *
EquiSolnAlgo::getLinearSOEptr() const
{
  return theSysOfEqn;
}
