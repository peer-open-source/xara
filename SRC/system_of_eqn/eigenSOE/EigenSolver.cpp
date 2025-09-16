//
// Written: Jun Peng
// Created: Sat Feb. 6, 1999
// Revision: A
//
// Description: This file contains the class definition of EigenSOE.
// EigenSOE is a subclass of Solver.
// This is an abstract base class and thus no objects of it's type
// can be instantiated. Instances of EigenSolver are used to solve
// a EigenSOE. (perform eigen analysis)
//

#include <EigenSolver.h>


EigenSolver::EigenSolver(int classTag)
  :MovableObject(classTag)
{

}


EigenSolver::~EigenSolver()
{

}

