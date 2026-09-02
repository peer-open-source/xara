// File: ~/system_of_eqn/linearSOE/LawSolver/SymSparseLinSolver.h
//
// Written: Jun Peng  (junpeng@stanford.edu)
//          Prof. Kincho H. Law
//          Stanford University
// Created: 12/98
// Revision: A
//
// Description: This file contains the class definition for 
// SymSparseinSolver. It solves the SymSparseLinSOEobject by calling
// some "C" functions. The solver used here is generalized sparse
// solver. The user can choose three different ordering schema.
//
// What: "@(#) SymSparseLinSolver.h, revA"
#pragma once

#include <LinearSOESolver.h>


class SymSparseLinSOE;

class SymSparseLinSolver : public LinearSOESolver
{
  public:
    SymSparseLinSolver();     
    ~SymSparseLinSolver();

    int solve() override;
    int solve(const Vector& B, Vector& X) override;
    int setSize() override;

    int setLinearSOE(SymSparseLinSOE &theSOE); 

  private:

    SymSparseLinSOE *theSOE;
    
};
