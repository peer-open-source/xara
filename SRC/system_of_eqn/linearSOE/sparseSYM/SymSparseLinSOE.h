// File: ~/system_of_eqn/linearSOE/LawSolver/SymSparseLinSOE.h
//
// Written: Jun Peng  (junpeng@stanford.edu)
//          Prof. Kincho H. Law
//          Stanford University
// Created: 12/98
// Revision: A
//
// Description: This file contains the class definition for 
// SymSparseLinSOE.h. It stores the sparse matrix A in a fashion
// that only store the none zero entries.
//
// What: "@(#) SymSparseLinSOE.h, revA"
//
// Almost all the information (Matrix A and Vector B) is stored as 
// global variables in the file "symbolic.h".
//
#pragma once
#include <LinearSOE.h>
#include <Vector.h>

extern "C" {
  #include <FeStructs.h>
}

class SymSparseLinSolver;

class SymSparseLinSOE : public LinearSOE
{
  public:
    SymSparseLinSOE(SymSparseLinSolver &theSolver, int numbering);        
    SymSparseLinSOE(int N, int NNZ, int *rowStartA, int *colA,
		    SymSparseLinSolver &theSolver, 
        int numbering
    );

    ~SymSparseLinSOE();

    int getNumEqn() const;
    int setSize(Graph &theGraph);
    int addA(const Matrix &, const ID &, double fact = 1.0);
    int addB(const Vector &, const ID &, double fact = 1.0);    
    int setB(const Vector &, double fact = 1.0);        
    
    void zeroA();
    void zeroB();
    
    const Vector &getX();
    const Vector &getB();

    void setX(int loc, double value);        
    void setX(const Vector &x);  

    friend class SymSparseLinSolver;

  private:
    int size;            // order of A
    int nnz;             // number of non-zeros in A
    double *B, *X;       // 1d arrays containing coefficients of B and X
    int *colA, *rowStartA;  //These are (ADJNCY, XADJ) pair.

    Vector *vectX;
    Vector *vectB;
    int Bsize;
    bool factored;

    int      LSPARSE;
    int      nblks;
    int      *xblk,  *invp;
    double   *diag, **penv;
    int      *rowblks;
    OFFDBLK  **begblk;
    OFFDBLK  *first;

};
