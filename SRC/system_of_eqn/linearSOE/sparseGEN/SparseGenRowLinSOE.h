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

// $Revision: 1.2 $
// $Date: 2006-01-10 00:42:51 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/SparseGenRowLinSOE.h,v $

#ifndef SparseGenRowLinSOE_h
#define SparseGenRowLinSOE_h

// Written: fmk
// Created: 04/05
// Revision: A
//
// Description: This file contains the class definition for SparseGenRowLinSOE
// SparseGenRowLinSOE is a subclass of LinearSOE. It stores the matrix equation
// Ax=b using the sparse row-compacted storage scheme for storing the
// matrix A.
//
// What: "@(#) SparseGenRowLinSOE.h, revA"

#include <LinearSOE.h>
#include <Vector.h>

class SparseGenRowLinSolver;

class SparseGenRowLinSOE : public LinearSOE {
public:
  SparseGenRowLinSOE(SparseGenRowLinSolver& theSolver);

  ~SparseGenRowLinSOE();

    int getNumEqn() const;
    int setSize(Graph &theGraph);

    void zeroA();
    void zeroB();

    const Vector& getX();
    const Vector& getB();
    double normRHS();

    int addA(const Matrix&, const ID&, double fact = 1.0);
    int addB(const Vector&, const ID&, double fact = 1.0);
    int setB(const Vector&, double fact = 1.0);

    int sendSelf(int commitTag, Channel &);
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);

#if 0
    void setX(int loc, double value);        
    void setX(const Vector &x);
    int setSparseGenRowSolver(SparseGenRowLinSolver &newSolver);    
    friend class PetscSparseSeqSolver;    
    friend class CulaSparseSolverS4;
    friend class CulaSparseSolverS5;    
	  friend class CuSPSolver;
#endif

protected:
private:
  int size;  // order of A
  int nnz;   // number of non-zeros in A
  double* A; // 1d arrays containing coefficients of A, B and X
  Vector X, B;
  int *colA, *rowStartA; // int arrays containing info about coeficientss in A
  int Asize, Bsize;      // size of the 1d array holding A
  bool factored;
};


#endif
