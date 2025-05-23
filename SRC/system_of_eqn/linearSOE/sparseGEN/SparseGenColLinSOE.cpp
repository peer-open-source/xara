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

// $Revision: 1.8 $
// $Date: 2009-05-11 20:57:49 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/SparseGenColLinSOE.cpp,v $


// Written: fmk
// Created: 04/98
// Revision: A
//
// Description: This file contains the implementation for SparseGenColLinSOE
//
#include <assert.h>
#include <stdlib.h>
#include <SparseGenColLinSOE.h>
#include <SparseGenColLinSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <math.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <iostream>
using std::nothrow;

SparseGenColLinSOE::SparseGenColLinSOE(SparseGenColLinSolver& the_Solver)
 : LinearSOE(the_Solver, LinSOE_TAGS_SparseGenColLinSOE),
   size(0),
   nnz(0),
   A(0),
   B(0),
   X(0),
   rowA(0),
   colStartA(0),
   Asize(0),
   Bsize(0),
   factored(false)
{
  the_Solver.setLinearSOE(*this);
}


SparseGenColLinSOE::SparseGenColLinSOE()
 : LinearSOE(LinSOE_TAGS_SparseGenColLinSOE),
   size(0),
   nnz(0),
   A(0),
   B(0),
   X(0),
   rowA(0),
   colStartA(0),
   Asize(0),
   Bsize(0),
   factored(false)
{
}

SparseGenColLinSOE::SparseGenColLinSOE(int classTag)
 : LinearSOE(classTag),
   size(0),
   nnz(0),
   A(0),
   B(0),
   X(0),
   rowA(0),
   colStartA(0),
   Asize(0),
   Bsize(0),
   factored(false)
{
}


SparseGenColLinSOE::SparseGenColLinSOE(SparseGenColLinSolver& the_Solver, int classTag)
 : LinearSOE(the_Solver, classTag),
   size(0),
   nnz(0),
   A(0),
   B(0),
   X(0),
   rowA(0),
   colStartA(0),
   Asize(0),
   Bsize(0),
   factored(false)
{
  //    the_Solver.setLinearSOE(*this);
}


SparseGenColLinSOE::SparseGenColLinSOE(int N, int NNZ, int* ColStartA, int* RowA,
                                       SparseGenColLinSolver& the_Solver)
 : LinearSOE(the_Solver, LinSOE_TAGS_SparseGenColLinSOE),
   size(N),
   nnz(NNZ),
   A(0),
   B(N),
   X(N),
   rowA(RowA),
   colStartA(ColStartA),
   Asize(0),
   Bsize(0),
   factored(false)
{

  A = new double[NNZ]{};

  Asize = NNZ;
  Bsize = size;


  the_Solver.setLinearSOE(*this);

  // invoke setSize() on the Solver
  if (the_Solver.setSize() < 0) {
    // opserr << "WARNING :SparseGenColLinSOE::SparseGenColLinSOE :";
    // opserr << " solver failed setSize() in constructor\n";
  }
}


SparseGenColLinSOE::~SparseGenColLinSOE()
{
  if (A != 0)
    delete[] A;
  if (colStartA != 0)
    delete[] colStartA;
  if (rowA != 0)
    delete[] rowA;
}


int
SparseGenColLinSOE::getNumEqn(void) const
{
  return size;
}

int
SparseGenColLinSOE::setSize(Graph& theGraph)
{

  int result  = 0;
  int oldSize = size;
  size        = theGraph.getNumVertex();

  // fist itearte through the vertices of the graph to get nnz
  Vertex* theVertex;
  int newNNZ              = 0;
  VertexIter& theVertices = theGraph.getVertices();
  while ((theVertex = theVertices()) != 0) {
    const ID& theAdjacency = theVertex->getAdjacency();
    newNNZ += theAdjacency.Size() + 1; // the +1 is for the diag entry
  }
  nnz = newNNZ;

  if (newNNZ > Asize) { // we have to get more space for A and rowA
    if (A != 0)
      delete[] A;
    if (rowA != 0)
      delete[] rowA;

    A    = new double[newNNZ]{};
    rowA = new int[newNNZ];

    Asize = newNNZ;
  }

  factored = false;
  B.resize(size);
  B.Zero();
  X.resize(size);
  X.Zero();

  if (size > Bsize) { // we have to get space for the vectors

    if (colStartA != 0)
      delete[] colStartA;
    colStartA = new int[size + 1];

    Bsize = size;
  }

  // fill in colStartA and rowA
  if (size != 0) {
    colStartA[0] = 0;
    int startLoc = 0;
    int lastLoc  = 0;
    for (int a = 0; a < size; a++) {

      theVertex = theGraph.getVertexPtr(a);
      if (theVertex == 0) {
        // opserr << "WARNING:SparseGenColLinSOE::setSize :";
        // opserr << " vertex " << a << " not in graph! - size set to 0\n";
        size = 0;
        return -1;
      }

      rowA[lastLoc++]        = theVertex->getTag(); // place diag in first
      const ID& theAdjacency = theVertex->getAdjacency();
      int idSize             = theAdjacency.Size();

      // now we have to place the entries in the ID into order in rowA
      for (int i = 0; i < idSize; i++) {

        int row         = theAdjacency(i);
        bool foundPlace = false;
        // find a place in rowA for current col
        for (int j = startLoc; j < lastLoc; j++)
          if (rowA[j] > row) {
            // move the entries already there one further on
            // and place col in current location
            for (int k = lastLoc; k > j; k--)

              rowA[k] = rowA[k - 1];
            rowA[j]    = row;
            foundPlace = true;
            j          = lastLoc;
          }
        if (foundPlace == false) // put in at the end
          rowA[lastLoc] = row;

        lastLoc++;
      }
      colStartA[a + 1] = lastLoc;
      ;
      startLoc = lastLoc;
    }
  }


  // invoke setSize() on the Solver
  LinearSOESolver* the_Solver = this->getSolver();
  int solverOK                = the_Solver->setSize();
  if (solverOK < 0) {
    return solverOK;
  }

  return result;
}

int
SparseGenColLinSOE::addA(const Matrix& m, const ID& id, double fact)
{

  // check that m and id are of similar size
  assert(id.Size() == m.noRows() && id.Size() == m.noCols());

  // check for a quick return
  if (fact == 0.0)
    return 0;

  int idSize = id.Size();

  if (fact == 1.0) { // do not need to multiply
    for (int i = 0; i < idSize; i++) {
      int col = id(i);
      if (col < size && col >= 0) {
        int startColLoc = colStartA[col];
        int endColLoc   = colStartA[col + 1];
        for (int j = 0; j < idSize; j++) {
          int row = id(j);
          if (row < size && row >= 0) {
            // find place in A using rowA
            for (int k = startColLoc; k < endColLoc; k++)
              if (rowA[k] == row) {
                A[k] += m(j, i);
                k = endColLoc;
              }
          }
        } // for j
      }
    } // for i
  } else {
    for (int i = 0; i < idSize; i++) {
      int col = id(i);
      if (col < size && col >= 0) {
        int startColLoc = colStartA[col];
        int endColLoc   = colStartA[col + 1];
        for (int j = 0; j < idSize; j++) {
          int row = id(j);
          if (row < size && row >= 0) {
            // find place in A using rowA
            for (int k = startColLoc; k < endColLoc; k++)
              if (rowA[k] == row) {
                A[k] += fact * m(j, i);
                k = endColLoc;
              }
          }
        } // for j
      }
    } // for i
  }
  return 0;
}


int
SparseGenColLinSOE::addB(const Vector& v, const ID& id, double fact)
{
  // check that m and id are of similar size
  assert(id.Size() == v.Size());

  // check for a quick return
  if (fact == 0.0)
    return 0;

  int idSize = id.Size();

  double* b = &B[0];
  if (fact == 1.0) { // do not need to multiply if fact == 1.0
    for (int i = 0; i < idSize; i++) {
      int pos = id(i);
      if (pos < size && pos >= 0)
        b[pos] += v(i);
    }
  } else if (fact == -1.0) { // do not need to multiply if fact == -1.0
    for (int i = 0; i < idSize; i++) {
      int pos = id(i);
      if (pos < size && pos >= 0)
        b[pos] -= v(i);
    }
  } else {
    for (int i = 0; i < idSize; i++) {
      int pos = id(i);
      if (pos < size && pos >= 0)
        b[pos] += v(i) * fact;
    }
  }

  return 0;
}


int
SparseGenColLinSOE::setB(const Vector& v, double fact)
{
  assert(v.Size() == size);

  // check for a quick return
  if (fact == 0.0)
    return 0;

  double* b = &B[0];
  if (fact == 1.0) { // do not need to multiply if fact == 1.0
    for (int i = 0; i < size; i++) {
      b[i] = v(i);
    }
  } else if (fact == -1.0) {
    for (int i = 0; i < size; i++) {
      b[i] = -v(i);
    }
  } else {
    for (int i = 0; i < size; i++) {
      b[i] = v(i) * fact;
    }
  }
  return 0;
}

void
SparseGenColLinSOE::zeroA()
{
  double* Aptr = A;
  for (int i = 0; i < Asize; i++)
    *Aptr++ = 0;

  factored = false;
}

void
SparseGenColLinSOE::zeroB()
{
  B.Zero();
}

void
SparseGenColLinSOE::setX(int loc, double value)
{
  if (loc < size && loc >= 0)
    X[loc] = value;
}

void
SparseGenColLinSOE::setX(const Vector& x)
{
  if (x.Size() == size)
    X = x;
}

const Vector&
SparseGenColLinSOE::getX()
{
  return X;
}

const Vector&
SparseGenColLinSOE::getB()
{
  return B;
}

double
SparseGenColLinSOE::normRHS()
{
  double norm = 0.0;
  for (int i = 0; i < size; i++) {
    double Yi = B[i];
    norm += Yi * Yi;
  }
  return sqrt(norm);
}


int
SparseGenColLinSOE::setSparseGenColSolver(SparseGenColLinSolver& newSolver)
{
  newSolver.setLinearSOE(*this);

  if (size != 0) {
    int solverOK = newSolver.setSize();
    if (solverOK < 0) {
      // opserr << "WARNING:SparseGenColLinSOE::setSolver :";
      // opserr << "the new solver could not setSeize() - staying with old\n";
      return -1;
    }
  }

  return this->LinearSOE::setSolver(newSolver);
}


int
SparseGenColLinSOE::sendSelf(int cTag, Channel& theChannel)
{
  return 0;
}

int
SparseGenColLinSOE::recvSelf(int cTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
  return 0;
}
