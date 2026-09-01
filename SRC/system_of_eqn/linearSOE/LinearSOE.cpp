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
// Written: fmk 
// Created: 11/96
//
// Description: This file contains the implementation of LinearSOE.
//
// What: "@(#) LinearSOE.C, revA"

#include <LinearSOE.h>
#include <LinearSOESolver.h>
#include <OPS_Stream.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

LinearSOE::LinearSOE(LinearSOESolver &theLinearSOESolver, int classtag)
  : MovableObject(classtag)
  , theSolver(&theLinearSOESolver)
  , m_fwd_update(nullptr)
  , m_inv_update(nullptr)
{

}

LinearSOE::LinearSOE(int classtag)
:MovableObject(classtag)
, theSolver(0)
, m_fwd_update(nullptr)
, m_inv_update(nullptr)
{

}


LinearSOE::~LinearSOE()
{
  if (theSolver != nullptr)
    delete theSolver;
}

int 
LinearSOE::solve()
{
  if (theSolver != 0)
    return (theSolver->solve());
  else 
    return -1;
}

#if 0
int
LinearSOE::reset()
{
  m_inv_update = 0;
  m_fwd_update = 0;
  return 0;
}
#endif


int
LinearSOE::setInverseUpdate(LinearAction * update)
{
  int ok = 0;
  if (m_inv_update != nullptr && update != nullptr)
    ok = -1;

  m_inv_update = update;
  return ok;
}


int
LinearSOE::setForwardUpdate(LinearAction * update)
{
  int ok = 0;
  if (m_fwd_update != nullptr && update != nullptr)
    ok = -1;

  m_fwd_update = update;
  return ok;
}


int
LinearSOE::solve(const Vector& b, Vector& x)
{
  LinearSOESolver* solver = this->getSolver();
  assert(solver != nullptr);

  x.Zero();
  int ok = solver->solve(b, x);
  if (ok < 0)
    return ok;

  if (m_fwd_update)
   if ((ok = m_fwd_update->solve(b,x)) < 0)
     return ok;

  if (m_inv_update)
    if ((ok = m_inv_update->apply(b,x)) < 0)
      return ok;

  return 0;
}


int
LinearSOE::formAp(const Vector &p, Vector &Ap)
{
  return -1;
}


double
LinearSOE::getDeterminant()
{
  if (theSolver != nullptr)
    return theSolver->getDeterminant();
  else 
    return 0;
}


int 
LinearSOE::setSolver(LinearSOESolver &newSolver)
{
  theSolver = &newSolver;
  return 0;
}


LinearSOESolver *
LinearSOE::getSolver()
{
  return theSolver;
}


int
LinearSOE::addA(const Matrix &)
{
  return -1;
}


int
LinearSOE::saveSparseA(OPS_Stream& output, int baseIndex)
{
  const Matrix* A = this->getA();
  
  if (A == nullptr) {
    return -1;
  }
  
  int rows = A->noRows();
  int cols = A->noCols();
  int nnz = rows * cols;
  
  // Assume the header is already written to output stream
  
  output << rows << " " << cols << " " << nnz << "\n";
  
  // Write all elements with base index
  int nnz_written = 0;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      double val = (*A)(i,j);
      output << i + baseIndex << " " << j + baseIndex << " " << val << "\n";
      nnz_written++;
    }
  }
  if (nnz_written != nnz) {
    return -1;
  }
  return 0;
}

int
LinearSOE::getSparseA(ID& rowIndices, ID& colIndices, Vector& values, int baseIndex) 
{
  const Matrix* A = this->getA();
  
  if (A == nullptr) {
    return -1;
  }
  
  int rows = A->noRows();
  int cols = A->noCols();
  int nnz = rows * cols;
  
  // Resize vectors to hold all elements
  rowIndices.resize(nnz);
  colIndices.resize(nnz);
  values.resize(nnz);
  
  // Fill vectors with all elements
  int idx = 0;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      rowIndices(idx) = i + baseIndex;
      colIndices(idx) = j + baseIndex;
      values(idx) = (*A)(i,j);
      idx++;
    }
  }
  
  return 0;
}

int
LinearSOE::getSparseA(std::vector<int>& rowIndices, std::vector<int>& colIndices, std::vector<double>& values, int baseIndex) {
  const Matrix* A = this->getA();
  
  if (A == nullptr) {
    return -1;
  }

  int rows = A->noRows();
  int cols = A->noCols();
  int nnz = rows * cols;
  
  // Resize vectors to hold all elements
  rowIndices.resize(nnz);
  colIndices.resize(nnz);
  values.resize(nnz);
  
  // Fill vectors with all elements
  int idx = 0;
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      rowIndices[idx] = i + baseIndex;
      colIndices[idx] = j + baseIndex;
      values[idx] = (*A)(i,j);
      idx++;
    }
  }
  
  return 0;
}