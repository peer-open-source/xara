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
// Description: This file contains the implementation for BandGenLinSOE
// Written: fmk 
//
#include <stdlib.h>

#include <BandGenLinSOE.h>
#include <BandGenLinSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <math.h>
#include <iostream>
#include <assert.h>

BandGenLinSOE::BandGenLinSOE(BandGenLinSolver &theSolvr)
:LinearSOE(theSolvr, LinSOE_TAGS_BandGenLinSOE),
 theSolver(&theSolvr),
 size(0), numSuperD(0), numSubD(0), A(0), B(0), X(0), 
 vectX(0), vectB(0), Asize(0), Bsize(0), factored(false)
{
  theSolvr.setLinearSOE(*this);
}


BandGenLinSOE::BandGenLinSOE(int classTag)
:LinearSOE(classTag),
 size(0), numSuperD(0), numSubD(0), A(0), 
 B(0), X(0), vectX(0), vectB(0), 
 Asize(0), Bsize(0), 
 factored(false)
{

}


#if 0
BandGenLinSOE::BandGenLinSOE(int N, int numSuperDiag, int numSubDiag,
                             BandGenLinSolver &theSolvr)
:LinearSOE(theSolvr, LinSOE_TAGS_BandGenLinSOE),
 theSolver(&theSolvr),
 size(N), numSuperD(numSuperDiag), numSubD(numSubDiag), A(0), 
 B(0), X(0), vectX(0), vectB(0), 
 Asize(0), Bsize(0), 
 factored(false)
{
  Asize = N * (2*numSubD + numSuperD +1);
  A = new double[Asize];

  // zero the matrix
  for (int i=0; i<Asize; i++)
    A[i] = 0;

  B = new double[size];
  X = new double[size];

  Bsize = size;
  // zero the vectors
  for (int j=0; j<size; j++) {
    B[j] = 0;
    X[j] = 0;
  }

  vectX.setData(X,size);
  vectB.setData(B,size);

  theSolvr.setLinearSOE(*this);        
  
  int solverOK = theSolvr.setSize();
  if (solverOK < 0) {
      // opserr << "WARNING BandGenLinSOE::BandGenLinSOE :";
      // opserr << " solver failed setSize() in constructor\n";
  }
}
#endif

int
BandGenLinSOE::getNumEqn() const
{
  return size;
}


BandGenLinSOE::~BandGenLinSOE()
{
  if (A != nullptr) delete [] A;
  if (B != nullptr) delete [] B;
  if (X != nullptr) delete [] X; 
}



int 
BandGenLinSOE::setSize(Graph &theGraph)
{
  int result = 0;
  size = theGraph.getNumVertex();
  factored = false;

  /*
    * determine the number of superdiagonals and subdiagonals
    */

  numSubD = 0;
  numSuperD = 0;

  Vertex *vertexPtr;
  VertexIter &theVertices = theGraph.getVertices();

  while ((vertexPtr = theVertices()) != nullptr) {
    int vertexNum = vertexPtr->getTag();
    const ID &theAdjacency = vertexPtr->getAdjacency();
    for (int i=0; i<theAdjacency.Size(); i++) {
      int otherNum = theAdjacency(i);
      int diff = vertexNum - otherNum;
      if (diff > 0) {
        if (diff > numSuperD)
          numSuperD = diff;
      } else 
        if (diff < numSubD)
          numSubD = diff;
    }
  }
  numSubD *= -1;

  int newSize = size * (2*numSubD + numSuperD +1);
  if (newSize > Asize) { 
    // we have to get another space for A
    if (A != 0) 
      delete [] A;

    A = new double[newSize];
    Asize = newSize;
  }

  // zero the matrix
  for (int i=0; i<Asize; i++)
    A[i] = 0;

  if (size > Bsize) {
    // we have to get space for the vectors
    // delete the old
    if (B != nullptr) delete [] B;
    if (X != nullptr) delete [] X;

    // create the new
    B = new double[size];
    X = new double[size];

    Bsize = size;
  }

  // get new Vector objects if size has changes
  vectX.setData(X,size);
  vectB.setData(B,size);

  // zero the vectors
  vectX.Zero();
  vectB.Zero();

  // invoke setSize() on the Solver
  LinearSOESolver *theSolvr = this->getSolver();
  int solverOK = theSolvr->setSize();
  if (solverOK < 0) {
    return solverOK;
  }

  return result;
}


int 
BandGenLinSOE::addA(const Matrix &m, const ID &id, double fact)
{
  assert(id.Size() == m.noRows() && id.Size() == m.noCols());

  // check for a quick return 
  if (fact == 0.0)  
    return 0;

  factored = false;
  
  // check that m and id are of similar size
  const int idSize = id.Size();    

  int ldA = 2*numSubD + numSuperD + 1;


  if (fact == 1.0) { // do not need to multiply 
    for (int i=0; i<idSize; i++) {
      int col = id(i);
      if (col < size && col >= 0) {
        double *coliiPtr = A + col*ldA + numSubD + numSuperD;
        for (int j=0; j<idSize; j++) {
          int row = id(j);
          if (row < size && row >= 0) {
            int diff = col - row;
            if (diff > 0) {
              if (diff <= numSuperD) {
                double *APtr = coliiPtr - diff;
                *APtr += m(j,i);
              }
            }
            else {
              diff *= -1;
              if (diff <= numSubD) {
                double *APtr = coliiPtr + diff;
                *APtr += m(j,i);
              }
            }
          }
        }  // for j
      }
    }  // for i
  } else {
    for (int i=0; i<idSize; i++) {
      int col = id(i);
      if (col < size && col >= 0) {
        double *coliiPtr = A + col*ldA + numSubD + numSuperD;
        for (int j=0; j<idSize; j++) {
          int row = id(j);
          if (row <size && row >= 0) {                    
            int diff = col - row;
            if (diff > 0) {
              if (diff <= numSuperD) {
                double *APtr = coliiPtr - diff;
                *APtr += m(j,i) *fact;
              }
            } else {
              diff *= -1;
              if (diff <= numSubD) {
                double *APtr = coliiPtr + diff;
                *APtr += m(j,i) *fact;
              }
            }
          }
        }  // for j
      }
    }  // for i
  }

  return 0;
}


    
int 
BandGenLinSOE::addB(const Vector &v, const ID &id, double fact)
{
  assert(id.Size() == v.Size() );

  // check for a quick return 
  if (fact == 0.0)
    return 0;

  // check that m and id are of similar size
  const int idSize = id.Size();        

  if (fact == 1.0) { // do not need to multiply if fact == 1.0
    for (int i=0; i<idSize; i++) {
      int pos = id(i);
      if (pos <size && pos >= 0)
          B[pos] += v(i);
    }
  } else if (fact == -1.0) {
    for (int i=0; i<idSize; i++) {
      int pos = id(i);
      if (pos <size && pos >= 0)
        B[pos] -= v(i);
    }
  } else {
    for (int i=0; i<idSize; i++) {
      int pos = id(i);
      if (pos <size && pos >= 0)
        B[pos] += v(i) * fact;
    }
  }
  return 0;
}


int
BandGenLinSOE::setB(const Vector &v, double fact)
{
  assert(v.Size() == size);

  // check for a quick return 
  if (fact == 0.0)
    return 0;

  else if (fact == 1.0) { // do not need to multiply if fact == 1.0
    for (int i=0; i<size; i++) {
      B[i] = v(i);
    }

  } else if (fact == -1.0) {
    for (int i=0; i<size; i++) {
      B[i] = -v(i);
    }
  } else {
    for (int i=0; i<size; i++) {
      B[i] = v(i) * fact;
    }
  }        
  return 0;
}


void 
BandGenLinSOE::zeroA()
{
  double *Aptr = A;
  int theSize = Asize;
  for (int i=0; i<theSize; i++)
    *Aptr++ = 0;
  
  factored = false;
}
        
void 
BandGenLinSOE::zeroB()
{
  double *Bptr = B;
  for (int i=0; i<size; i++)
    *Bptr++ = 0;
}


const Vector &
BandGenLinSOE::getX()
{
  return vectX;
}


const Vector &
BandGenLinSOE::getB()
{
  return vectB;
}



void 
BandGenLinSOE::setX(int loc, double value)
{
  if (loc < size && loc >= 0)
    X[loc] = value;
}

void 
BandGenLinSOE::setX(const Vector &x)
{
  if (x.Size() == vectX.Size())
    vectX = x;
}

