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
// Revision: A
//
// Description: This file contains the implementation for ProfileSPDLinSOE
//
#include <assert.h>
#include <ProfileSPDLinSOE.h>
#include <ProfileSPDLinSolver.h>
#include <Matrix.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <math.h>
#include <Logging.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <iostream>
using std::nothrow;


ProfileSPDLinSOE::ProfileSPDLinSOE(ProfileSPDLinSolver &the_Solver)
:LinearSOE(the_Solver, LinSOE_TAGS_ProfileSPDLinSOE),
 size(0), profileSize(0), A(0), B(0), X(0), vectX(0), vectB(0),
 iDiagLoc(0), Asize(0), Bsize(0), isAfactored(false), isAcondensed(false),
 numInt(0) 
{
    the_Solver.setLinearSOE(*this);
}


ProfileSPDLinSOE::ProfileSPDLinSOE(int classTag)
:LinearSOE(classTag),
 size(0), profileSize(0), A(0), B(0), X(0), vectX(0), vectB(0),
 iDiagLoc(0), Asize(0), Bsize(0), isAfactored(false), isAcondensed(false),
 numInt(0) 
{

}


ProfileSPDLinSOE::ProfileSPDLinSOE(ProfileSPDLinSolver &the_Solver, int classTag)
:LinearSOE(the_Solver, classTag),
 size(0), profileSize(0), A(0), B(0), X(0), vectX(0), vectB(0),
 iDiagLoc(0), Asize(0), Bsize(0), isAfactored(false), isAcondensed(false),
 numInt(0) 
{
    the_Solver.setLinearSOE(*this);
}


ProfileSPDLinSOE::ProfileSPDLinSOE(int N, int *iLoc,
				   ProfileSPDLinSolver &the_Solver)
:LinearSOE(the_Solver, LinSOE_TAGS_ProfileSPDLinSOE),
 size(0), profileSize(0), A(0), B(0), X(0), vectX(0), vectB(0),
 iDiagLoc(0), Asize(0), Bsize(0), isAfactored(false), isAcondensed(false),
 numInt(0)
{
    size = N;
    profileSize = iLoc[N-1];
    
    A = new  double[iLoc[N-1]]{};
    // zero the matrix
    Asize = iLoc[N-1];

    B = new  double[size]{};
    X = new  double[size]{};
    iDiagLoc = new  int[size]{};

    Bsize = size;

    vectX = new Vector(X,size);
    vectB = new Vector(B,size);
    
    the_Solver.setLinearSOE(*this);        
    
    int solverOK = the_Solver.setSize();
    if (solverOK < 0) {
	opserr << "WARNING ProfileSPDLinSOE::ProfileSPDLinSOE :";
	opserr << " solver failed setSize() in constructor\n";
    }
}

    
ProfileSPDLinSOE::~ProfileSPDLinSOE()
{
    if (A != 0) delete [] A;
    if (B != 0) delete [] B;
    if (X != 0) delete [] X;
    if (iDiagLoc != 0) delete [] iDiagLoc;
    if (vectX != 0) delete vectX;    
    if (vectB != 0) delete vectB;    
}


int 
ProfileSPDLinSOE::getNumEqn(void) const
{
    return size;
}

int 
ProfileSPDLinSOE::setSize(Graph &theGraph)
{
    int oldSize = size;
    int result = 0;
    size = theGraph.getNumVertex();

    // check we have enough space in iDiagLoc and iLastCol
    // if not delete old and create new
    if (size > Bsize) { 
	if (iDiagLoc != 0) delete [] iDiagLoc;
	iDiagLoc = new int[size]{};

    } else {
      // zero out iDiagLoc 
      for (int i=0; i<size; i++) {
          iDiagLoc[i] = 0;
      }
    }

    // now we go through the vertices to find the height of each col and
    // width of each row from the connectivity information.
    
    Vertex *vertexPtr;
    VertexIter &theVertices = theGraph.getVertices();

    while ((vertexPtr = theVertices()) != 0) {
	int vertexNum = vertexPtr->getTag();
	const ID &theAdjacency = vertexPtr->getAdjacency();
	int iiDiagLoc = iDiagLoc[vertexNum];
	int *iiDiagLocPtr = &(iDiagLoc[vertexNum]);

	for (int i=0; i<theAdjacency.Size(); i++) {
	    int otherNum = theAdjacency(i);
	    int diff = vertexNum-otherNum;
	    if (diff > 0) {
		if (iiDiagLoc < diff) {
		    iiDiagLoc = diff;
		    *iiDiagLocPtr = diff;
		}
	    } 
	}
    }


    // now go through iDiagLoc, adding 1 for the diagonal element
    // and then adding previous entry to give current location.
    if (iDiagLoc != 0)
	iDiagLoc[0] = 1; // NOTE FORTRAN ARRAY LOCATION

    for (int j=1; j<size; j++)
	iDiagLoc[j] = iDiagLoc[j] + 1 + iDiagLoc[j-1];

    if (iDiagLoc != 0)       
    	profileSize = iDiagLoc[size-1];

    // check if we need more space to hold A
    // if so then go get it
    if (profileSize > Asize) {
	// delete old space
	if (A != 0)
	    delete [] A;

	// get new zeroed-out space
	A = new double[profileSize]{};

	Asize = profileSize;

    } else {
      // zero the matrix
      for (int k=0; k<profileSize; k++)
          A[k] = 0;
    }

    isAfactored = false;
    isAcondensed = false;    

    if (size > Bsize) { // we have to get another space for A

	// delete the old
	if (B != 0) delete [] B;
	if (X != 0) delete [] X;

	// allocate zeroed-out memory 
	B = new double[size]{};
	X = new double[size]{};

    } else {
      // zero the vectors
      for (int l=0; l<size; l++) {
          B[l] = 0;
          X[l] = 0;
      }
    }
    
    if (size != oldSize) {
	if (vectX != 0)
	    delete vectX;
	if (vectB != 0)
	    delete vectB;
	
	vectX = new Vector(X,size);
	vectB = new Vector(B,size);
	
	if (size > Bsize)
	    Bsize = size;
    }
    
    // invoke setSize() on the Solver
    LinearSOESolver *the_Solver = this->getSolver();
    int solverOK = the_Solver->setSize();
    if (solverOK < 0) {
	// opserr << "WARNING ProfileSPDLinSOE::setSize :";
	// opserr << " solver failed setSize()\n";
	return solverOK;
    }    

    return result;
}

int 
ProfileSPDLinSOE::addA(const Matrix &m, const ID &id, double fact)
{
    assert(id.Size() == m.noRows() && id.Size() == m.noCols());

    // check for a quick return 
    if (fact == 0.0)  return 0;
    
    // check that m and id are of similar size
    int idSize = id.Size();    

    if (fact == 1.0) { // do not need to multiply 
	for (int i=0; i<idSize; i++) {
	    int col = id(i);
	    if (col < size && col >= 0) {
		double *coliiPtr = &A[iDiagLoc[col] -1]; // -1 as fortran indexing 
		int minColRow;
		if (col == 0)
		    minColRow = 0;
		else
		    minColRow = col - (iDiagLoc[col] - iDiagLoc[col-1]) +1;
		for (int j=0; j<idSize; j++) {
		    int row = id(j);
		    if (row <size && row >= 0 && 
			row <= col && row >= minColRow) { 

			// we only add upper and inside profile
			double *APtr = coliiPtr + (row-col);
			 *APtr += m(j,i);
		     }
		}  // for j
	    } 
	}  // for i
    } else {
	for (int i=0; i<idSize; i++) {
	    int col = id(i);
	    if (col < size && col >= 0) {
		double *coliiPtr = &A[iDiagLoc[col] -1]; // -1 as fortran indexing 		
		int minColRow;
		if (col == 0)
		    minColRow = 0;
		else
		    minColRow = col - (iDiagLoc[col] - iDiagLoc[col-1]) +1;

		for (int j=0; j<idSize; j++) {
		    int row = id(j);
		    if (row < size && row >= 0 && 
			row <= col && row >= minColRow) { 

			// we only add upper and inside profile
			double *APtr = coliiPtr + (row-col);
			 *APtr += m(j,i) * fact;
		     }
		}  // for j
	    } 
	}  // for i
    
    }
    return 0;
}


int 
ProfileSPDLinSOE::addColA(const Vector &colData, int col, double fact)
{
  // check for a quick return 
  if (fact == 0.0)  return 0;
  
  if (colData.Size() != size) {
    opserr << "ProfileSPDLinSOE::addColA() - colData size not equal to n\n";
    return -1;
  }

  if (col > size && col < 0) {
    opserr << "ProfileSPDLinSOE::addColA() - col " << col << "outside range 0 to " << size << endln;
    return -1;
  }

  
  if (fact == 1.0) { // do not need to multiply 

    double *coliiPtr = &A[iDiagLoc[col] -1]; // -1 as fortran indexing 
    int minColRow;
    if (col == 0)
      minColRow = 0;
    else
      minColRow = col - (iDiagLoc[col] - iDiagLoc[col-1]) +1;
    for (int row=0; row<size; row++) {
      double data = colData(row);
      if (data != 0) {
	if (row < size && row >= 0 && 
	    row <= col && row >= minColRow) { 
	  
	  // we only add upper and inside profile
	  double *APtr = coliiPtr + (row-col);
	  *APtr += data;
	}
      }
    }

  } else {

    double *coliiPtr = &A[iDiagLoc[col] -1]; // -1 as fortran indexing 
    int minColRow;
    if (col == 0)
      minColRow = 0;
    else
      minColRow = col - (iDiagLoc[col] - iDiagLoc[col-1]) +1;

    for (int row=0; row<size; row++) {
      double data = colData(row);
      if (data != 0) {
	if (row <size && row >= 0 && 
	    row <= col && row >= minColRow) { 
	  
	  // we only add upper and inside profile
	  double *APtr = coliiPtr + (row-col);
	  *APtr += data * fact;
	}
      }
    }
  }
  return 0;
}

    
int 
ProfileSPDLinSOE::addB(const Vector &v, const ID &id, double fact)
{
    
    // check for a quick return 
    if (fact == 0.0)  return 0;

    // check that m and id are of similar size
    int idSize = id.Size();        
    if (idSize != v.Size() ) {
	opserr << "ProfileSPDLinSOE::addB() -";
	opserr << " Vector and ID not of similar sizes\n";
	return -1;
    }    
    
    if (fact == 1.0) { // do not need to multiply if fact == 1.0
	for (int i=0; i<id.Size(); i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] += v(i);
	}
    } else if (fact == -1.0) { // do not need to multiply if fact == -1.0
	for (int i=0; i<id.Size(); i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] -= v(i);
	}
    } else {
	for (int i=0; i<id.Size(); i++) {
	    int pos = id(i);
	    if (pos <size && pos >= 0)
		B[pos] += v(i) * fact;
	}
    }	
    return 0;
}


int
ProfileSPDLinSOE::setB(const Vector &v, double fact)
{
    // check for a quick return 
    if (fact == 0.0)  return 0;


    if (v.Size() != size) {
	opserr << "WARNING BandGenLinSOE::setB() -";
	opserr << " incompatible sizes " << size << " and " << v.Size() << endln;
	return -1;
    }
    
    if (fact == 1.0) { // do not need to multiply if fact == 1.0
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
ProfileSPDLinSOE::zeroA(void)
{
    double *Aptr = A;
    for (int i=0; i<Asize; i++)
	*Aptr++ = 0;
    
    isAfactored = false;
}
	
void 
ProfileSPDLinSOE::zeroB(void)
{
    double *Bptr = B;
    for (int i=0; i<size; i++)
	*Bptr++ = 0;
}


void 
ProfileSPDLinSOE::setX(int loc, double value)
{
    if (loc < size && loc >=0)
	X[loc] = value;
}

void 
ProfileSPDLinSOE::setX(const Vector &x)
{
  if (x.Size() == size && vectX != 0)
    *vectX = x;
}

const Vector &
ProfileSPDLinSOE::getX(void)
{
  assert(vectX != nullptr);
  return *vectX;
}

const Vector &
ProfileSPDLinSOE::getB(void)
{
  assert(vectB != nullptr);
  return *vectB;
}

double 
ProfileSPDLinSOE::normRHS(void)
{
    double norm =0.0;
    for (int i=0; i<size; i++) {
	double Yi = B[i];
	norm += Yi*Yi;
    }
    return sqrt(norm);
    
}    


int
ProfileSPDLinSOE::setProfileSPDSolver(ProfileSPDLinSolver &newSolver)
{
    newSolver.setLinearSOE(*this);

    if (size != 0) {
	int solverOK = newSolver.setSize();
	if (solverOK < 0) {
	    opserr << "WARNING:ProfileSPDLinSOE::setSolver :";
	    opserr << "the new solver could not setSeize() - staying with old\n";
	    return -1;
	}
    }
    
    return this->setSolver(newSolver);
}


int 
ProfileSPDLinSOE::sendSelf(int cTag, Channel &theChannel)
{
    return 0;
}


int 
ProfileSPDLinSOE::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return 0;
}
