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
// Created: Mar 1998
// Revision: A
//
// Description: This file contains the class definition for 
// ProfileSPDLinDirectThreadSolver. ProfileSPDLinDirectThreadSolver will solve
// a linear system of equations stored using the profile scheme using threads.
// It solves a ProfileSPDLinSOE object using the LDL^t factorization and a block approach.
//
#include <ProfileSPDLinDirectThreadSolver.h>
#include <ProfileSPDLinSOE.h>
#include <math.h>
#include <assert.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Logging.h>
#include <pthread.h>
#include <threads/solaris2posix.h>
// #include <thread.h>

// global data that will be needed by the threads
struct ProfileTCB {
  mutex_t start_mutex;
  mutex_t end_mutex;
  mutex_t startBlock_mutex;
  mutex_t endBlock_mutex;

  cond_t start_cond;
  cond_t end_cond;
  cond_t startBlock_cond;
  cond_t endBlock_cond;

  double *A,*X,*B;
  int *iDiagLoc;
  int size;
  int blockSize;
  int  maxColHeight;
  double minDiagTol;
  int *RowTop;
  double **topRowPtr, *invD;

  int currentBlock;
  int workToDo;
  int info;
  int threadID;
  int numThreads;
  int threadsDone;
  int threadsDoneBlock;
} ; //TCB_ProfileSPDDirectThreadSolver;


void *ProfileSPDLinDirectThreadSolver_Worker(void *arg);


ProfileSPDLinDirectThreadSolver::ProfileSPDLinDirectThreadSolver
         (int numProcessors, int blckSize, double tol) 
:ProfileSPDLinSolver(SOLVER_TAGS_ProfileSPDLinDirectThreadSolver),
 NP(numProcessors), running(false),
 minDiagTol(tol), blockSize(blckSize), maxColHeight(0), 
 size(0), RowTop(0), topRowPtr(0), invD(0), tcb(new ProfileTCB)
{
    mutex_init(&tcb->start_mutex, USYNC_THREAD, 0);
    mutex_init(&tcb->end_mutex, USYNC_THREAD, 0);
    mutex_init(&tcb->startBlock_mutex, USYNC_THREAD, 0);
    mutex_init(&tcb->endBlock_mutex, USYNC_THREAD, 0);
    cond_init(&tcb->start_cond, USYNC_THREAD, 0);
    cond_init(&tcb->end_cond, USYNC_THREAD, 0);
    cond_init(&tcb->startBlock_cond, USYNC_THREAD, 0);
    cond_init(&tcb->endBlock_cond, USYNC_THREAD, 0);
    tcb->workToDo = 0;

#if 1
    // pthread_attr_t* attr = nullptr;
    pthread_attr_t  attr;
    pthread_attr_init(&attr);
    // pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
    // PTHREAD_CREATE_DETACHED
    pthread_t thread;
    for (int j = 0; j < NP; j++) {
      pthread_create(&thread, &attr, 
                 ProfileSPDLinDirectThreadSolver_Worker,
                 (void*)tcb);
      pthread_detach(thread);
    }
    pthread_attr_destroy(&attr);
#else 
    // THR_DAEMON implies THR_DETACHED
    int attr = THR_BOUND|THR_DAEMON;
    // Create the threads - Bound daemon threads
    for (int j = 0; j < NP; j++) 
      thr_create(NULL,0, ProfileSPDLinDirectThreadSolver_Worker, 
                 (void*)&TCB_ProfileSPDDirectThreadSolver, attr, NULL);
#endif 

}

    
ProfileSPDLinDirectThreadSolver::~ProfileSPDLinDirectThreadSolver()
{
    if (RowTop != 0) delete [] RowTop;
    if (topRowPtr != 0) free((void *)topRowPtr);
    if (invD != 0) delete [] invD;
    // TODO: free TCB
}

int
ProfileSPDLinDirectThreadSolver::setSize(void)
{
    assert(theSOE != nullptr);

    // check for quick return 
    if (theSOE->size == 0)
	return 0;

    if (size != theSOE->size) {    
      size = theSOE->size;

      if (RowTop != 0) delete [] RowTop;
      if (topRowPtr != 0) delete [] topRowPtr;
      if (invD != 0) delete [] invD;

      RowTop = new int[size];

      // we cannot use topRowPtr = new (double *)[size] with the cxx compiler
      topRowPtr = (double **)malloc(size *sizeof(double *));

      invD = new double[size]; 
    }


    // set some pointers
    double *A = theSOE->A;
    int *iDiagLoc = theSOE->iDiagLoc;

    // set RowTop and topRowPtr info

    maxColHeight = 0;
    RowTop[0] = 0;
    topRowPtr[0] = A;
    for (int j=1; j<size; j++) {
	int icolsz = iDiagLoc[j] - iDiagLoc[j-1];
        if (icolsz > maxColHeight) maxColHeight = icolsz;
	RowTop[j] = j - icolsz +  1;
	topRowPtr[j] = &A[iDiagLoc[j-1]]; // FORTRAN array indexing in iDiagLoc
    }

    size = theSOE->size;
    return 0;
}


int 
ProfileSPDLinDirectThreadSolver::solve(void)
{
    assert(theSOE != nullptr);
    
    // check for quick returns
    if (theSOE->size == 0)
	return 0;

    // set some pointers
    double *A = theSOE->A;
    double *B = &theSOE->B[0];
    double *X = &theSOE->X[0];
    int *iDiagLoc = theSOE->iDiagLoc;
    int size = theSOE->size;

    // copy B into X
    for (int ii=0; ii<size; ii++)
	X[ii] = B[ii];
    
    if (theSOE->isAfactored == false)  {

      mutex_lock(&tcb->start_mutex); 

      invD[0] = 1.0/A[0];	

      // assign the global variables the threads will need
      tcb->A = A;
      tcb->X = X;
      tcb->B = B;
      tcb->size = size;
      tcb->minDiagTol = minDiagTol;
      tcb->maxColHeight = maxColHeight;
      tcb->RowTop = RowTop;
      tcb->topRowPtr = topRowPtr;
      tcb->invD = invD;
      tcb->blockSize    = blockSize;
      tcb->iDiagLoc     = iDiagLoc;
      
      tcb->info = 0;
      tcb->workToDo     = NP;
      tcb->currentBlock = -1;
      tcb->threadID     = 0;
      tcb->numThreads   = NP;
      tcb->threadsDone  = 0;
      tcb->threadsDoneBlock = NP;

      // wake up the threads
      cond_broadcast(&tcb->start_cond);
      mutex_unlock(&tcb->start_mutex); 

      // yield the LWP
      thr_yield();
 
      // wait for all the threads to finish
      mutex_lock(&tcb->end_mutex);

      while (tcb->threadsDone < NP) 
        cond_wait(&tcb->end_cond, &tcb->end_mutex);

      mutex_unlock(&tcb->end_mutex);

      theSOE->isAfactored = true;
	
      // divide by diag term 
      double *bjPtr = X; 
      double *aiiPtr = invD;
      for (int j=0; j<size; j++) 
	*bjPtr++ = *aiiPtr++ * X[j];

      // now do the back substitution storing result in X
      for (int k=(size-1); k>0; k--) {
      
	int rowktop = RowTop[k];
	double bk = X[k];
	double *ajiPtr = topRowPtr[k]; 		
      
	for (int j=rowktop; j<k; j++) 
	  X[j] -= *ajiPtr++ * bk;
      }   	 
    } else { // just do forward and back substitution
      
      // do forward substitution 
      for (int i=1; i<size; i++) {
	    
	int rowitop = RowTop[i];	    
	double *ajiPtr = topRowPtr[i];
	double *bjPtr  = &X[rowitop];  
	double tmp = 0;	    
	    
	for (int j=rowitop; j<i; j++) 
	  tmp -= *ajiPtr++ * *bjPtr++; 
	    
	X[i] += tmp;
      }

      // divide by diag term 
      double *bjPtr = X; 
      double *aiiPtr = invD;
      for (int j=0; j<size; j++) 
	*bjPtr++ = *aiiPtr++ * X[j];

    
      // now do the back substitution storing result in X
      for (int k=(size-1); k>0; k--) {
      
	int rowktop = RowTop[k];
	double bk = X[k];
	double *ajiPtr = topRowPtr[k]; 		

	for (int j=rowktop; j<k; j++) 
	  X[j] -= *ajiPtr++ * bk;
      }   	 
    }    
    return 0;
}

int 
ProfileSPDLinDirectThreadSolver::setProfileSOE(ProfileSPDLinSOE &theNewSOE)
{
    assert(theSOE == nullptr); 
    theSOE = &theNewSOE;
    return 0;
}
	
int
ProfileSPDLinDirectThreadSolver::sendSelf(int cTag,
					  Channel &theChannel)
{
//     if (size != 0)
// 	opserr << "ProfileSPDLinDirectThreadSolver::sendSelf - does not send itself YET\n"; 
    return 0;
}


int 
ProfileSPDLinDirectThreadSolver::recvSelf(int cTag,
					  Channel &theChannel, 
					  FEM_ObjectBroker &theBroker)
{
    return 0;
}



void *ProfileSPDLinDirectThreadSolver_Worker(void *arg)
{
  struct ProfileTCB* tcb = (struct ProfileTCB*)arg;
  // Do this loop forever - or until all the Non-Daemon threads have exited
  while (true) {

      // Wait for some work to do
      mutex_lock(&tcb->start_mutex);

      while (tcb->workToDo == 0)
        cond_wait(&tcb->start_cond, 
                  &tcb->start_mutex);

      tcb->workToDo--;

      // get an id;
      int myID = tcb->threadID++;

      mutex_unlock(&tcb->start_mutex);

      // do some initialisation for the factorization
//    double *A     = tcb->A;
//    double *B     = tcb->B;
      double *X     = tcb->X;
      int *iDiagLoc = tcb->iDiagLoc;
      int size      = tcb->size;
      double minDiagTol = tcb->minDiagTol;
      int maxColHeight  = tcb->maxColHeight;
      int *RowTop = tcb->RowTop;
      double **topRowPtr = tcb->topRowPtr;
      double *invD = tcb->invD;
      int blockSize = tcb->blockSize;

      int currentBlock = 0;
      int numThreads = tcb->numThreads;

      // FACTOR 
      int startRow        = 0;
      int lastRow         = startRow + blockSize-1;
      int lastColEffected = lastRow  + maxColHeight -1;
      int nBlck = size/blockSize;
      if ((size % blockSize) != 0)
	nBlck++;

      // for every block across      
      for (int i=0; i<nBlck; i++) {

        int currentDiagThread = i%numThreads;

        if (myID == currentDiagThread) {

	  // first factor the diagonal block int Ui,i and Di
          int j;
	  for (j=0; j<blockSize; j++) {
            int currentRow = startRow + j;

	    if (currentRow < size) { // this is for case when size%blockSize != 0

	      int rowjTop = RowTop[currentRow];
	      //  double *akjPtr = &A[iDiagLoc[currentRow-1]];
	      double *akjPtr = topRowPtr[currentRow];
	      int maxRowijTop;
	      if (rowjTop < startRow) {
		akjPtr += startRow-rowjTop; // pointer to start of block row
		maxRowijTop = startRow;
	      } else
		maxRowijTop = rowjTop;

	      int k;
	      for (k=maxRowijTop; k<currentRow; k++) {
		double tmp = *akjPtr;
		int rowkTop = RowTop[k];
		int maxRowkjTop;
		double *alkPtr, *aljPtr;
		if (rowkTop < rowjTop) {
		  alkPtr = topRowPtr[k] + (rowjTop - rowkTop);
		  aljPtr = topRowPtr[currentRow];
		  maxRowkjTop = rowjTop;
		} else {
		  alkPtr = topRowPtr[k];
		  aljPtr = topRowPtr[currentRow] + (rowkTop - rowjTop);
		  maxRowkjTop = rowkTop;
		}

		for (int l = maxRowkjTop; l<k; l++) 
		  tmp -= *alkPtr++ * *aljPtr++;
		
		*akjPtr++ = tmp;
	      }

	      double ajj = *akjPtr;
	      akjPtr = topRowPtr[currentRow];
	      //	      akjPtr = &A[iDiagLoc[currentRow-1]];
	      double *bjPtr  = &X[rowjTop];  
	      double tmp = 0;	    

	      for (k=rowjTop; k<currentRow; k++){
		double akj = *akjPtr;
		double lkj = akj * invD[k];
		tmp -= lkj * *bjPtr++; 		
		*akjPtr++ = lkj;
		ajj = ajj -lkj * akj;
	      }

	      X[currentRow] += tmp;

	      // check that the diag > the tolerance specified
	      if (ajj <= 0.0) {
		// opserr << "ProfileSPDLinDirectThreadSolver::solve() - ";
		// opserr << " aii < 0 (i, aii): (" << currentRow << ", " << ajj << ")\n"; 
		tcb->info = -2;
		return NULL;
	      }
	      if (ajj <= minDiagTol) {
		// opserr << "ProfileSPDLinDirectThreadSolver::solve() - ";
		// opserr << " aii < minDiagTol (i, aii): (" << currentRow;
		// opserr << ", " << ajj << ")\n"; 
		tcb->info = -2;
		return NULL;
	      }		
	      invD[currentRow] = 1.0/ajj; 

	    } else 
	      j = blockSize;

	  }

          // allow other threads to now proceed
          mutex_lock(&tcb->startBlock_mutex);
	  tcb->currentBlock = i;

          cond_broadcast(&tcb->startBlock_cond);
          mutex_unlock(&tcb->startBlock_mutex);        


	  // now do rest of i'th block row belonging to thread
	  // doing a block of columns at a time forming Ui,j*Di
	  int currentCol = startRow + blockSize;
	  for (j=i+1; j<nBlck; j++) {

	    if (j%numThreads == myID) {
	      for (int k=0; k<blockSize; k++) {

		if (currentCol < size) { // this is for case when size%blockSize != 0

		  int rowkTop = RowTop[currentCol];
		  double *alkPtr = topRowPtr[currentCol];
		  // double *alkPtr = &A[iDiagLoc[currentCol-1]];
		  int maxRowikTop;
		  if (rowkTop < startRow) {
		    alkPtr += startRow-rowkTop; // pointer to start of block row
		    maxRowikTop = startRow;
		  } else
		    maxRowikTop = rowkTop;

		  for (int l=maxRowikTop; l<=lastRow; l++) {
		    double tmp = *alkPtr;
		    int rowlTop = RowTop[l];
		    int maxRowklTop;
		    double *amlPtr, *amkPtr;
		    if (rowlTop < rowkTop) {
		      amlPtr = topRowPtr[l] + (rowkTop - rowlTop);
		      amkPtr = topRowPtr[currentCol];
		      maxRowklTop = rowkTop;
		    } else {
		      amlPtr = topRowPtr[l];
		      amkPtr = topRowPtr[currentCol] + (rowlTop - rowkTop);
		      maxRowklTop = rowlTop;
		    }
		  
		    for (int m = maxRowklTop; m<l; m++) 
		      tmp -= *amkPtr++ * *amlPtr++;
		  
		    *alkPtr++ = tmp;
		  }
		  currentCol++;
		  if (currentCol > lastColEffected) {
		    k = blockSize;
		    j = nBlck;
		  }
		  
		} else
		  k = blockSize;
	      }
	    } else
	      currentCol += blockSize;
	  }
	} else {
          // wait till diag i is done 
          mutex_lock(&tcb->startBlock_mutex);
          while (tcb->currentBlock < i)
             cond_wait(&tcb->startBlock_cond, &tcb->startBlock_mutex);

          mutex_unlock(&tcb->startBlock_mutex);

          // FACTOR REST OF BLOCK ROW BELONGING TO THREAD
	  // now do rest of i'th block row belonging to thread
	  // doing a block of columns at a time forming Ui,j*Di
	  int currentCol = startRow + blockSize;
	  for (int j=i+1; j<nBlck; j++) {

	    if (j%numThreads == myID) {

	      for (int k=0; k<blockSize; k++) {

		if (currentCol < size) { // this is for case when size%blockSize != 0

		  int rowkTop = RowTop[currentCol];
		  double *alkPtr = topRowPtr[currentCol];
		  // double *alkPtr = &A[iDiagLoc[currentCol-1]];
		  int maxRowikTop;
		  if (rowkTop < startRow) {
		    alkPtr += startRow-rowkTop; // pointer to start of block row
		    maxRowikTop = startRow;
		  } else
		    maxRowikTop = rowkTop;

		  for (int l=maxRowikTop; l<=lastRow; l++) {
		    double tmp = *alkPtr;
		    int rowlTop = RowTop[l];
		    int maxRowklTop;
		    double *amlPtr, *amkPtr;
		    if (rowlTop < rowkTop) {
		      amlPtr = topRowPtr[l] + (rowkTop - rowlTop);
		      amkPtr = topRowPtr[currentCol];
		      maxRowklTop = rowkTop;
		    } else {
		      amlPtr = topRowPtr[l];
		      amkPtr = topRowPtr[currentCol] + (rowlTop - rowkTop);
		      maxRowklTop = rowlTop;
		    } 

		    for (int m = maxRowklTop; m<l; m++) 
		      tmp -= *amkPtr++ * *amlPtr++;
		  
		    *alkPtr++ = tmp;
		  }
		  currentCol++;
		  if (currentCol > lastColEffected) {
		    k = blockSize;
		    j = nBlck;
		  }
		  
		} else
		  k = blockSize;
	      }
	    } else
	      currentCol += blockSize;  
	  }
	}

	// update the data for the next block
	startRow += blockSize;
	lastRow = startRow + blockSize -1;
	lastColEffected = lastRow + maxColHeight -1;
      }

      // signal the main thread that this thread is done with the work
      mutex_lock(&tcb->end_mutex);
      tcb->threadsDone++;

      cond_signal(&tcb->end_cond);
      mutex_unlock(&tcb->end_mutex);
    }
}
