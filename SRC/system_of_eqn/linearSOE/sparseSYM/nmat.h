/*
 * File:  nmat.h
 * =============
 * altered to improve data access
 *
 * Originally written by:  David R. Mackay
 *
 * Modified by:
 *  Jun Peng (junpeng@stanford.edu)
 *  Prof. Kincho H. Law
 *  Stanford University
 * --------------------
 */


#ifndef nmat_h
#define nmat_h
#include "FeStructs.h" // OFFDBLK

int  pfsfct(int neqns, double *diag, double **penv, int nblks, 
            int *xblk, OFFDBLK **begblk, OFFDBLK *first, int *rowblks);

void pfsslv(int neqns, double *diag, double **penv, int nblks, 
            int *xblk, double *rhs, OFFDBLK **begblk);

#endif
