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
#ifndef FullGenEigenSolver_h
#define FullGenEigenSolver_h
//
// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 11/07
// Revision: A
//
// Description: This file contains the class definition for 
// FullGenEigenSolver. It computes the generalized eigenvalues
// and eigenvectors of a pair of real nonsymmetric matrices using
// the LAPACK subroutine DGGEV.

#include <EigenSolver.h>
#include <FullGenEigenSOE.h>

class FullGenEigenSolver : public EigenSolver
{
public:
    FullGenEigenSolver();    
    virtual ~FullGenEigenSolver();

    int solve(int numEigen, bool generalized, bool findSmallest = true) override;    
    int setSize() override;
    virtual int setEigenSOE(FullGenEigenSOE &);

    const Vector &getEigenvector(int mode) override;
    double getEigenvalue(int mode) override;

    int sendSelf(int commitTag, Channel &);
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);

private:
    void sort(int length, double *x, int *id);

    FullGenEigenSOE *theSOE;
    int numEigen;

    double *eigenvalue;
    double *eigenvector;
    int *sortingID;
    Vector *eigenV;
};

#endif
