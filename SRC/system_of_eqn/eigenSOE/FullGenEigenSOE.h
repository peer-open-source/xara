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
#ifndef FullGenEigenSOE_h
#define FullGenEigenSOE_h
//
// Description: This file contains the class definition for
// FullGenEigenSOE, which stores full nonsymmetric matrices,
// A and M, for generalized eigenvalue computations.
//
// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 11/07
// Revision: A
//
#include <EigenSOE.h>
#include <Vector.h>

class AnalysisModel;
class FullGenEigenSolver;

class FullGenEigenSOE : public EigenSOE
{
public:
    FullGenEigenSOE(FullGenEigenSolver &, AnalysisModel &);

    virtual ~FullGenEigenSOE();

    virtual int getNumEqn() const;
    virtual int setSize(Graph &);

    int addA(const Matrix &, const ID &, double fact = 1.0) override;
    int addM(const Matrix &, const ID &, double fact = 1.0) override;    

    void zeroA() override;
    void zeroM() override;

    int sendSelf(int commitTag, Channel &) override;
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &) override;

    friend class FullGenEigenSolver;

private:
    int size;
    double *A;
    int Asize;
    double *M;
    int Msize;
    bool factored;
    AnalysisModel *theModel;
};

#endif
