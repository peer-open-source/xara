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

#ifndef SymmGeneralizedEigenSOE_h
#define SymmGeneralizedEigenSOE_h

#include <EigenSOE.h>

class Vector;
class AnalysisModel;
class SymmGeneralizedEigenSolver;


class SymmGeneralizedEigenSOE : public EigenSOE
{
public:
    SymmGeneralizedEigenSOE(SymmGeneralizedEigenSolver &, AnalysisModel &);

    virtual ~SymmGeneralizedEigenSOE();

    virtual int getNumEqn() const;
    virtual int setSize(Graph &);

    virtual int addA(const Matrix &, const ID &, double fact = 1.0);
    virtual int addM(const Matrix &, const ID &, double fact = 1.0);    

    virtual void zeroA();
    virtual void zeroM();

    int sendSelf(int commitTag, Channel &);
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);

    friend class SymmGeneralizedEigenSolver;

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
