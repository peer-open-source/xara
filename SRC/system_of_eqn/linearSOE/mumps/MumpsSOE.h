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
                                                                        
// $Revision: 1.7 $
// $Date: 2009-05-11 20:56:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/mumps/MumpsSOE.h,v $
                                                                        
#ifndef MumpsSOE_h
#define MumpsSOE_h

// Written: fmk 
// Created: 02/06
//
// Description: This file contains the class definition for MumpsSOE
// MumpsSOE is a subclass of LinearSOE. It stores the matrix equation
// Ax=b using the sparse column-compacted storage scheme for storing the 
// matrix A. 
//
// matrix types (matType): 0 Unsymmetrc
//                         1 Symmetrix positive definite
//                         2 General Symmetric
//
// What: "@(#) MumpsSOE.h, revA"

#include <LinearSOE.h>
#include <Vector.h>

class MumpsSolver;
class MumpsParallelSolver;
class LinearSOESolver;

enum class MumpsMatrixType: int {
  General = 0, 
  SPD = 1, 
  SID = 2
};

class MumpsSOE : public LinearSOE
{
  public:
    MumpsSOE(MumpsSolver &, int matType=0);        
    MumpsSOE();       
    MumpsSOE(int classTag);        
    MumpsSOE(LinearSOESolver &, int classTag, int matType = 0);        

    virtual ~MumpsSOE();

    int getNumEqn() const override;
    int setSize(Graph &theGraph) override;
    int addA(const Matrix &, const ID &, double fact = 1.0) override;
    int addB(const Vector &, const ID &, double fact = 1.0) override;    
    int setB(const Vector &, double fact = 1.0) override; 

    void zeroA() override;
    void zeroB() override;
    
    const Vector &getX() override;
    const Vector &getB() override;    
    double normRHS() override;

    void setX(int loc, double value) override;        
    void setX(const Vector &x) override;        
    virtual int setMumpsSolver(MumpsSolver &);    

    int sendSelf(int commitTag, Channel &) override;
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &) override;    

    friend class MumpsSolver;    
    friend class MumpsParallelSolver;    

  protected:
    int size;            // order of A
    int nnz;             // number of non-zeros in A
    double *A, *B, *X;   // 1d arrays containing coefficients of A, B and X
    int *colA, *rowA, *rowB, *colStartA; // int arrays containing info about coeficientss in A
    Vector *vectX;
    Vector *vectB;    
    int Asize, Bsize;    // size of the 1d array holding A
    bool factored;
    int matType;

  private:
};


#endif

