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
// Description: This file contains the class definition for FullGenLinSOE
// FullGenLinSOE is a subclass of LinearSOE. It stores all the components
// of the linear system of equation in 1d arrays.
//
//
// Written: fmk 
// Created: 02/97
// Revision: A
//
#pragma once
#include <LinearSOE.h>
#include <Vector.h>

class FullGenLinSolver;

class FullGenLinSOE : public LinearSOE
{
  public:
    FullGenLinSOE(FullGenLinSolver &);        
    FullGenLinSOE(int N, FullGenLinSolver &);

    ~FullGenLinSOE();

    int getNumEqn() const;
    int setSize(Graph &theGraph);
    int addA(const Matrix &, const ID &, double fact = 1.0);
    int addB(const Vector &, const ID &, double fact = 1.0);    
    int setB(const Vector &, double fact = 1.0);
    
    void zeroA();
    void zeroB();
    
    int formAp(const Vector &p, Vector &Ap);

    const Vector &getX();
    const Vector &getB();    
    const Matrix *getA();

    void setX(int loc, double value);        
    void setX(const Vector &x);        

    friend class FullGenLinLapackSolver;

  private:
    int size;    
    double *A;
    Vector B, X;  
    Matrix *matA;
    bool factored;
};
