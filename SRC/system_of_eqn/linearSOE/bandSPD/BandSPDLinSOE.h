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
// Created: February 1997
//
// Description: This file contains the class definition for BandSPDLinSOE
// BandSPDLinSOE is a subclass of LinearSOE. It uses the LAPACK Upper storage
// scheme to store the components of the A matrix.
//
//
#pragma once
#include <LinearSOE.h>
#include <Vector.h>

class BandSPDLinSolver;

class BandSPDLinSOE : public LinearSOE
{
  public:
    BandSPDLinSOE(BandSPDLinSolver &theSolver);    
    BandSPDLinSOE(int classTag);
    BandSPDLinSOE(BandSPDLinSolver &theSolver, int classTag);

    virtual ~BandSPDLinSOE();

    virtual int getNumEqn() const;
    virtual int setSize(Graph &) override;

    virtual int addA(const Matrix &, const ID &, double fact = 1.0) override;

    virtual int addB(const Vector &, const ID &, double fact = 1.0);    
    virtual int setB(const Vector &, double fact = 1.0);        
    
    virtual void zeroA();
    virtual void zeroB();

    const Vector &getX() override;
    virtual const Vector &getB() override;

    virtual void setX(int loc, double value);    
    virtual void setX(const Vector &x);
    
    friend class BandSPDLinSolver;
    friend class BandSPDLinLapackSolver;    
    friend class BandSPDLinThreadSolver;        
    
  protected:
    int size, half_band;    
    double *A;
    Vector B, X; 
    int Asize;
    int aFactored;
    bool factored;
    
  private:
};
