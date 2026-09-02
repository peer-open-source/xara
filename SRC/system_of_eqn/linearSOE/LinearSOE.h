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
// Created: 11/96
//
// Description: This file contains the class definition for LinearSOE.
// LinearSOE is an abstract base class and thus no objects of it's type
// can be instantiated. It has pure virtual functions which must be
// implemented in it's derived classes. LinearSystemOfEqn is an abstraction 
// of the linear system of equation given by : [A]{X} = {B} - {C},
// where A is a matrix and B,C and X are vectors. To solve the equation means
// given A, B and C to find the unknown X such that the equation is satisfied.
//
// What: "@(#) LinearSOE.h, revA"
#pragma once
#include <LinearAction.h>
#include <MovableObject.h>
#include <vector>

class LinearSOESolver;
class Graph;
class Matrix;
class Vector;
class ID;
class AnalysisModel;
class OPS_Stream;

class LinearSOE : public MovableObject
{
  public:
    LinearSOE(LinearSOESolver &, int classTag);    
    LinearSOE(int classTag);    
    virtual ~LinearSOE();

    virtual int solve();
    virtual int solve(const Vector& B, Vector& X) final;

    virtual int setInverseUpdate(LinearAction * update);
    virtual int setForwardUpdate(LinearAction * update);

    // pure virtual functions
    virtual int setSize(Graph &) =0;    
    virtual int getNumEqn() const =0;

    // virtual bool symmetric() const =0;
    
    virtual int addA(const Matrix &, const ID &, double fact = 1.0) =0;
    virtual int addB(const Vector &, const ID &, double fact = 1.0) =0;
    virtual int addX(const Vector &) { return -1; }
    virtual int setB(const Vector &, double fact = 1.0) =0;        

    virtual int addA(const Matrix &);

    virtual void zeroA() =0;
    virtual void zeroB() =0;

    virtual int formAp(const Vector &p, Vector &Ap) {return -1;};

    virtual const Vector &getX() = 0;
    virtual const Vector &getB() = 0;    
    virtual const Matrix *getA() {return nullptr;}
    virtual int   getA(int row, int col, double &value) const {return -2;}
  
    virtual double getDeterminant();
    // turn on determinant calculation in the solve; return false if not supported
    virtual bool   requireDeterminant();
  
    virtual double normRHS() {return this->getB().Norm();}

    virtual void setX(int loc, double value) =0;
    virtual void setX(const Vector &X) =0;

    virtual int saveSparseA(OPS_Stream& output, int baseIndex = 0); 
    virtual int getSparseA(ID& rowIndices, ID& colIndices, Vector& values, int baseIndex = 0);
    virtual int getSparseA(std::vector<int>& rowIndices, std::vector<int>& colIndices, std::vector<double>& values, int baseIndex = 0);
    
    
    virtual int setProcessID(int processTag) {return -1;}
    virtual int setChannels(int numChannels, Channel **theChannels) {return -1;}

  protected:
    LinearSOESolver *getSolver();
    int setSolver(LinearSOESolver &);

  private:
    LinearSOESolver *theSolver;
    LinearAction *m_fwd_update;
    LinearAction *m_inv_update;
};
