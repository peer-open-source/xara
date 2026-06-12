//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the class definition for SubdomainFE.
//
// Written: fmk 
// Created: 11/96
// Revision: A
//
#ifndef SubdomainFE_h
#define SubdomainFE_h

#include <ID.h>
#include <FE_Element.h>

class Vector;
class Matrix;
class Element;
class Integrator;
class StaticIntegrator;
class AnalysisModel;


class SubdomainFE: public FE_Element
{
  public:
    SubdomainFE(int tag, Element *theElement);
    ~SubdomainFE() override;

    // public methods for setting/obtaining mapping information
    virtual int setID(AnalysisModel &);
    virtual const ID &getID() const;

    // methods to form and obtain the tangent and residual
    virtual const Matrix &getTangent(Integrator *);
    virtual const Vector &getResidual(Integrator *);

    // methods called by integrator to build tangent
    virtual void  zeroTangent()                  ;
    virtual void  addKtToTang(double fact = 1.0) ;
    virtual void  addKiToTang(double fact = 1.0) ;
    virtual void  addCtoTang (double fact = 1.0) ;
    virtual void  addMtoTang (double fact = 1.0) ;
    virtual void  addKpToTang(double fact = 1.0, int numP = 0);
    virtual int   storePreviousK(int numP);

    // methods used by integrator to build residual    
    virtual void  zeroResidual();    
    virtual void  addRtoResidual(double fact = 1.0);
    virtual void  addRIncInertiaToResidual(double fact);    

    // methods for ele-by-ele strategies
    virtual const Vector &getTangForce(const Vector &x, double fact = 1.0);
    virtual const Vector &getK_Force(const Vector &x, double fact = 1.0);
    virtual const Vector &getKi_Force(const Vector &x, double fact = 1.0);
    virtual const Vector &getC_Force(const Vector &x, double fact = 1.0);
    virtual const Vector &getM_Force(const Vector &x, double fact = 1.0);
    virtual void  addM_Force(const Vector &accel, double fact = 1.0);    
    virtual void  addD_Force(const Vector &vel, double fact = 1.0);    
    virtual void  addK_Force(const Vector &disp, double fact = 1.0);

    virtual int updateElement();

    virtual Integrator   *getLastIntegrator();
    virtual const Vector &getLastResponse();
    virtual Element *getElement();

    virtual void  Print(OPS_Stream&, int flag) {return;}

    virtual void addResistingForceSensitivity(int gradNumber, double fact = 1.0);
    virtual void addM_ForceSensitivity       (int gradNumber, const Vector &vect, double fact = 1.0);
    virtual void addD_ForceSensitivity       (int gradNumber, const Vector &vect, double fact = 1.0);
    virtual int  commitSensitivity           (int gradNum, int numGrads);
   
  protected:
    // protected variables
    ID myID;

  private:
    // private variables    
    int numDOF;
    Element       *myEle;
    Vector        *theResidual;
    Matrix        *theTangent;
    Integrator    *theIntegrator; // need for Subdomain

    //
    // static variables
    //
    static Matrix **theMatrices; // array of pointers to class wide matrices
    static Vector **theVectors;  // array of pointers to class widde vectors
    static int numFEs;           // number of objects

};

#endif

