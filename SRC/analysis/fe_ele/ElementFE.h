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
// Description: This file contains the class definition for ElementFE.
//
// Written: cmp 
// Created: Jan 2026
// Revision: A
//
#ifndef ElementFE_h
#define ElementFE_h

#include <ID.h>
#include <FE_Element.h>

class Vector;
class Matrix;
class Element;
class Integrator;
class StaticIntegrator;
class AnalysisModel;


class ElementFE: public FE_Element
{
  public:
    ElementFE(int tag, Element *);
    ~ElementFE() override;
    const char* getClassName() const override {return "ElementFE";}

    // public methods for setting/obtaining mapping information
    int  setID(AnalysisModel &) override;
    const ID &getID() const override;

    // methods to form and obtain the tangent and residual
    const Matrix &getTangent(Integrator *) override;
    const Vector &getResidual(Integrator *) override;

    // methods called by integrator to build tangent
    void  zeroTangent() override;
    void  addKtToTang(double fact = 1.0) override;
    void  addKiToTang(double fact = 1.0) override;
    void  addCtoTang (double fact = 1.0) override;
    void  addMtoTang (double fact = 1.0) override;
    virtual void  addKpToTang(double fact = 1.0, int numP = 0);
    virtual int   storePreviousK(int numP);

    // methods used by integrator to build residual    
    virtual void  zeroResidual();    
    virtual void  addRtoResidual(double fact = 1.0);
    void  addRIncInertiaToResidual(double fact) override;    

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
    Element *getElement() override;

    virtual void  Print(OPS_Stream&, int flag) {return;}

    virtual void addResistingForceSensitivity(int gradNumber, double fact = 1.0);
    virtual void addM_ForceSensitivity       (int gradNumber, const Vector &vect, double fact = 1.0);
    virtual void addD_ForceSensitivity       (int gradNumber, const Vector &vect, double fact = 1.0);
    virtual int  commitSensitivity           (int gradNum, int numGrads);
   
  protected:
    void  addLocalM_Force(const Vector &accel, double fact = 1.0);
    void  addLocalD_Force(const Vector &vel, double fact = 1.0);
    void  addLocalM_ForceSensitivity(int gradNumber, const Vector &accel, double fact = 1.0);
    void  addLocalD_ForceSensitivity(int gradNumber, const Vector &vel, double fact = 1.0);

    //
    ID myID;

  private:
    // private variables
    int numDOF;
    Element       &myEle;
    Vector        *theResidual, *vecY;
    Matrix        *theTangent;
    Integrator    *theIntegrator; // need for Subdomain
    bool          own_workspace; // flag to indicate if this object owns the workspace or is using class wide objects

    //
    // static variables
    //
    static Matrix **theMatrices; // array of pointers to class wide matrices
    static Vector **VecsX;       // array of pointers to class widde vectors
    static Vector **VecsY;       // array of pointers to class widde vectors
    static int numFEs;           // number of objects
};

#endif

