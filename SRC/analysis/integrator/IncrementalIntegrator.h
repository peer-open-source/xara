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
// Description: This file contains the interface for IncrementalIntegrator. 
// IncrementalIntegrator is an algorithmic class for setting up the finite 
// element equations in an incremental analysis and for updating the nodal
// response quantities based on the values in the soln vector.
//
// Written: fmk 
// Created: Tue Sept 17 15:54:47: 1996
// Revision: A
//
#ifndef IncrementalIntegrator_h
#define IncrementalIntegrator_h


#include <Integrator.h>

class LinearSOE;
// class EigenSOE;
class AnalysisModel;
class ConvergenceTest;
class FE_Element;
class DOF_Group;
class Vector;

enum TangentFlag {
 CURRENT_TANGENT               =0,
 INITIAL_TANGENT               =1,
 CURRENT_SECANT                =2,
 INITIAL_THEN_CURRENT_TANGENT  =3,
 NO_TANGENT                    =4,
 SECOND_TANGENT                =5,
 HALL_TANGENT                  =6
};

class IncrementalIntegrator : public Integrator
{
  public:
    typedef int TangentFlagType; 

    IncrementalIntegrator(int classTag);
    virtual ~IncrementalIntegrator();

    void setLinks(AnalysisModel &theModel,
                  LinearSOE &theSOE,
                  ConvergenceTest *theTest);

    // methods to set up the system of equations.
    // These are invoked by the Algorithm, and therefore must be 
    // in the IncrementalIntegrator interface.
    virtual int  update(const Vector &deltaU) =0;
    virtual int  formUnbalance() = 0;

    virtual int  formTangent(int statusFlag = CURRENT_TANGENT);
    virtual int  formTangent(int statusFlag, 
                             double iFactor,
                             double cFactor);    

    // methods to update the domain
//  virtual int newStep(double deltaT) =0;
    virtual int commit();
    virtual int revertToLastStep();
    virtual int initialize();

    virtual double getCFactor();

    virtual const Vector &getVel();
    int doMv(const Vector &v, Vector &res);


    // pure virtual methods to define the FE_ELe and DOF_Group contributions
//  virtual int formEleTangent(FE_Element *theEle)  =0;
//  virtual int formNodTangent(DOF_Group *theDof)   =0;    
//  virtual int formEleResidual(FE_Element *theEle) =0;
//  virtual int formNodUnbalance(DOF_Group *theDof) =0;    

// AddingSensitivity:BEGIN //////////////////////////////////
    virtual int revertToStart();
    virtual int formIndependentSensitivityLHS(int statusFlag = CURRENT_TANGENT);
// AddingSensitivity:END ////////////////////////////////////
 
  protected:
    // TODO(cmp): Move to TransientIntegrator
    // int setModalDampingFactors(const Vector &);
    int setupModal(const Vector *modalDampingValues);
    int addModalDampingForce(const Vector *modalDampingValues);
    int addModalDampingMatrix(const Vector *modalDampingValues);

    virtual int  formNodalUnbalance();
    virtual int  formElementResidual();

    LinearSOE       *getLinearSOE() const;
    AnalysisModel   *getAnalysisModel() const;
    ConvergenceTest *getConvergenceTest() const;

    int statusFlag;
    double iFactor;
    double cFactor;

    //    Vector *modalDampingValues;
    double   *eigenVectors;
    Vector   *eigenValues;
    Vector   *dampingForces;
    bool      isDiagonal;
    double   *diagMass;
    Vector   *mV;
    Vector   *tmpV1;
    Vector   *tmpV2;
    
  private:
    LinearSOE *theSOE;
    AnalysisModel *theAnalysisModel;
    ConvergenceTest *theTest;

    // method introduced for domain decomposition
    // This is private here because it should only be called by
    // classes using the `Integrator` interface (where it is public), 
    // not `IncrementalIntegrator`. The only reason it is not implemented
    // in Integrator is that it needs the LinearSOE
    // - cmp
    virtual int getLastResponse(Vector &result, const ID &id) final;

};

#endif

