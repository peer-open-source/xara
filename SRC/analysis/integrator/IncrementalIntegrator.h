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
#include <MovableObject.h>

class LinearSOE;
class AnalysisModel;
class ConvergenceTest;
class FE_Element;
class DOF_Group;
class Vector;
class OPS_Stream;

enum TangentFlag {
  CURRENT_TANGENT               =0,
  INITIAL_TANGENT               =1,
  CURRENT_SECANT                =2,
  INITIAL_THEN_CURRENT_TANGENT  =3,
  NO_TANGENT                    =4,
  SECOND_TANGENT                =5,
  HALL_TANGENT                  =6
};


class IncrementalIntegrator : protected Integrator, public MovableObject
{
  public:
    typedef int TangentFlagType; 

    IncrementalIntegrator(int classTag);
    virtual ~IncrementalIntegrator();

    void setLinks(AnalysisModel &,
                  LinearSOE &,
                  ConvergenceTest *);

    // methods to update the domain

    virtual int initialize() final;
//  virtual int solve(Algorithm&, ConvergenceTest&) =0;
    virtual int commit();
    virtual int revertToLastStep();
    virtual int revertToStart();
    virtual int domainChanged() {return 0;};
//  virtual int advance(double dT) =0;

    //
    // Invoked by the Algorithm
    //
    virtual int  update(const Vector &dx) =0;// correct;
    // Proxies for AnalysisModel
    virtual int  formUnbalance() = 0;

    virtual int  formTangent(int statusFlag = CURRENT_TANGENT);
    virtual int  formTangent(int statusFlag, 
                             double iFactor,
                             double cFactor);

    virtual double getCFactor();
    virtual const Vector &getVel();

    
    // Sensitivity interface
    // virtual int formResidualGradient(int grad);
    virtual int computeSensitivities(); // solveTangent
    virtual int formIndependentSensitivityLHS(int statusFlag = CURRENT_TANGENT);
    virtual  bool computeSensitivityAtEachIteration();
    int  setGradientType(int flag);
    bool shouldComputeAtEachStep(); // TODO(cmp) remove
    void activateSensitivityKey() {SensitivityKey=true;}
    bool activateSensitivity()    {return SensitivityKey;};


    virtual void Print(OPS_Stream &s, int flag =0) =0;
 
  protected:
    // These are almost final; they are invoked by 
    // [Transient,Static]Integrator::formUnbalance()
    virtual int  formNodalUnbalance() final;
    virtual int  formElementResidual() final;

    int addModalDampingForce(const Vector *modalDampingValues);
    int addModalDampingMatrix(const Vector *modalDampingValues);

  
    LinearSOE       *getLinearSOE() const;
    AnalysisModel   *getAnalysisModel() const;
    ConvergenceTest *getConvergenceTest() const;

    int statusFlag;
    double iFactor;
    double cFactor;


  private:
    // TODO(cmp): Move to TransientIntegrator
    // int setModalDampingFactors(const Vector &);
    int setupModal(const Vector *modalDampingValues);
    int doMv(const Vector &v, Vector &res);
    double   *eigenVectors;
    Vector   *eigenValues;
    Vector   *dampingForces;
    bool      isDiagonal;
    double   *diagMass;
    Vector   *mV;
    Vector   *tmpV1;
    Vector   *tmpV2;

    LinearSOE *theSOE;
    AnalysisModel *theAnalysisModel;
    ConvergenceTest *theTest;


    bool SensitivityKey;
    int analysisTypeTag;

    // method introduced for domain decomposition
    // This is private here because it should only be called by
    // classes using the `Integrator` interface (where it is public), 
    // not `IncrementalIntegrator`. The only reason it is not implemented
    // in Integrator is that it needs the LinearSOE
    // - cmp
    virtual int getLastResponse(Vector &result, const ID &id) final;
};

#endif

