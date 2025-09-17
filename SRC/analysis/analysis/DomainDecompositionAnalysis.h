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
// Created: Tue Sept 17 16:34:47: 1996
// Revision: A
//
// Description: This file contains the class definition for 
// DomainDecompositionAnalysis. DomainDecompositionAnalysis is a subclass 
// of AnalysisAnalysis, it is used when performing a domain decomposition
// analysis. It provides methods which can be invoked by a subdomain to
// perform the numerical computations required.
//
// What: "@(#) DomainDecompositionAnalysis.h, revA"

#ifndef DomainDecompositionAnalysis_h
#define DomainDecompositionAnalysis_h

#include <Analysis.h>
#include <Matrix.h>
#include <Vector.h>
#include <MovableObject.h>

class ConstraintHandler;
class DOF_Numberer;
class AnalysisModel;
class IncrementalIntegrator;
class LinearSOE;
class EigenSOE;
class DomainSolver;
class DomainDecompAlgo;
class Subdomain;
class Vector;
class EquiSolnAlgo;
class ConvergenceTest;

class DomainDecompositionAnalysis: public Analysis, public MovableObject
{
  public:
    DomainDecompositionAnalysis(Subdomain &theDomain);

    DomainDecompositionAnalysis(int classTag, 
                                Subdomain &theDomain);    
    
    DomainDecompositionAnalysis(Subdomain &theDomain,
                                ConstraintHandler &theHandler,
                                DOF_Numberer &theNumberer,
                                AnalysisModel &theModel,
                                DomainDecompAlgo &theSolnAlgo,                   
                                IncrementalIntegrator &theIntegrator,        
                                LinearSOE &theSOE,
                                DomainSolver &theSolver,
                                ConvergenceTest *theTest);



    virtual ~DomainDecompositionAnalysis();
    virtual void clearAll();            
    virtual int initialize();
    virtual int domainChanged() =0;

    // methods for non standard domain deomposition analysis
    virtual bool doesIndependentAnalysis();    
    virtual int analyze(double dT);

    // methods for standard domain deomposition analysis
    // that do some form of condensation to the tangent
            int  getNumExternalEqn();
            int  getNumInternalEqn();

//  virtual int  newStep(double dT);
    virtual int  analysisStep(double dT) = 0;
    virtual int  eigenAnalysis(int numMode, bool generalized, bool findSmallest);
    virtual int  computeInternalResponse();
    virtual int  formTangent() final;
    virtual int  formResidual() final;
    virtual int  formTangVectProduct(Vector &force) final;
    const Matrix& getTangent();
    const Vector& getResidual();
    const Vector& getTangVectProduct();

    virtual int sendSelf(int commitTag, Channel &);
    virtual int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);

    // methods to change the analysis aggregates
    virtual int setAlgorithm(EquiSolnAlgo &theAlgorithm);
    virtual int setIntegrator(IncrementalIntegrator &theIntegrator);
    virtual int setLinearSOE(LinearSOE &theSOE);
    virtual int setEigenSOE(EigenSOE &theSOE);
    virtual int setConvergenceTest(ConvergenceTest &theTest);
    
  protected: 
    Subdomain                *getSubdomainPtr() const;
    ConstraintHandler         *getConstraintHandlerPtr() const;
    DOF_Numberer         *getDOF_NumbererPtr() const;
    AnalysisModel         *getAnalysisModelPtr() const;
    DomainDecompAlgo    *getDomainDecompAlgoPtr() const;        
    IncrementalIntegrator  *getIncrementalIntegratorPtr() const;    
    LinearSOE                 *getLinSOEPtr() const;
    DomainSolver        *getDomainSolverPtr() const;

    Channel *myChannel;
    int checkAllResult(int);
    
  private:
    Subdomain                      *theSubdomain;
    ConstraintHandler              *theHandler;    
    DOF_Numberer              *theNumberer;
    AnalysisModel              *theModel;
    DomainDecompAlgo              *theAlgorithm;
    IncrementalIntegrator    *theIntegrator;        

    LinearSOE                      *theSOE;
    DomainSolver             *theSolver;
    ConvergenceTest          *theTest;

    Vector                      *theResidual;
    int numEqn;
    int numExtEqn;

    // the following 2 variables are used to allow formResidual()
    // and formTangVectProduct() to be called before formTangent()
    // - this must be allowed as typical elements will not have to fromTangent
    // before being asked to form Residual(). 
    bool tangFormed;
    int tangFormedCount; // saves the expense of computing formTangent() 
                               // for same state of Subdomain.
    int domainStamp;                           
};

#endif

