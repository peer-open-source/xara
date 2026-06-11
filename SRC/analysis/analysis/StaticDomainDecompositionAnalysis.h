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
// Revision: A
//
// Description: This file contains the class definition for 
// StaticDomainDecompositionAnalysis.StaticDomainDecompositionAnalysis is a subclass 
// of DomainDecompositionAnalysis used to perform a static analysis step on a subdomain.

//
// What: "@(#) StaticDomainDecompositionAnalysis.h, revA"

#ifndef StaticDomainDecompositionAnalysis_h
#define StaticDomainDecompositionAnalysis_h

#include <DomainDecompositionAnalysis.h>
#include <Matrix.h>
#include <Vector.h>
#include <MovableObject.h>

class ConstraintHandler;
class DOF_Numberer;
class AnalysisModel;
class StaticIntegrator;
class ConvergenceTest;
class LinearSOE;
class EigenSOE;
class LinearSolver;
class EquiSolnAlgo;
class Subdomain;
class ConvergenceTest;

class StaticDomainDecompositionAnalysis: public DomainDecompositionAnalysis
{
  public:
    StaticDomainDecompositionAnalysis(Subdomain &theDomain);
    StaticDomainDecompositionAnalysis(Subdomain &theDomain,
				      ConstraintHandler &theHandler,
				      DOF_Numberer &theNumberer,
				      AnalysisModel &theModel,
				      EquiSolnAlgo &theSolnAlgo,		   
				      LinearSOE &theSOE,
				      StaticIntegrator &theIntegrator,
				      ConvergenceTest *theTest,
				      bool setLinks = true);

    ~StaticDomainDecompositionAnalysis();
    void clearAll();	    
    int initialize();
    int domainChanged();

    // methods for non standard domain deomposition analysis
    int analyze(double dT);
    int eigen(int numMode, bool generalized, bool findSmallest);
    bool doesIndependentAnalysis();    

    // methods for standard domain deomposition analysis
    // that do some form of condensation to the tangent
    virtual int  analysisStep(double dT) final;
    int  eigenAnalysis(int numMode, bool generalized, bool findSmallest);
    int  computeInternalResponse();

    // methods to change the analysis aggregates
    int setAlgorithm(EquiSolnAlgo &theAlgorithm);
    int setIntegrator(IncrementalIntegrator &theIntegrator);
    int setLinearSOE(LinearSOE &theSOE);
    int setEigenSOE(EigenSOE &theSOE);
    int setConvergenceTest(ConvergenceTest &theTest);

    // methods to send/receive
    int sendSelf(int commitTag, Channel &);
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);

  private:
    ConstraintHandler 	*theConstraintHandler;    
    DOF_Numberer 	*theDOF_Numberer;
    AnalysisModel 	*theAnalysisModel;
    EquiSolnAlgo 	*theAlgorithm;
    LinearSOE 		*theSOE;
    EigenSOE 		*theEigenSOE;
    StaticIntegrator    *theIntegrator;
    ConvergenceTest    *theTest;
    int domainStamp;			   
};

#endif











