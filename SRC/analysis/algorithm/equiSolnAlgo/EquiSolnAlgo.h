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
#ifndef EquiSolnAlgo_h
#define EquiSolnAlgo_h
//
// Description: This file contains the class definition for 
// EquiSolnAlgo. EquiSolnAlgo is an abstract base class, 
// i.e. no objects of it's type can be created.  Its subclasses deifine
// the sequence of operations to be performed in the analysis by static
// equilibrium of a finite element model.  
// 
// Written: fmk 
// Created: 11/96 
// Revision: A 
//
#include <SolutionAlgorithm.h>
#include <logging/Logging.h> // TODO
#include <IncrementalIntegrator.h>
class AnalysisModel;
class LinearSOE;
class ConvergenceTest;
class IncrementalIntegrator;

class EquiSolnAlgo: public SolutionAlgorithm
{
  public:
    EquiSolnAlgo(int classTag);
    virtual ~EquiSolnAlgo();

    // public functions defined for subclasses
    void 
    setLinks( IncrementalIntegrator &,
              LinearSOE &,
              ConvergenceTest *);

    virtual int solveCurrentStep() =0;

    virtual int setConvergenceTest(ConvergenceTest *theNewTest) final;    
    virtual ConvergenceTest *getConvergenceTest() final;     


    virtual void Print(OPS_Stream &, int flag) const =0;    

    virtual int    getNumFactorizations() const {return 0;}
    virtual int    getNumIterations()  const {return 0;}
    virtual double getTotalTimeCPU()   const {return 0.0;}
    virtual double getTotalTimeReal()  const {return 0.0;}
    virtual double getSolveTimeCPU()   const {return 0.0;}
    virtual double getSolveTimeReal()  const {return 0.0;}
    virtual double getAccelTimeCPU()   const {return 0.0;}
    virtual double getAccelTimeReal()  const {return 0.0;}
 
    // the following are not protected as convergence test
    // may need access to them
    // AnalysisModel           *getAnalysisModelPtr() const;
    IncrementalIntegrator   *getIncrementalIntegratorPtr() const;
    LinearSOE               *getLinearSOEptr() const;

  protected:
    ConvergenceTest *theTest;
    
  private:
    IncrementalIntegrator *theIntegrator;
    LinearSOE             *theSysOfEqn;
};

#endif


