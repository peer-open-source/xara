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
// Description: This file contains the class definition for SolutionAlgorithm.
// SolutionAlgorithm is an abstract base class, i.e. no objects of it's
// type can be created. 
//
// Written: fmk 
// Created: 11/96
//
#ifndef SolutionAlgorithm_h
#define SolutionAlgorithm_h

#include <MovableObject.h>

class Channel;
class FEM_ObjectBroker;
class Recorder;

extern int SOLUTION_ALGORITHM_tangentFlag;

class SolutionAlgorithm: public MovableObject
{
  public:
    enum SolveError {
      BadAlgorithm    = -5, // Algo::solveCurrentStep          -5
      BadTestStart    = -5, // test->start()                   -5
      TestFailed      = -9, // result == CTest::Failure,       -3
      BadLinearSolve  = -3, // theSOE->solve(),                -3
      BadStepUpdate   = -4, // theIntegrator->update,          -4
      BadFormResidual = -2, // theIntegrator->formUnbalance(), -2
      BadFormTangent  = -1, // theIntegrator->formTangent(),   -1
    };

    SolutionAlgorithm(int classTag);
    virtual ~SolutionAlgorithm();

    // methods for monitoring the analysis during an algorithm
    virtual int  addRecorder(Recorder &theRecorder);    	
    virtual int  record(int track);    

    
  private:
    Recorder **theRecorders;
    int numRecorders;    
};

#endif


