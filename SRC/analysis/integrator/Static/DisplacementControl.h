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
// Description: This file contains the class definition for DisplacementControl.
// DisplacementControl is an algorithmic class for performing a static analysis
// using the arc length scheme, that is within a load step the following
// constraint is enforced: 
//  i=1        delta U^T delta U + alpha^2 delta lambda^2 = delta s^2
//  i>1        dU^T delta U + alpha^2 dLambda delta lambda = 0
// where dU is change in nodal displacements for step, dLambda is
// change in applied load and DisplacementControl is a control parameter.
//
// Written: fmk 
// Created: 07/98
// Revision: A
//
#ifndef DisplacementControl_h
#define DisplacementControl_h

#include <StaticIntegrator.h>
#include <Vector.h>
class LinearSOE;
class AnalysisModel;
class FE_Element;
class Vector;
class Domain;

class DisplacementControl : public StaticIntegrator
{
   public:
      DisplacementControl(int node, int dof, double increment, Domain *theDomain,
			  int numIncrStep, double minIncrement, double maxIncrement,
			  int tangFlag = 0);

      ~DisplacementControl();

      int newStep() final;    
      int update(const Vector &deltaU) final;
      int domainChanged() final;
      int revertToLastStep() final {return 0;}

      int sendSelf(int commitTag, Channel &) final;
      int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &) final;

      void Print(OPS_Stream &s, int flag) final;   

      // Sensitivity
//    int formEleResidual(FE_Element *theEle);
      int saveSensitivity(const Vector &v, int gradNum, int numGrads);
      int saveLambdaSensitivity(double dlambdadh, int gradNum, int numGrads);
      int commitSensitivity(int gradNum, int numGrads);
      int computeSensitivities();

      /////////////////////////////// Abbas //////////////////////////////
      Vector *formTangDispSensitivity(Vector *dUhatdh,int gradNumber); //Obtain *dUhatdh 
     // Vector *formResidualDispSensitivity(Vector *dUIJdh,int gradNumber);// Obtain dKdh*deltaUbar
      double formdLambdaDh(int gradNumber);//calculate dLambdadh for J=1
      double getLambdaSensitivity(int gradNumber);// update the dLambdadh for J>1
      bool computeSensitivityAtEachIteration();// A key that return 1 for loadControl and 2 for DisplacementControl


   protected:

      double dlambdadh; // deltaLambdaI1 for the first iteration J=1
      double Dlambdadh;// deltaLambdaIJ: for J>1
      double dLambda;
      double calldLambda1dh;//Abbas
      int CallParam;

   private:
      int formSensitivityRHS(int gradNum);
      int formIndependentSensitivityRHS();
   
      int theNode;          // the node that is being followed
      int theDof;           // the dof at the node being followed
      double theIncrement;  // deltaU at step (i)
      Domain *theDomain;    // the domain containing the node being followed
      int theDofID;         // the system level id of the dof being followed

      Vector *deltaUhat, *deltaUbar, *deltaU,*phat,*deltaUstep,*dphatdh;
      //  Vector *deltaUhat_newStep ;
      Vector *dLAMBDAdh; 


      // Pointers used for sensitivity analysis
      Vector
         *dUhatdh     , // The derivative of the tangent displacement w/r to parameter h
         *sensU       , // Displacement sensitivity using displacement control scheme
         *d_deltaU_dh , // The derivative of the residual displacement
         *dUIJdh      , // The sensitivity of the residual displacement
         *Residual    , // the residual forces that are required to obtain the dLambdadh
         *Residual2   ; // the residual forces required to obtain dUdh


      double deltaLambdaStep, currentLambda;  // dLambda(i) & current value of lambda  
      double dLambdaStepDh ;//Abbas
      double specNumIncrStep, numIncrLastStep; // Jd & J(i-1) 
      double minIncrement, maxIncrement; // min/max values of deltaU at (i)

      int tangFlag;

      // adding sensitivity
//    int gradNumber;
//    int sensitivityFlag;
//    FE_Element *theEle;
};

#endif

