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
                                                                       
// Created: Pedro Arduino, UW, 11.2011
//
// Description: This file contains the class definition for SAniSandMSPlaneStrain.

#ifndef SAniSandMSPS_h
#define SAniSandMSPS_h

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>

#include "SAniSandMS.h"

class SAniSandMSPlaneStrain : public SAniSandMS
{
  public : 

  //full constructor
  SAniSandMSPlaneStrain(int tag, double G0, double nu, double e_init, double Mc, double c, double lambda_c, double e0, double ksi,
        double P_atm, double m, double h0, double ch, double nb, double A0, double nd, double zeta, double mu0, 
          double beta,
          double mDen, 
        int integrationScheme = 2, int tangentType = 2, int JacoType = 1, double TolF = 1.0e-7, double TolR = 1.0e-7);
    
  //null constructor
  SAniSandMSPlaneStrain();

  //destructor
  ~SAniSandMSPlaneStrain();

  NDMaterial *getCopy(void);
  const char *getType(void) const;
  int getOrder(void) const;

  int setTrialStrain(const Vector &strain_from_element);

  // Unused trialStrain functions
  int setTrialStrain(const Vector &v, const Vector &r);
    
  const Vector& getStrain();
  const Vector& getEStrain();

  const Vector& getStress();
  const Vector& getStressToRecord();

  const Matrix& getTangent();
  const Matrix& getInitialTangent();
  
  //send back the state parameters
  const Vector  getState();

  private :

  // static vectors and matrices
  static Vector mSigma_M  ; // mSigma with continuum mechanic sign convention
  static Vector mEpsilon_M; // mEpsilon with continuum mechanic sign convention
  static Vector mEpsilonE_M; // mEpsilon with continuum mechanic sign convention
  static Vector rSigma;     // Stress for the recorders
  static Matrix mTangent;
  static Matrix mTangent_init;

};

#endif
