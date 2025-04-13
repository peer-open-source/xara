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
// Written: MHS
// Created: May 2000
//
// Description: This file contains the class definition for 
// HardeningMaterial2.  HardeningMaterial2 provides the abstraction
// for a one-dimensional rate-independent plasticity model
// with combined isotropic and kinematic hardening.
//
#ifndef HardeningMaterial2_h
#define HardeningMaterial2_h

#include <UniaxialMaterial.h>
class Matrix;

class HardeningMaterial2 : public UniaxialMaterial
{
 public:
  HardeningMaterial2(int tag, double E, double sigmaY,
                     double K, double H, double eta = 0.0);
  HardeningMaterial2();
  ~HardeningMaterial2();
  
  int setTrialStrain(double strain, double strainRate = 0.0); 
  double getStrain();          
  double getStress();
  double getTangent();
  double getInitialTangent() {return E;};
  
  int commitState();
  int revertToLastCommit();    
  int revertToStart();        
  
  UniaxialMaterial *getCopy();
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, 
               FEM_ObjectBroker &theBroker);    
  
  void Print(OPS_Stream &s, int flag =0);
  
  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);
  
  // AddingSensitivity:BEGIN //////////////////////////////////////////
  int activateParameter(int parameterID);
  double getStressSensitivity(int gradNumber, bool conditional);
  double getInitialTangentSensitivity(int gradNumber);
  int commitSensitivity(double strainGradient, int gradNumber, int numGrads);
  // AddingSensitivity:END ///////////////////////////////////////////
  
 protected:
  
 private:
  // Material parameters
  double E;           // Elastic modulus
  double sigmaY;      // Yield stress
  double Hiso;        // Isotropic hardening parameter
  double Hkin;        // Kinematic hardening parameter
  double eta;
  
  // Committed history variables
  double CplasticStrain;        // Committed plastic strain
  double CbackStress;           // Committed back stress;
  double Chardening;            // Committed internal hardening variable
  
  // Trial history variables
  double TplasticStrain;        // Trial plastic strain
  double TbackStress;           // Trial back stress
  double Thardening;            // Trial internal hardening variable
  
  // Trial state variables
  double Tstrain;               // Trial strain
  double Tstress;               // Trial stress
  double Ttangent;              // Trial tangent
  
  // AddingSensitivity:BEGIN //////////////////////////////////////////
  int parameterID;
  Matrix *SHVs;
  // AddingSensitivity:END ///////////////////////////////////////////
  
  double getStressGradient(int gradNumber);
  int setStrainGradient(int gradNumber, double depsilondh);
};

#endif
