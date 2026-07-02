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
                                                                        
#ifndef Steel01_h
#define Steel01_h

// Written: MHS 
// Created: 06/99
// Revision: A
//
// Description: This file contains the class definition for 
// Steel01.h
// 
//
#include <UniaxialMaterial.h>


class Steel01 : public UniaxialMaterial
{
  public:
    Steel01(int tag, double fy, double E0, double b,
       double a1, 
       double a2,
       double a3, 
       double a4,
       double density = 0.0
    );
    Steel01();
    ~Steel01();

    const char *getClassType() const {return "Steel01";}

    int setTrialStrain(double strain, double strainRate); 
    int setTrial (double strain, double &stress, double &tangent, double strainRate = 0.0);
    double getStrain();              
    double getStress();
    double getTangent();
    double getInitialTangent() {return E0;};

    int commitState();
    int revertToLastCommit();    
    int revertToStart();        

    UniaxialMaterial *getCopy();
    double getRho() override { return density; }
    
    int sendSelf(int commitTag, Channel &);  
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);    
    
    void Print(OPS_Stream &s, int flag);
    
    // Sensitivity
    int setParameter(const char **argv, int argc, Parameter &param);
    int    updateParameter          (int parameterID, Information &info);
    int    activateParameter        (int parameterID);
    double getStressSensitivity     (int gradIndex, bool conditional);
    double getInitialTangentSensitivity(int gradIndex);
    int    commitSensitivity        (double strainGradient, int gradIndex, int numGrads);

	// by SAJalali
	virtual double getEnergy() { return Energy; }

 private:
	// Material Properties
    double fy;  // Yield stress
    double E0;  // Initial stiffness
    double b;   // Hardening ratio (b = Esh/E0)
    double a1;
    double a2;
    double a3;
    double a4;      // a1 through a4 are coefficients for isotropic hardening
    double density; // mass per unit volume
    
	double Energy;	//by SAJalali

    /*** CONVERGED History Variables ***/
    double CminStrain;  // Minimum strain in compression
    double CmaxStrain;  // Maximum strain in tension
    double CshiftP;     // Shift in hysteresis loop for positive loading
    double CshiftN;     // Shift in hysteresis loop for negative loading
    int Cloading;       // Flag for loading/unloading
                        // 1 = loading (positive strain increment)
                        // -1 = unloading (negative strain increment)
                        // 0 initially

    /*** CONVERGED State Variables ***/    
    double Cstrain;
    double Cstress;
    double Ctangent;    

    /*** TRIAL History Variables ***/
    double TminStrain;
    double TmaxStrain;
    double TshiftP;
    double TshiftN;
    int Tloading;
    
    /*** TRIAL State Variables ***/
    double Tstrain;
    double Tstress;
    double Ttangent; // Not really a state variable, but declared here
                     // for convenience

    // Calculates the trial state variables based on the trial strain
    void determineTrialState (double dStrain);

    // Determines if a load reversal has occurred based on the trial strain
    void detectLoadReversal (double dStrain);

// AddingSensitivity:BEGIN //////////////////////////////////////////
    int parameterID;
	Matrix *SHVs;
// AddingSensitivity:END ///////////////////////////////////////////
};

#endif
