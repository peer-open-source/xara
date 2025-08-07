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
                                                                        
// $Revision: 1.7 $
// $Date: 2008-08-26 16:30:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ElasticMaterial.h,v $
                                                                        
                                                                        
#ifndef ElasticMaterial_h
#define ElasticMaterial_h

// Written: fmk 
// Created: 07/98
// Revision: A
//
// Description: This file contains the class definition for 
// ElasticMaterial. ElasticMaterial provides the abstraction
// of an viscoelastic uniaxial material,
// i.e. stress = E*strain + eta*strainrate


#include <UniaxialMaterial.h>

class ElasticMaterial : public UniaxialMaterial
{
  public:
    ElasticMaterial(int tag, double E);
    ElasticMaterial(int tag, double Epos, double eta, double Eneg, double density = 0.0);
    ElasticMaterial();
    ~ElasticMaterial();

    const char *getClassType() const {return "ElasticMaterial";};

    int setTrialStrain(double strain, double strainRate = 0.0); 
    int setTrial(double strain, double &stress, double &tangent, double strainRate = 0.0); 
    double getStrain() {return trialStrain;};
    double getStrainRate() {return trialStrainRate;};
    double getStress();
    double getTangent();
    double getDampTangent() {return eta;};
    double getInitialTangent();
    double getRho() final;

    int commitState();
    int revertToLastCommit();    
    int revertToStart();        

    UniaxialMaterial *getCopy();
    
    int sendSelf(int commitTag, Channel &) final;  
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &) final;
    
    void Print(OPS_Stream &s, int flag) final;
    
    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int activateParameter(int parameterID);
    double getStressSensitivity(int gradIndex, bool conditional);
    double getTangentSensitivity(int gradIndex);
    double getInitialTangentSensitivity(int gradIndex);
    int commitSensitivity(double strainGradient, int gradIndex, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////

  protected:
    
  private:
    double trialStrain;
    double trialStrainRate;
    double committedStrain;
    double committedStrainRate;
    double Epos;
    double Eneg;
    double eta;
    double density;
    int parameterID;
};


#endif

