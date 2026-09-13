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
                                                                        
// $Revision: 1.1 $
// $Date: 2008/12/09 20:00:16 $
// $Source: /usr/local/cvs/OpenSees/PACKAGES/NewMaterial/cpp/ElasticPPcpp.h,v $
                                                                        
#ifndef CompressionDamage1d_h
#define CompressionDamage1d_h

#include <UniaxialMaterial.h>

class CompressionDamage1d : public UniaxialMaterial
{
  public:
    CompressionDamage1d(int _tag, double _E, double _fc, double _maxDuctility=-1.0, double r=1, bool triangular=false);    
    CompressionDamage1d();    

    ~CompressionDamage1d();

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain(void);          
    double getStress(void);
    double getTangent(void);

    double getInitialTangent(void);

    int commitState(void);
    int revertToLastCommit(void);    
    int revertToStart(void);    

    UniaxialMaterial *getCopy(void);
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
    
    void Print(OPS_Stream &s, int flag =0);
    
  protected:
    
  private:
    double fc;      // compressive strength
    double E;		// elastic modulus
	double maxDuctility; // max ductlity 
	double r;       // residual strength ratio
	bool triangular;// flag for triangular postpeak response 
	double ey;      // strain at compressive strength
    double D;		// damage variable
	double Dcommit; // committed damage variable

	double trialStrain;	     // trial strain
    double trialStress;      // current trial stress
    double trialTangent;     // current trial tangent

    double commitStrain;     // last commited strain
    double commitStress;     // last commited stress
    double commitTangent;    // last committed  tangent
};


#endif



