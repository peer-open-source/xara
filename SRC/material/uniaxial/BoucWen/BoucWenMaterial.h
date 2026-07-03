/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.5 $
// $Date: 2008-08-26 16:22:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/BoucWenMaterial.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu) 
//

#ifndef BoucWenMaterial_h
#define BoucWenMaterial_h

#include <UniaxialMaterial.h>
#include <Matrix.h>

class BoucWenMaterial : public UniaxialMaterial
{
  public:
    BoucWenMaterial(int tag, 
		    double alpha,
		    double ko,
		    double n,
		    double gamma,
		    double beta,
		    double Ao,
		    double deltaA,
		    double deltaNu,
		    double deltaEta,
		    double tolerance,
		    int maxNumIter);
    BoucWenMaterial();
    ~BoucWenMaterial();

    const char *getClassType() const {return "BoucWenMaterial";}

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain();          
    double getStress();
    double getTangent();

    int commitState();
    int revertToLastCommit();    
    int revertToStart();
        
    UniaxialMaterial *getCopy();
    int sendSelf(int commitTag, Channel &);  
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);    
    void Print(OPS_Stream &s, int flag);
    
	// Reliability and sensitivity stuff
    double getInitialTangent();
    int setParameter (const char **argv, int argc, Parameter &param);
    int    updateParameter          (int parameterID, Information &info);
	int    activateParameter        (int parameterID);
	double getStressSensitivity     (int gradIndex, bool conditional);
	double getStrainSensitivity     (int gradIndex);
	double getTangentSensitivity    (int gradIndex);
	double getDampTangentSensitivity(int gradIndex);
	double getRhoSensitivity        (int gradIndex);
	int    commitSensitivity        (double strainGradient, int gradIndex, int numGrads);
    double  getInitialTangentSensitivity(int gradIndex);

  private:
    // Material parameters
    double alpha;
    double ko;
	double n;
    double gamma;
    double beta;
    double Ao;
    double deltaA;
    double deltaNu;
    double deltaEta;

    // History variables (trial and committed)
    double Tstrain, Cstrain;
	double Tz, Cz;
	double Te, Ce;

	// Ohter variables
	double Tstress, Ttangent;

	double tolerance;
	int maxNumIter;

	// Sensitivit stuff
    int parameterID;
	Matrix *SHVs;
};


#endif

