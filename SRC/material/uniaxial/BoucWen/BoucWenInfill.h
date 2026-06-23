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
#pragma once

// Description: This file contains the class definition for 
// BoucWenInfill material. BoucWenInfill material provides 
// a hysteretic uniaxial material described by the Bouc-Wen law 
// with stiffness, strength degradation and pinching. 
// Particularly suitable to simulate hysteresis of infill panels.
//
// Reference: Sirotti, S., Pelliciari, M., Di Trapani, F.,
// Briseghella, B., Carlo Marano, G., Nuti, C., & Tarantino, A. M. (2021).
// Development and validation of new bouc�wen data-driven hysteresis model 
// for masonry infilled rc frames. 
// Journal of Engineering Mechanics, 147(11), 04021092.
//
// Written by Stefano Sirotti (stefano.sirotti@unimore.it)
// Created on January 2022
//

#include <UniaxialMaterial.h>

class BoucWenInfill : public UniaxialMaterial
{
  public:
    BoucWenInfill(int tag,
     double mass,
	 double alpha,
	 double beta0,
	 double eta0,
	 double n,
	 double k,
	 double xy,
	 double deltak,
	 double deltaf,
	 double psi,
	 double Zs,
	 double As,
	 double epsp,
	 double tolerance,
	 int maxNumIter,
     double density);
	
    BoucWenInfill();
    ~BoucWenInfill();  

    const char *getClassType() const {return "BoucWenInfill";}

    int setTrialStrain(double strain, double strainRate = 0.0); 
    double getStrain();          
    double getStress();
    double getTangent();
    double getInitialTangent();

    int commitState();
    int revertToLastCommit();    
    int revertToStart();
    
    UniaxialMaterial *getCopy();
    int sendSelf(int commitTag, Channel &);  
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);    
    void Print(OPS_Stream &s, int flag);

  private:

    double signum(double);

    // Material parameters
    double mass;
	double alpha;
	double beta0;
	double eta0;
	double n;
	double k;
	double xy;
	double deltak;
	double deltaf;
	double psi;
	double Zs;
	double As;
	double epsp;
    
    double tolerance;
    int maxNumIter;
    double density;
    // History variables (trial and committed)
	double xmaxp;
    double xmax;
	double Tstrain, Cstrain;
    double Tz, Cz;
    double Te, Ce;
    
    // Other variables
    double Tstress, Ttangent;
};
