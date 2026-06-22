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
                                                                        
// $Revision: 1.6 $
// $Date: 2006-08-04 18:18:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticIsotropicPlaneStress2D.h,v $

#ifndef ElasticIsotropic3DThermal_h
#define ElasticIsotropic3DThermal_h

// Written: fmk
// Created: 10/11
//
// Description: 
//
// What: "@(#) ElasticIsotropicThreeDimesnional.h, revA"

//Modified for SIF modelling by Liming Jiang [http://openseesforfire.github.io]

#include <ElasticIsotropicMaterialThermal.h>
#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class ElasticIsotropic3DThermal : public ElasticIsotropicMaterialThermal
{
  public:
    ElasticIsotropic3DThermal(int tag, double e, double nu, double rho=0,double alpha=0.0, int softindex = 0);
    ElasticIsotropic3DThermal();
    ~ElasticIsotropic3DThermal();

    const char *getClassType(void) const {return "ElasticIsotropic3DThermal";}

    int setTrialStrain(const Vector &v);
    int setTrialStrain(const Vector &v, const Vector &r);
    int setTrialStrainIncr(const Vector &v);
    int setTrialStrainIncr(const Vector &v, const Vector &r);
    const Matrix &getTangent();
    const Matrix &getInitialTangent();

    double setThermalTangentAndElongation(double &TempT, double &, double &);//J.Jiang add
    const Vector& getTempAndElong();

    const Vector &getStress();
    const Vector &getStrain();
    
    int commitState();
    int revertToLastCommit();
    int revertToStart();
    
    NDMaterial *getCopy();
    const char *getType() const;
    int getOrder() const;

    int sendSelf(int commitTag, Channel &);  
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);    
    
 protected:
	
  private:
    static Vector sigma;	// Stress vector ... class-wide for returns
    static Matrix D;		// Elastic constants
    Vector epsilon;	        // Trial strains
    Vector Cepsilon;	        // Committed strain
	
    int softIndex;
	double Temp;  //Temperature
	double ThermalElong;  // eps(theata) = alpha *temperature
	double E0T;//Elasticity modulus at temperature T
	double E;
	double Alpha; //Coefficient of thermal exmapnsion
	double* redfactors;
};

#endif
