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
#pragma once
//
// Written: fmk
// Created: 10/11
//

#include <ElasticOrthotropicMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class ElasticOrthotropicThreeDimensional : public ElasticOrthotropicMaterial
{
  public:
    ElasticOrthotropicThreeDimensional(int tag, double Ex, double Ey, 
      double Ez, double vxy, double vyz, double vzx, double Gxy, double Gyz, 
      double Gzx, double rho = 0.0);
    ElasticOrthotropicThreeDimensional();
    ~ElasticOrthotropicThreeDimensional();

    const char *getClassType() const {return "ElasticOrthotropicThreeDimensional";}

    int setTrialStrain (const Vector &v);
    int setTrialStrain (const Vector &v, const Vector &r);
    int setTrialStrainIncr (const Vector &v);
    int setTrialStrainIncr (const Vector &v, const Vector &r);
    const Matrix &getTangent();
    const Matrix &getInitialTangent();
    
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

    const Vector& getStressSensitivity(int gradIndex, bool conditional);

  private:
    static Vector sigma;	// Stress vector ... class-wide for returns
    static Matrix D;		  // Elastic constants
    Vector epsilon;	      // Trial strains
    Vector Cepsilon;	    // Committed strain
};
