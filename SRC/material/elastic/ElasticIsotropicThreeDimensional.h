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
#ifndef ElasticIsotropicThreeDimensional_h
#define ElasticIsotropicThreeDimensional_h

// Written: fmk
// Created: 10/11
//
// Description: 
//
// What: "@(#) ElasticIsotropicThreeDimesnional.h, revA"

#include <ElasticIsotropicMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

namespace OpenSees {

class ElasticIsotropicThreeDimensional : public ElasticIsotropicMaterial
{
  public:
    ElasticIsotropicThreeDimensional(int tag, double E, double nu, double rho);
    ElasticIsotropicThreeDimensional();
    ~ElasticIsotropicThreeDimensional();

    const char *getClassType() const {return "ElasticIsotropicThreeDimensional";}

    virtual int setTrialStrain (const Vector &v) final;
    virtual int setTrialStrain (const Vector &v, const Vector &r) final;
    virtual int setTrialStrainIncr (const Vector &v) final;
    virtual int setTrialStrainIncr (const Vector &v, const Vector &r) final;
    virtual const Matrix &getTangent() final;
    virtual const Matrix &getInitialTangent() final;
    
    const Vector &getStress() override;
    const Vector &getStrain() override;
    
    int commitState() override;
    int revertToLastCommit() override;
    int revertToStart() override;
    
    NDMaterial *getCopy() override;
    const char *getType() const override;
    int getOrder() const override;
    int sendSelf(int commitTag, Channel &);
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);
    
    const Vector& getStressSensitivity(int gradIndex, bool conditional) override;

 protected:

  private:
    static Vector sigma;	// Stress vector ... class-wide for returns
    static Matrix D;		// Elastic constants
    Vector epsilon;	        // Trial strains
    Vector Cepsilon;	        // Committed strain
};
} // namespace
#endif
