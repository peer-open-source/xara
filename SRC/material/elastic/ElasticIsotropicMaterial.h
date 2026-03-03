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
                       
#ifndef ElasticIsotropicMaterial_h
#define ElasticIsotropicMaterial_h

// Written: MHS
// Created: Feb 2000
// Revision: A
//
// Description: This file contains the class definition for ElasticIsotropicMaterialModel.
// ElasticIsotropicMaterialModel is an abstract base class and thus no objects of it's type
// can be instantiated. It has pure virtual functions which must be
// implemented in it's derived classes. 
//
// What: "@(#) ElasticIsotropicMaterial.h, revA"

#include <NDMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

namespace OpenSees {
class ElasticIsotropicMaterial : public NDMaterial
{
  public:
    // Only called by subclasses to pass their tags to NDMaterialModel
    ElasticIsotropicMaterial (int tag, int classTag, double E, double nu, double rho = 0.0);

    // Called by clients
    ElasticIsotropicMaterial (int tag, double E, double nu, double rho = 0.0);

    // For parallel processing
    ElasticIsotropicMaterial (void);

    virtual ~ElasticIsotropicMaterial();

    virtual const char *getClassType() const {return "ElasticIsotropicMaterial";}

    virtual double getRho( ) ;

    virtual int setTrialStrain (const Vector &v);
    virtual int setTrialStrain (const Vector &v, const Vector &r);
    virtual int setTrialStrainIncr (const Vector &v);
    virtual int setTrialStrainIncr (const Vector &v, const Vector &r);
    virtual const Matrix &getTangent(void);
    virtual const Matrix &getInitialTangent(void);
    virtual const Vector &getStress();
    virtual const Vector &getStrain();

    virtual int commitState(void);
    virtual int revertToLastCommit (void);
    virtual int revertToStart();
    
    // Create a copy of material parameters AND state variables
    // Called by GenericSectionXD
    NDMaterial *getCopy() override;

    // Create a copy of just the material parameters
    // Called by the continuum elements
    NDMaterial *getCopy (const char *type) override;

    // Return a string indicating the type of material model
    virtual const char *getType() const;

    virtual int getOrder() const;
    
    virtual int sendSelf(int commitTag, Channel &);  
    virtual int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);    
    
    void Print(OPS_Stream &s, int flag = 0);

    virtual int setParameter(const char **argv, int argc, Parameter &);
    virtual int updateParameter(int parameterID, Information &);
    virtual int activateParameter(int paramID);

    friend class ElasticIsotropicBeamThread;

  protected:
    double E;	// Elastic modulus
    double v;	// Poisson ratio
    double rho ; //mass per unit 3D volume

    int parameterID;
  private:

};
} //

#endif
