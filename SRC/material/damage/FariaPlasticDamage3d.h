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
#ifndef PlasticDamageConcrete3d_h
#define PlasticDamageConcrete3d_h

// Written: Thanh Do
// Created: 07/16
//
// Description: 
//
// What: "@(#) ElasticIsotropicThreeDimesnional.h, revA"

#include <ElasticIsotropicMaterial.h>

#include <Matrix.h>
#include <MatrixND.h>
#include <Vector.h>
#include <ID.h>

class PlasticDamageConcrete3d : public NDMaterial
{
  public:
  PlasticDamageConcrete3d(int tag, 
                          double E, 
                          double nu, 
                          double ft,
                          double fc, 
                          double beta = 0.6, 
                          double Ap = 0.5, 
                          double An = 2.0, 
                          double Bn = 0.75);
    PlasticDamageConcrete3d();
    ~PlasticDamageConcrete3d();

    const char *getClassType() const {return "PlasticDamageConcrete3d";}

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

    NDMaterial*getCopy(const char *type);
    NDMaterial *getCopy();
    const char *getType() const;
    int getOrder() const;

    int sendSelf(int commitTag, Channel &);  
    int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &);    
    void Print(OPS_Stream &s, int flag =0);       

  protected:

  private:
    // parameters
    double E;     // elastic modulus
    double nu;    // Poisson ratio 
    double ft;    // tensile yield strength
    double fc;    // compressive yield strength
    double beta;  // plastic deformation rate
    double Ap;    // damage parameter
    double An;    // damage parameter
    double Bn;    // damage parameter

    // current state variables
    double rp;    // positive damage threshold
    double rn;    // negative damage threshold
    double dp;    // positive damage variable
    double dn;    // negative damage variable

    OpenSees::VectorND<6> eps;   // strain
    OpenSees::VectorND<6> sig;   // stress
    OpenSees::VectorND<6> sige;  // effective stress
    OpenSees::VectorND<6> eps_p; // plastic strain
    // OpenSees::VectorND<6> sigeP; // effective stress

    // committed state variables
    double rpCommit; 
    double rnCommit; 
    double dpCommit; 
    double dnCommit; 

    OpenSees::VectorND<6> epsCommit;
    OpenSees::VectorND<6> sigCommit;
    OpenSees::VectorND<6> sigeCommit;  
    OpenSees::VectorND<6> eps_pCommit; 
    // Vector sigePCommit;

    // tangent matrices
    OpenSees::MatrixND<6,6> Ce, C, Ccommit; 
};

#endif
