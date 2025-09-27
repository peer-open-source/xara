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
** ****************************************************************** */

#ifndef ElasticPlaneStress_h
#define ElasticPlaneStress_h

#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>

#include <NDMaterial.h>

class ElasticPlaneStress : public NDMaterial {

  public : 
  ElasticPlaneStress();
  ElasticPlaneStress(int tag, 
                   double E,
                   double nu,
                   double rho);
  ~ElasticPlaneStress();

  const char *getClassType() const {return "ElasticPlaneStress";}

    NDMaterial* getCopy( ) ;

  const char* getType( ) const override;
  int getOrder( ) const override;
  double getRho() override;

  int setTrialStrain( const Vector &strain_from_element) override;

  // unsupported trial strain functions
  int setTrialStrainIncr( const Vector &v );

  const Vector& getStrain() override;
  const Vector& getStress() override;
  const Matrix& getTangent() override;
  const Matrix& getInitialTangent() override;

  int commitState() override; 
  int revertToLastCommit() override;
  int revertToStart() override;

  int sendSelf(int commitTag, Channel &) override;  
  int recvSelf(int commitTag, Channel &, FEM_ObjectBroker & ) override;

  void Print(OPS_Stream &s, int flag) override;

private: 

  //static vectors and matrices
  Vector strain_vec ;     //strain in vector notation
  static Vector stress_vec ;     //stress in vector notation
  static Matrix tangent_matrix ; //material tangent in matrix notation

  double E, nu, rho;
};

#endif
