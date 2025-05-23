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

#ifndef ElasticOrthotropicPlaneStress_h
#define ElasticOrthotropicPlaneStress_h

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>

#include <NDMaterial.h>

class ElasticOrthotropicPlaneStress : public NDMaterial {

  public : 

  //null constructor
  ElasticOrthotropicPlaneStress( ) ;

  //full constructor
  ElasticOrthotropicPlaneStress(int tag, 
                   double E1,
                   double E2,
                   double nu12,
                   double nu21,
                   double G12,
                   double rho) ;


  //destructor
  ~ElasticOrthotropicPlaneStress( ) ;

  const char *getClassType(void) const {return "ElasticOrthotropicPlaneStress";};

    NDMaterial* getCopy( ) ;

  //send back type of material
  const char* getType( ) const ;

    int getOrder( ) const ;

  //mass per unit volume
  double getRho();

  //get the strain and integrate plasticity equations
  int setTrialStrain( const Vector &strain_from_element) ;

  //unused trial strain functions
  int setTrialStrain( const Vector &v, const Vector &r ) ;
  int setTrialStrainIncr( const Vector &v ) ;
  int setTrialStrainIncr( const Vector &v, const Vector &r ) ;

  const Vector& getStrain( ) ;

  const Vector& getStress( ) ;

  const Matrix& getTangent( ) ;
  const Matrix& getInitialTangent( ) ;

  //swap history variables
  int commitState( ) ; 
  int revertToLastCommit( ) ;
  int revertToStart( ) ;

  //sending and receiving
  int sendSelf(int commitTag, Channel &theChannel) ;  
  int recvSelf(int commitTag, Channel &theChannel, 
               FEM_ObjectBroker &theBroker ) ;
  //print out material data
  void Print(OPS_Stream &s, int flag = 0) ;
  
  private : 
  
  //static vectors and matrices
  Vector strain_vec ;     //strain in vector notation
  static Vector stress_vec ;     //stress in vector notation
  static Matrix tangent_matrix ; //material tangent in matrix notation

  double E1, E2, nu12, nu21, G12, rho;

  //index mapping special for plane stress because of 
  // condensation on tangent
  void index_map( int matrix_index, int &i, int &j ) ;

} ; //end of ElasticOrthotropicPlaneStress declarations


#endif
