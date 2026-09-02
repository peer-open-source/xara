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
// Written: Ed "C++" Love
//
// J2PlaneStrain isotropic hardening material class
// 
//  Elastic Model
//  sigma = K*trace(epsilion_elastic) + (2*G)*dev(epsilon_elastic)
//
//  Yield Function
//  phi(sigma,q) = || dev(sigma) ||  - sqrt(2/3)*q(xi) 
//
//  Saturation Isotropic Hardening with linear term
//  q(xi) = simga_0 + (sigma_infty - sigma_0)*exp(-delta*xi) + H*xi 
//
//  Flow Rules
//  \dot{epsilon_p} =  gamma * d_phi/d_sigma
//  \dot{xi}        = -gamma * d_phi/d_q 
//
//  Linear Viscosity 
//  gamma = phi / eta  ( if phi > 0 ) 
//
//  Backward Euler Integration Routine 
//  Yield condition enforced at time n+1 
//
//  Send strains in following format :
// 
//     strain_vec = {   eps_00
//                      eps_11
//                    2 eps_01   }   <--- note the 2
// 
//  set eta := 0 for rate independent case
//
#ifndef J2PlaneStrain_h
#define J2PlaneStrain_h

#include <math.h> 

#include <Vector.h>
#include <Matrix.h>

#include <J2Plasticity.h>
using namespace OpenSees;

class J2PlaneStrain : public J2Plasticity {
  public : 

  J2PlaneStrain( ) ;

  //full constructor
  J2PlaneStrain(   int    tag, 
                   double K,
                   double G,
                   double yield0,
                   double yield_infty,
                   double d,
                   double H,
                   double viscosity = 0,
		   double rho = 0) ;


  ~J2PlaneStrain();

  NDMaterial* getCopy() final;
  const char* getType( ) const ;

  const char *getClassType() const final {return "J2PlaneStrain";}

  int getOrder( ) const ;


  //get the strain and integrate plasticity equations
  int setTrialStrain( const Vector &strain_from_element) ;

  // unused trial strain functions
  int setTrialStrainIncr(const Vector &v);

  const Vector& getStrain( );
  const Vector& getStress( );
  const Matrix& getTangent( );
  const Matrix& getInitialTangent( );

  int commitState(); 
  int revertToLastCommit( );
  int revertToStart();


protected:
  // matrix index to tensor index mapping
  void index_map(int matrix_index, int& i, int& j) const final;

private :
    
  // static vectors and matrices
  static Vector strain_vec ;     //strain in vector notation
  static Vector stress_vec ;     //stress in vector notation
  static Matrix tangent_matrix ; //material tangent in matrix notation
};

#endif
