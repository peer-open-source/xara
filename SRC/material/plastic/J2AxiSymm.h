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
                                                                        
// $Revision: 1.7 $
// $Date: 2008-10-20 22:23:03 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/J2AxiSymm.h,v $

// Written: Ed "C++" Love
//
// J2AxiSymmetric isotropic hardening material class
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
//	                eps_22 		      
//                    2 eps_01   }   <--- note the 2
// 
//  set eta := 0 for rate independent case
//

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>

#include <J2Plasticity.h>
using namespace OpenSees;


class J2AxiSymm : public J2Plasticity {

public:

  J2AxiSymm();
  J2AxiSymm(int    tag, 
            double K,
            double G,
            double yield0,
            double yield_infty,
            double d,
            double H,
            double viscosity = 0,
            double rho = 0 );

  ~J2AxiSymm();

  const char *getClassType() const override {return "J2AxiSymm";}

  NDMaterial* getCopy( ) override;
  const char* getType( ) const override;

  int getOrder( ) const override;

  //get the strain and integrate plasticity equations
  int setTrialStrain( const Vector &strain_from_element) override;

  //unused trial strain functions
  int setTrialStrain( const Vector &v, const Vector &r ) ;
  int setTrialStrainIncr( const Vector &v ) ;
  int setTrialStrainIncr( const Vector &v, const Vector &r ) ;

  const Vector& getStrain() override;
  const Vector& getStress() override;
  const Matrix& getTangent() override;
  const Matrix& getInitialTangent() override;

  // swap history variables
  int commitState() override;
  int revertToLastCommit() override;
  int revertToStart() override;

  // sending and receiving
  int sendSelf(int commitTag, Channel &) override;
  int recvSelf(int commitTag, Channel &, FEM_ObjectBroker &) override;


protected:
  // matrix index to tensor index mapping
  void index_map(int matrix_index, int& i, int& j) const final;

private:

  // static vectors and matrices
  static Vector strain_vec ;     //strain in vector notation
  static Vector stress_vec ;     //stress in vector notation
  static Matrix tangent_matrix ; //material tangent in matrix notation

  double commitEps00;
  double commitEps11;
  double commitEps01;
  double commitEps22;
 				     
} ; //end of J2AxiSymm declarations


