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
//
#include <ElasticPlaneStress.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>

// Vector ElasticPlaneStress :: strain_vec(3) ;
Vector ElasticPlaneStress :: stress_vec(3) ;
Matrix ElasticPlaneStress :: tangent_matrix(3,3) ;


//null constructor
ElasticPlaneStress ::  ElasticPlaneStress( ) : 
NDMaterial(0, ND_TAG_ElasticPlaneStress), 
strain_vec(3),
E(0),
nu(0),
rho(0)
{
}


//full constructor
ElasticPlaneStress :: 
ElasticPlaneStress(int tag, 
                   double E_,
                   double nu_,
                   double rho_) : 
NDMaterial(tag, ND_TAG_ElasticPlaneStress), 
strain_vec(3),
E(E_),
nu(nu_),
rho(rho_)
{ 

}


//destructor
ElasticPlaneStress :: ~ElasticPlaneStress() 
{  } 


NDMaterial* ElasticPlaneStress :: getCopy( ) 
{ 
  ElasticPlaneStress  *clone;
  clone = new ElasticPlaneStress( ) ;   //new instance of this class
  *clone = *this ;                 //asignment to make copy
  return clone ;
}


//send back type of material
const char* ElasticPlaneStress :: getType( ) const 
{
  return "PlaneStress" ;
}


int ElasticPlaneStress :: getOrder( ) const 
{ 
  return 3 ; 
} 

//mass per unit volume
double
ElasticPlaneStress::getRho( )
{
  return rho ;
}

//get the strain and integrate plasticity equations
int ElasticPlaneStress :: setTrialStrain( const Vector &strain_from_element ) 
{
  strain_vec = strain_from_element;
  
  return 0 ;
}


// unused trial strain functions
int ElasticPlaneStress :: setTrialStrain( const Vector &v, const Vector &r )
{ 
   opserr << "ElasticPlaneStress :: setTrialStrain( const Vector &v, const Vector &r ) -- should not be used! \n";
   return this->setTrialStrain( v ) ;
} 

int ElasticPlaneStress :: setTrialStrainIncr( const Vector &v ) 
{
   opserr << "ElasticPlaneStress :: setTrialStrainIncr( const Vector &v ) -- should not be used! \n";
   return -1 ;
}

int ElasticPlaneStress :: setTrialStrainIncr( const Vector &v, const Vector &r ) 
{
   opserr << "ElasticPlaneStress :: setTrialStrainIncr( const Vector &v, const Vector &r ) -- should not be used! \n";
    return this->setTrialStrainIncr(v);
}

const Vector& ElasticPlaneStress :: getStrain( ) 
{
  return strain_vec ;
}


const Vector& ElasticPlaneStress :: getStress( ) 
{

  //stress_vec = this->getTangent() * strain_vec;
  double den = 1 - nu*nu;
  double G = E / (2*(1+nu));

  stress_vec(0) = E/den*(strain_vec(0) + nu*strain_vec(1));
  stress_vec(1) = E/den*(strain_vec(1) + nu*strain_vec(0));
  stress_vec(2) = G*strain_vec(2);
  
  return stress_vec ;
}

const Matrix& 
ElasticPlaneStress :: getTangent( ) 
{
  // matrix to tensor mapping
  //  Matrix      Tensor
  // -------     -------
  //   0           0 0
  //   1           1 1
  //   2           0 1  ( or 1 0 ) 
  // 
       
  double den = 1 - nu*nu;
  double G = E / (2*(1+nu));

  tangent_matrix(0,0) = E / den ;
  tangent_matrix(1,1) = E / den ;
  tangent_matrix(2,2) = G ;

  tangent_matrix(0,1) = nu * E / den ;
  tangent_matrix(1,0) = nu * E / den ;

  tangent_matrix(0,2) = 0. ;
  tangent_matrix(2,0) = 0. ;

  tangent_matrix(1,2) = 0. ;
  tangent_matrix(2,1) = 0. ;

  return tangent_matrix ;
} 


const Matrix& ElasticPlaneStress::getInitialTangent( ) 
{
  // matrix to tensor mapping
  //  Matrix      Tensor
  // -------     -------
  //   0           0 0
  //   1           1 1
  //   2           0 1  ( or 1 0 ) 
  //
  return this->getTangent() ;
} 

int 
ElasticPlaneStress::commitState( ) 
{
  return 0;
}

int 
ElasticPlaneStress::revertToLastCommit( ) 
{
  return 0;
}


int 
ElasticPlaneStress::revertToStart()
{
  return 0;
}

int
ElasticPlaneStress::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int
ElasticPlaneStress::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  return -1;
}


void ElasticPlaneStress :: Print( OPS_Stream &s, int flag )
{
  s << endln ;
  s << "ElasticPlaneStress : " ; 
  s << this->getType( ) << endln ;
  s << "Elastic Modulus =   " << E        << endln ;
  s << "Poisson's ratio =  " << nu       << endln ;
  s << "mass density =        " << rho     << endln ;
  s << endln ;
}
