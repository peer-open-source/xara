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
                                                                        
// $Revision: 1.9 $
// $Date: 2003-02-14 23:01:33 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/ElasticMembranePlateSection.cpp,v $

// Ed "C++" Love
//
//  Elastic Plate Section with membrane
//


#include "OrthotropicMembraneSection.h"

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>

#include <ID.h>
#include <Domain.h>
#include <string.h>
#include <Information.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Parameter.h>
#include <math.h>
#include <algorithm>
#include <elementAPI.h>
#include <string>



//parameters
const double OrthotropicMembraneSection::five6 = 5.0/6.0 ; //shear correction

//static vector and matrices
Vector  OrthotropicMembraneSection::stress(8) ;
Matrix  OrthotropicMembraneSection::tangent(8,8) ;
ID      OrthotropicMembraneSection::array(8) ;


//null constructor
OrthotropicMembraneSection::OrthotropicMembraneSection( ) : SectionForceDeformation(0, 0), strain(8)  { 
}

//full constructor
OrthotropicMembraneSection::OrthotropicMembraneSection(int tag, double _E1,double _E2, double _v, double _G, double _thickness, double _r, double _angle):
	SectionForceDeformation(tag, 0), strain(8), toLocal(3,3), toMaterial(3,3)
{
  this->E1   = _E1;
  this->E2   = _E2;
  this->v    = _v;
  this->G    = _G;
  this->h    = _thickness;
  this->rhoH = _r*_thickness;
  this->angle = _angle; // in degrees

  double c = cos(_angle/180.* 3.14159);
  double s = sin(_angle / 180.* 3.14159);

  toMaterial(0, 0) = c*c;
  toMaterial(0, 1) = s*s;
  toMaterial(0, 2) = s*c;

  toMaterial(1, 0) = s*s;
  toMaterial(1, 1) = c*c;
  toMaterial(1, 2) = -s*c;

  toMaterial(2, 0) = -2.0*s*c;
  toMaterial(2, 1) = 2.0*s*c;
  toMaterial(2, 2) = c*c-s*s;


  c = cos(-_angle / 180.* 3.14159);
  s = sin(-_angle / 180.* 3.14159);

  toLocal(0, 0) = c*c;
  toLocal(0, 1) = s*s;
  toLocal(0, 2) = 2.0*s*c;

  toLocal(1, 0) = s*s;
  toLocal(1, 1) = c*c;
  toLocal(1, 2) = -2.0*s*c;

  toLocal(2, 0) = -s*c;
  toLocal(2, 1) = s*c;
  toLocal(2, 2) = c*c - s*s;

}

//destructor
OrthotropicMembraneSection::~OrthotropicMembraneSection() { 
} 

//make a clone of this material
SectionForceDeformation*  OrthotropicMembraneSection::getCopy( ) 
{
  OrthotropicMembraneSection *clone ;   

  clone = new OrthotropicMembraneSection(this->getTag(), E1, E2, v, G, h, rhoH/h, angle) ; //new instance of this class

  clone->rhoH = this->rhoH ;
  clone->strain = this->strain;

  return clone ;
}

//density per unit area
double
OrthotropicMembraneSection::getRho() {
  return rhoH ;
}

//send back order of strain in vector form
int 
OrthotropicMembraneSection::getOrder() const {
  return 8 ;
}

//send back order of strain in vector form
const ID& OrthotropicMembraneSection::getType() {
  return array ;
}

int 
OrthotropicMembraneSection::commitState()  {
  return 0 ;
}

int 
OrthotropicMembraneSection::revertToLastCommit() {
  return 0 ;
}

int 
OrthotropicMembraneSection::revertToStart() {
  return 0 ;
}

int 
OrthotropicMembraneSection::setTrialSectionDeformation(const Vector &strain_from_element) {
  this->strain = strain_from_element;
  return 0 ;
}


const Vector& 
OrthotropicMembraneSection::getSectionDeformation() {
  return this->strain ;
}

const Vector&  
OrthotropicMembraneSection::getStressResultant() {

  // membrane formulation from: "TREMURI program: An equivalent frame model for the nonlinear seismic analysis of masonry buildings",
  // Sergio Lagomarsino, Andrea Penna, Alessandro Galasco, Serena Cattari, Engineering Structures 56, 2013

  // For shell 3d elements, the first material axis has the direction of the mean of the orientation of the side 1-2 and 4-3 (defined counterclockwise).
  // The strains in local coordinates are transformed in material coordinates and the brought back to the local system

  Vector strainMat(3);
  Vector stressMat(3);
  Vector temp(3);

  for (int i = 0; i < 3; i++) {
	  temp(i) = strain(i);
  }

  strainMat = toMaterial*temp;

  double e = E2/E1;
  double stiffness1 = E1 / (1.0 - e*v*v) *h;

  //membrane resultants
  stressMat(0) =     stiffness1 * strainMat(0)   + v*e*stiffness1* strainMat(1);
  stressMat(1) = v*e*stiffness1 * strainMat(0)   +   e*stiffness1* strainMat(1);
  stressMat(2) =  G*h*strainMat(2) ;

  double GG = h*G*five6 ;
  double D  =  E1 * (h*h*h) / 12.0 / ( 1.0 - v*v ) ;  //bending modulus

  //bending resultants
  stress(3) = -( D*strain(3) + v*D*strain(4) ) ;
  stress(4) = -( v*D*strain(3) + D*strain(4) ) ;
  stress(5) = -0.5*D*( 1.0 - v )*strain(5) ;
  stress(6) = GG*strain(6) ;
  stress(7) = GG*strain(7) ;

  
  // zero bending components:
  for (int i=3; i<8; i++) {
	  stress(i) *= 0.0;
  }
  
  temp = toLocal*stressMat;
  for (int i = 0; i < 3; i++) {
	  stress(i) = temp(i);
  }

   //opserr << "stressUpdate = " << stress;
 
  return this->stress ;
}


//send back the tangent 
const Matrix&  
OrthotropicMembraneSection::getSectionTangent() {

  double e = E2/E1;
  double stiffness1 = E1 / (1.0 - e*v*v) *h;

  Matrix tangentMat(3, 3);
  Matrix tangentLoc(3, 3);

  tangent.Zero() ;

  //membrane tangent terms
  tangentMat(0,0) = stiffness1;
  tangentMat(0,1) = v*e*stiffness1;
  tangentMat(1,0) = tangentMat(0,1);
  tangentMat(1,1) = e*stiffness1;
  
  tangentMat(2,2) = G*h ;
  
  double GG = h*G*five6 ;
  double D  =  E1 * (h*h*h) / 12.0 / ( 1.0 - v*v ) ;  //bending modulus
  
  //bending tangent terms
  tangent(3,3) = -D ; //  *(0.0001);
  tangent(4,4) = -D ; //  *(0.0001);

  tangent(3,4) = -v*D ; //  *(0.0001);
  tangent(4,3) = tangent(3,4) ;

  tangent(5,5) = -0.5 * D * ( 1.0 - v ) ; //   *(0.0001);

  tangent(6,6) = GG ; //  *(0.0001);

  tangent(7,7) = GG ; //  *(0.0001);
  

  
  // zero bending components:
  for (int i=3; i<8; i++) {
	  tangent(i,i) = 0.0;
  }
  
  tangentLoc = toLocal*tangentMat*toMaterial;
  for (int i = 0; i<3; i++) {
	  for (int j = 0; j < 3; j++) {
		  tangent(i, j) = tangentLoc(i,j);
	  }
  }

  return this->tangent ;
}


const Matrix&  
OrthotropicMembraneSection::getInitialTangent() {
  return this->getSectionTangent();
}


//print out data
void  
OrthotropicMembraneSection::Print( OPS_Stream &s, int flag ) {
  s << "OrthotropicMembraneSection: \n " ;
  s <<  "  Young's Modulus E1 = "  <<  E1  <<  endln ;
  s <<  "  Young's Modulus E2 = "  <<  E2  <<  endln ;
  s <<  "  Poisson's Ratio nu = "  <<  v   <<  endln ;
  s <<  "  Shear Modulus = "       <<  G   <<  endln ;
  s <<  "  Thickness h = "         <<  h   <<  endln ;
  s <<  "  Density rho = "         <<  (rhoH/h)  <<  endln ;

  return ;
}

int 
OrthotropicMembraneSection::sendSelf(int cTag, Channel &theChannel) {
  int res = 0;
  static Vector data(7);
  data(0) = this->getTag();
  data(1) = E1;
  data(2) = E2;
  data(3) = v;
  data(4) = G;
  data(5) = h;
  data(6) = rhoH;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ElasticMembranePlateSection::sendSelf() - failed to send data\n";

  return res;
}


int 
OrthotropicMembraneSection::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker) {
  int res = 0;
  static Vector data(7);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ElasticMembranePlateSection::recvSelf() - failed to recv data\n";
  else {
    this->setTag(data(0));
    E1    = data(1);
	E2    = data(2);
    v     = data(3);
	G     = data(4);
    h     = data(5);
    rhoH  = data(6);
  }

  return res;
}
 
