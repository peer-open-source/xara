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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-02-14 23:01:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/PlaneStressMaterial.cpp,v $

//
// Ed "C++" Love
//
// Generic Plane Stress Material
//
#include <PlaneStressMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Vector.h>
#include <Vector3D.h>
#include <Matrix3D.h>
using namespace OpenSees;

// static vector and matrices
Vector  PlaneStressMaterial::stress(3) ;
Matrix  PlaneStressMaterial::tangent(3,3) ;

//      0  1  2  3  4  5
// ND: 11 22 33 12 23 31
// PS: 11 22 12 33 23 31


PlaneStressMaterial::PlaneStressMaterial( ) : 
NDMaterial(0, ND_TAG_PlaneStressMaterial ), 
strain(3) 
{ }


PlaneStressMaterial::PlaneStressMaterial(    
				   int tag, 
                                   NDMaterial &the3DMaterial ) :
NDMaterial( tag, ND_TAG_PlaneStressMaterial ),
strain(3)
{
  theMaterial = the3DMaterial.getCopy("ThreeDimensional") ;

  Tstrain22 = 0.0 ;
  Tgamma02 = 0.0 ;
  Tgamma12 = 0.0 ;
  Cstrain22 = 0.0 ;
  Cgamma02 = 0.0 ;
  Cgamma12 = 0.0 ;
}



PlaneStressMaterial::~PlaneStressMaterial( ) 
{ 
  delete theMaterial ;
} 



NDMaterial*
PlaneStressMaterial::getCopy( ) 
{
  PlaneStressMaterial *clone ;   //new instance of this class

  clone = new PlaneStressMaterial( this->getTag(), 
                                   *theMaterial ) ; //make the copy

  clone->Tstrain22 = this->Tstrain22 ;
  clone->Tgamma02  = this->Tgamma02 ;
  clone->Tgamma12  = this->Tgamma12 ;
  clone->Cstrain22 = this->Cstrain22 ;
  clone->Cgamma02  = this->Cgamma02 ;
  clone->Cgamma12  = this->Cgamma12 ;

  return clone ;
}


NDMaterial* 
PlaneStressMaterial::getCopy( const char *type ) 
{
  return this->getCopy( ) ;
}

int
PlaneStressMaterial::getOrder() const
{
  return 3 ;
}


const char*
PlaneStressMaterial::getType( ) const 
{
  return "PlaneStress" ; 
}



int 
PlaneStressMaterial::commitState( ) 
{
  Cstrain22 = Tstrain22;
  Cgamma02 = Tgamma02;
  Cgamma12 = Tgamma12;

  return theMaterial->commitState( ) ;
}



int 
PlaneStressMaterial::revertToLastCommit( )
{
  Tstrain22 = Cstrain22;
  Tgamma02 = Cgamma02;
  Tgamma12 = Cgamma12;

  return theMaterial->revertToLastCommit( )  ;
}


int
PlaneStressMaterial::revertToStart( )
{
  this->Tstrain22 = 0.0 ;
  this->Tgamma12  = 0.0 ;
  this->Tgamma02  = 0.0 ;
  this->Cstrain22 = 0.0 ;
  this->Cgamma12  = 0.0 ;
  this->Cgamma02  = 0.0 ;
  
  strain.Zero();

  return theMaterial->revertToStart( ) ;
}


//mass per unit volume
double
PlaneStressMaterial::getRho( )
{
  return theMaterial->getRho( ) ;
}


//receive the strain
int 
PlaneStressMaterial::setTrialStrain( const Vector &strainFromElement )
{
  static const double tolerance = 1.0e-08 ;

  this->strain(0) = strainFromElement(0) ;
  this->strain(1) = strainFromElement(1) ;
  this->strain(2) = strainFromElement(2) ;

  double norm ;
  Vector3D condensedStress;
  Vector3D strainIncrement;
  Vector threeDstrain(6);
  Matrix3D dd22;

  int count = 0;
  const int maxCount = 20;
  double norm0;

  //newton loop to solve for out-of-plane strains
  do {

    //set three dimensional strain
    threeDstrain(0) = this->strain(0) ;
    threeDstrain(1) = this->strain(1) ;
    threeDstrain(2) = this->Tstrain22 ;
    threeDstrain(3) = this->strain(2) ; 
    threeDstrain(4) = this->Tgamma12 ;
    threeDstrain(5) = this->Tgamma02 ;

    if (theMaterial->setTrialStrain( threeDstrain ) < 0) {
      opserr << "PlaneStressMaterial::setTrialStrain() - setTrialStrain in material failed with strain " << threeDstrain;
      return -1;
    }

    const Vector &threeDstress = theMaterial->getStress();
    const Matrix &threeDtangent = theMaterial->getTangent();

    //NDmaterial strain order          = 11, 22, 33, 12, 23, 31 
    //PlaneStressMaterial strain order = 11, 22, 12, 33, 23, 31 

    condensedStress[0] = threeDstress(2);
    condensedStress[1] = threeDstress(4);
    condensedStress[2] = threeDstress(5);

    dd22(0,0) = threeDtangent(2,2);
    dd22(1,0) = threeDtangent(4,2);
    dd22(2,0) = threeDtangent(5,2);

    dd22(0,1) = threeDtangent(2,4);
    dd22(1,1) = threeDtangent(4,4);
    dd22(2,1) = threeDtangent(5,4);

    dd22(0,2) = threeDtangent(2,5);
    dd22(1,2) = threeDtangent(4,5);
    dd22(2,2) = threeDtangent(5,5);

    //set norm
    norm = condensedStress.norm();
    if (count == 0)
      norm0 = norm;

    //condensation 
    dd22.solve(condensedStress, strainIncrement);

    // Update
    this->Tstrain22 -= strainIncrement[0];
    this->Tgamma12  -= strainIncrement[1];
    this->Tgamma02  -= strainIncrement[2];

  } while (count++ < maxCount && norm > tolerance);

  return 0;
}


const Vector& 
PlaneStressMaterial::getStrain( )
{
  return strain;
}


const Vector&  
PlaneStressMaterial::getStress( )
{
  //three dimensional stress
  const Vector &threeDstress = theMaterial->getStress();
  
  stress(0) = threeDstress(0);
  stress(1) = threeDstress(1);
  stress(2) = threeDstress(3);

  return stress;
}

const Vector& 
PlaneStressMaterial::getStressSensitivity(int gradIndex,
					  bool conditional)
{
  const Vector &threeDstress = theMaterial->getStressSensitivity(gradIndex, conditional);

  stress(0) = threeDstress(0);
  stress(1) = threeDstress(1);
  stress(2) = threeDstress(3);

  const Matrix &threeDtangent = theMaterial->getTangent();

  static Matrix dd12(3,3);
  static Matrix dd22(3,3);
  //
  dd12(0,0) = threeDtangent(0,2);
  dd12(1,0) = threeDtangent(1,2);
  dd12(2,0) = threeDtangent(3,2);

  dd12(0,1) = threeDtangent(0,4);
  dd12(1,1) = threeDtangent(1,4);
  dd12(2,1) = threeDtangent(3,4);

  dd12(0,2) = threeDtangent(0,5);
  dd12(1,2) = threeDtangent(1,5);
  dd12(2,2) = threeDtangent(3,5);

  //
  dd22(0,0) = threeDtangent(2,2);
  dd22(1,0) = threeDtangent(4,2);
  dd22(2,0) = threeDtangent(5,2);
  
  dd22(0,1) = threeDtangent(2,4);
  dd22(1,1) = threeDtangent(4,4);
  dd22(2,1) = threeDtangent(5,4);
  
  dd22(0,2) = threeDtangent(2,5);
  dd22(1,2) = threeDtangent(4,5);
  dd22(2,2) = threeDtangent(5,5);

  
  static Vector sigma2(3);
  sigma2(0) = threeDstress(2);
  sigma2(1) = threeDstress(4);
  sigma2(2) = threeDstress(5);

  static Vector dd22sigma2(3);
  dd22.Solve(sigma2,dd22sigma2);

  stress.addMatrixVector(1.0, dd12, dd22sigma2, -1.0);

  return stress;
}

const Matrix&  
PlaneStressMaterial::getTangent( )
{
  const Matrix &C = theMaterial->getTangent();

//Matrix3D dd11, dd12, dd21, dd22;
  static Matrix dd11(3,3);
  static Matrix dd12(3,3);
  static Matrix dd21(3,3);
  static Matrix dd22(3,3);
  //
  dd11(0,0) = C(0,0);
  dd11(1,0) = C(1,0);
  dd11(2,0) = C(3,0);

  dd11(0,1) = C(0,1);
  dd11(1,1) = C(1,1);
  dd11(2,1) = C(3,1);

  dd11(0,2) = C(0,3);
  dd11(1,2) = C(1,3);
  dd11(2,2) = C(3,3);

  //
  dd12(0,0) = C(0,2);
  dd12(1,0) = C(1,2);
  dd12(2,0) = C(3,2);

  dd12(0,1) = C(0,4);
  dd12(1,1) = C(1,4);
  dd12(2,1) = C(3,4);

  dd12(0,2) = C(0,5);
  dd12(1,2) = C(1,5);
  dd12(2,2) = C(3,5);

  //
  dd21(0,0) = C(2,0);
  dd21(1,0) = C(4,0);
  dd21(2,0) = C(5,0);

  dd21(0,1) = C(2,1);
  dd21(1,1) = C(4,1);
  dd21(2,1) = C(5,1);

  dd21(0,2) = C(2,3);
  dd21(1,2) = C(4,3);
  dd21(2,2) = C(5,3);

  //
  dd22(0,0) = C(2,2);
  dd22(1,0) = C(4,2);
  dd22(2,0) = C(5,2);

  dd22(0,1) = C(2,4);
  dd22(1,1) = C(4,4);
  dd22(2,1) = C(5,4);

  dd22(0,2) = C(2,5);
  dd22(1,2) = C(4,5);
  dd22(2,2) = C(5,5);


  //int Solve(const Vector &V, Vector &res) const;
  //int Solve(const Matrix &M, Matrix &res) const;
  //condensation 
  static Matrix dd22invdd21(3,3);
  dd22.Solve(dd21, dd22invdd21);

  //this->tangent   = dd11 ; 
  //this->tangent  -= ( dd12*dd22invdd21 ) ;

  dd11.addMatrixProduct(1.0, dd12, dd22invdd21, -1.0);
  tangent = dd11;  

  return tangent;
}



const Matrix&  
PlaneStressMaterial::getInitialTangent( )
{
  const Matrix &threeDtangent = theMaterial->getInitialTangent();

  static Matrix dd11(3,3);
  static Matrix dd12(3,3);
  static Matrix dd21(3,3);
  static Matrix dd22(3,3);

  dd11(0,0) = threeDtangent(0,0);
  dd11(1,0) = threeDtangent(1,0);
  dd11(2,0) = threeDtangent(3,0);

  dd11(0,1) = threeDtangent(0,1);
  dd11(1,1) = threeDtangent(1,1);
  dd11(2,1) = threeDtangent(3,1);

  dd11(0,2) = threeDtangent(0,3);
  dd11(1,2) = threeDtangent(1,3);
  dd11(2,2) = threeDtangent(3,3);

  //
  dd12(0,0) = threeDtangent(0,2);
  dd12(1,0) = threeDtangent(1,2);
  dd12(2,0) = threeDtangent(3,2);

  dd12(0,1) = threeDtangent(0,4);
  dd12(1,1) = threeDtangent(1,4);
  dd12(2,1) = threeDtangent(3,4);

  dd12(0,2) = threeDtangent(0,5);
  dd12(1,2) = threeDtangent(1,5);
  dd12(2,2) = threeDtangent(3,5);

  //
  dd21(0,0) = threeDtangent(2,0);
  dd21(1,0) = threeDtangent(4,0);
  dd21(2,0) = threeDtangent(5,0);

  dd21(0,1) = threeDtangent(2,1);
  dd21(1,1) = threeDtangent(4,1);
  dd21(2,1) = threeDtangent(5,1);

  dd21(0,2) = threeDtangent(2,3);
  dd21(1,2) = threeDtangent(4,3);
  dd21(2,2) = threeDtangent(5,3);

  //
  dd22(0,0) = threeDtangent(2,2);
  dd22(1,0) = threeDtangent(4,2);
  dd22(2,0) = threeDtangent(5,2);

  dd22(0,1) = threeDtangent(2,4);
  dd22(1,1) = threeDtangent(4,4);
  dd22(2,1) = threeDtangent(5,4);

  dd22(0,2) = threeDtangent(2,5);
  dd22(1,2) = threeDtangent(4,5);
  dd22(2,2) = threeDtangent(5,5);


  // condensation 
  static Matrix dd22invdd21(3,3);
  dd22.Solve(dd21, dd22invdd21);

  //this->tangent   = dd11 ; 
  //this->tangent  -= ( dd12*dd22invdd21 ) ;
  dd11.addMatrixProduct(1.0, dd12, dd22invdd21, -1.0);
  tangent = dd11;

  return tangent;
}


void  
PlaneStressMaterial::Print( OPS_Stream &s, int flag )
{
  if (flag ==  OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_ELEM_INDENT << "{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"PlaneStressMaterial\", ";
    s << "\"material\": \"" << theMaterial->getTag() << "\"";
    s << "}";

  } else if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
    s << "General Plane Stress Material \n" ;
    s << " Tag: " << this->getTag() << "\n" ; 
    s << "using the 3D material : \n" ;

    theMaterial->Print( s, flag ) ;
  }

  return ;
}


int 
PlaneStressMaterial::sendSelf(int commitTag, Channel &theChannel) 
{
  int res = 0;

  // put tag and associated materials class and database tags into an id and send it
  static ID idData(3);
  idData(0) = this->getTag();
  idData(1) = theMaterial->getClassTag();
  int matDbTag = theMaterial->getDbTag();
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    theMaterial->setDbTag(matDbTag);
  }
  idData(2) = matDbTag;

  res = theChannel.sendID(this->getDbTag(), commitTag, idData);
  if (res < 0) {
    opserr << "PlaneStressMaterial::sendSelf() - failed to send id data\n";
    return res;
  }

  // put the strains in a vector and send it
  static Vector vecData(3);
  vecData(0) = Cstrain22;
  vecData(1) = Cgamma02;
  vecData(2) = Cgamma12;

  res = theChannel.sendVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) {
    opserr << "PlaneStressMaterial::sendSelf() - failed to send vector data\n";
    return res;
  }

  // now send the materials data
  res = theMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) 
    opserr << "PlaneStressMaterial::sendSelf() - failed to send vector material\n";

  return res;
}

int 
PlaneStressMaterial::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;

  // recv an id containing the tag and associated materials class and db tags
  static ID idData(3);
  res = theChannel.recvID(this->getDbTag(), commitTag, idData);
  if (res < 0) {
    opserr << "PlaneStressMaterial::sendSelf() - failed to send id data\n";
    return res;
  }

  this->setTag(idData(0));
  int matClassTag = idData(1);
  
  // if the associated material has not yet been created or is of the wrong type
  // create a new material for recvSelf later
  if (theMaterial == 0 || theMaterial->getClassTag() != matClassTag) {
    if (theMaterial != 0)
      delete theMaterial;
    theMaterial = theBroker.getNewNDMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << "PlaneStressMaterial::recvSelf() - failed to get a material of type: " << matClassTag << endln;
      return -1;
    }
  }
  theMaterial->setDbTag(idData(2));

  // recv a vector containing strains and set the strains
  static Vector vecData(3);
  res = theChannel.recvVector(this->getDbTag(), commitTag, vecData);
  if (res < 0) {
    opserr << "PlaneStressMaterial::sendSelf() - failed to send vector data\n";
    return res;
  }

  Cstrain22 = vecData(0);
  Cgamma02 = vecData(1);
  Cgamma12  = vecData(2);

  Tstrain22 = Cstrain22;
  Tgamma02 = Cgamma02;
  Tgamma12  = Cgamma12;

  // now receive the materials data
  res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) 
    opserr << "PlaneStressMaterial::sendSelf() - failed to send vector material\n";
  
  return res;
}
 
int
PlaneStressMaterial::setParameter(const char **argv, int argc,
				  Parameter &param)
{
  return theMaterial->setParameter(argv, argc, param);
}

Response* PlaneStressMaterial::setResponse(const char** argv, int argc, OPS_Stream& s)
{
    // for strain, stress and tangent use the base class implementation
        // so that the output will be that of the adapter
    if (strcmp(argv[0], "Tangent") == 0 ||
        strcmp(argv[0], "tangent") == 0 ||
        strcmp(argv[0], "stress") == 0 ||
        strcmp(argv[0], "stresses") == 0 ||
        strcmp(argv[0], "strain") == 0 ||
        strcmp(argv[0], "strains") == 0
        ) {
        return NDMaterial::setResponse(argv, argc, s);
    }
    // otherwise, for other custom results, forward the call to the adaptee
    return theMaterial->setResponse(argv, argc, s);
}
