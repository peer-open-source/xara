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
// $Revision: 1.2 $
// $Date: 2009-05-20 22:03:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/PlaneStrainMaterial.cpp,v $
//
// Antonios Vytiniotis
//
// Generic Plane Strain Material
//
#include <VectorND.h>
#include <PlaneStrainMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MaterialResponse.h>
using namespace OpenSees;

Vector  PlaneStrainMaterial::stress(3) ;
Matrix  PlaneStrainMaterial::tangent(3,3) ;


PlaneStrainMaterial::PlaneStrainMaterial( ) : 
NDMaterial(0, ND_TAG_PlaneStrainMaterial ), 
strain(3) 
{
}



PlaneStrainMaterial::PlaneStrainMaterial(    
				   int tag, NDMaterial &the3DMaterial ) :
NDMaterial( tag, ND_TAG_PlaneStrainMaterial ),
strain(3)
{
  theMaterial = the3DMaterial.getCopy("ThreeDimensional");
  if (theMaterial == 0) {
    theMaterial = the3DMaterial.getCopy( ) ;
  }
}


//destructor
PlaneStrainMaterial::~PlaneStrainMaterial( ) 
{ 
  delete theMaterial ;
} 



NDMaterial*
PlaneStrainMaterial::getCopy( ) 
{
  PlaneStrainMaterial *clone ;   //new instance of this class

  clone = new PlaneStrainMaterial( this->getTag(), 
                                   *theMaterial ) ; //make the copy

  return clone ;
}


NDMaterial* 
PlaneStrainMaterial::getCopy( const char *type ) 
{
  return this->getCopy( ) ;
}


int 
PlaneStrainMaterial::getOrder( ) const
{
  return 3 ;
}


const char*
PlaneStrainMaterial::getType( ) const 
{
  return "PlaneStrain" ; 
}



//swap history variables
int 
PlaneStrainMaterial::commitState( ) 
{
  return theMaterial->commitState( ) ;
}



int 
PlaneStrainMaterial::revertToLastCommit( )
{
  return theMaterial->revertToLastCommit( )  ;
}


int
PlaneStrainMaterial::revertToStart( )
{

  strain.Zero();

  return theMaterial->revertToStart( ) ;
}


//mass per unit volume
double
PlaneStrainMaterial::getRho( )
{
  return theMaterial->getRho( ) ;
}


//receive the strain
int 
PlaneStrainMaterial::setTrialStrain( const Vector &strainFromElement )
{
  this->strain(0) = strainFromElement(0) ;
  this->strain(1) = strainFromElement(1) ;
  this->strain(2) = strainFromElement(2) ;

//static Vector threeDstrain(6) ;
  VectorND<6> threeDstrain;

  //set three dimensional strain
  threeDstrain[0] = this->strain(0) ;
  threeDstrain[1] = this->strain(1) ;
  threeDstrain[2] = 0.0 ;
  threeDstrain[3] = this->strain(2) ;
  threeDstrain[4] = 0.0 ;
  threeDstrain[5] = 0.0 ;

  if (theMaterial->setTrialStrain( threeDstrain ) < 0) {
    opserr << "PlaneStrainMaterial::setTrialStrain() - setTrialStrain in material failed with strain " ; // << threeDstrain;
    return -1;
  }

  return 0;
}


const Vector& 
PlaneStrainMaterial::getStrain( )
{
  return this->strain ;
}


const Vector&  
PlaneStrainMaterial::getStress( )
{
  //three dimensional stress
  const Vector &threeDstress = theMaterial->getStress();

  stress(0) = threeDstress(0);
  stress(1) = threeDstress(1);
  stress(2) = threeDstress(3);

  return this->stress ;
}


const Matrix&  
PlaneStrainMaterial::getTangent( )
{

//static Matrix threeDtangentCopy(6,6);

  //three dimensional tangent 
  const Matrix &threeDtangent = theMaterial->getTangent( ) ;

  tangent(0,0) = threeDtangent(0,0);
  tangent(1,0) = threeDtangent(1,0);
  tangent(2,0) = threeDtangent(3,0);
  tangent(0,1) = threeDtangent(0,1);
  tangent(1,1) = threeDtangent(1,1);
  tangent(2,1) = threeDtangent(3,1);
  tangent(0,2) = threeDtangent(0,3);
  tangent(1,2) = threeDtangent(1,3);
  tangent(2,2) = threeDtangent(3,3);

  return this->tangent ;
}

// AV not sure if it actually works
const Matrix&  
PlaneStrainMaterial::getInitialTangent( )
{
//static Matrix dd11(3,3) ;

//static Matrix threeDtangentCopy(6,6);

  //three dimensional tangent 
  const Matrix &threeDtangent = theMaterial->getInitialTangent( ) ;

  tangent(0,0) = threeDtangent(0,0);
  tangent(1,0) = threeDtangent(1,0);
  tangent(2,0) = threeDtangent(3,0);
  tangent(0,1) = threeDtangent(0,1);
  tangent(1,1) = threeDtangent(1,1);
  tangent(2,1) = threeDtangent(3,1);
  tangent(0,2) = threeDtangent(0,3);
  tangent(1,2) = threeDtangent(1,3);
  tangent(2,2) = threeDtangent(3,3);

  return this->tangent ;
}
// End of un-known part


// print out data
void  
PlaneStrainMaterial::Print(OPS_Stream &s, int flag)
{
  if (flag ==  OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_ELEM_INDENT << "{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"PlaneStrainMaterial\", ";
    s << "\"material\": \"" << theMaterial->getTag() << "\"";
    s << "}";

  } else if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
    s << "General Plane Strain Material \n" ;
    s << " Tag: " << this->getTag() << "\n" ; 
    s << "using the 3D material : \n" ;

    theMaterial->Print( s, flag ) ;
  }

  return ;
}


int 
PlaneStrainMaterial::sendSelf(int commitTag, Channel &theChannel) 
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
    opserr << "PlaneStrainMaterial::sendSelf() - failed to send id data\n";
    return res;
  }

  // now send the materials data
  res = theMaterial->sendSelf(commitTag, theChannel);
  if (res < 0) 
    opserr << "PlaneStrainMaterial::sendSelf() - failed to send vector material\n";

  return res;
}

int 
PlaneStrainMaterial::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;

  // recv an id containing the tag and associated materials class and db tags
  static ID idData(3);
  res = theChannel.recvID(this->getDbTag(), commitTag, idData);
  if (res < 0) {
    opserr << "PlaneStrainMaterial::sendSelf() - failed to send id data\n";
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
      opserr << "PlaneStrainMaterial::recvSelf() - failed to get a material of type: " << matClassTag << endln;
      return -1;
    }
  }
  theMaterial->setDbTag(idData(2));

  // now receive the materials data
  res = theMaterial->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) 
    opserr << "PlaneStrainMaterial::sendSelf() - failed to send vector material\n";
  
  return res;
}
 
int
PlaneStrainMaterial::setParameter(const char **argv, int argc,
				  Parameter &param)
{
  return theMaterial->setParameter(argv, argc, param);
}

Response* PlaneStrainMaterial::setResponse(const char** argv, int argc, OPS_Stream& s)
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
