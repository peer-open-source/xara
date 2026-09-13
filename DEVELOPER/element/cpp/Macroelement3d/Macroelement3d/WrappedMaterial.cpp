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

/*
Source: /DEVELOPER/element/cpp/Macroelement3d/Macroelement3d/WrappedMaterial.cpp
Written by Francesco Vanin (francesco.vanin@epfl.ch)
Ecole Polytechnique Federale de Lausanne, Switzerland,
Earthquake Engineering and Structural Dynamics laboratory, 2019

Builds a wrapped nDMaterial (dim: 2) that uses a linear elastic model for the 
first (axial) response and a generic material model for the second (shear) response.
No coupling between the two directions.
Scope: defining the shear response of the macroelement through a uniaxial material
model, and apply it to a nDMaterial to be attached to the macroelement.

Last edit: 27 Feb 2019
*/


#include <elementAPI.h>
#include "WrappedMaterial.h"

#include <Channel.h>
#include <math.h>
#include <float.h>
#include <Information.h>
#include <Channel.h>
#include <Message.h>
#include <NDMaterial.h>
#include <UniaxialMaterial.h>
#include <MaterialResponse.h>


#ifdef _USRDLL
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif

static int numWrappedMaterial = 0;

OPS_Export void *
OPS_WrappedMaterial()
{
  // print out some KUDO's
  if (numWrappedMaterial == 0) {
    opserr << "WrappedMaterial - Loaded from external library. Written by Francesco Vanin, EPFL, 2019.\n";
    numWrappedMaterial =1;
  }

  // Pointer to the nD material that will be returned
  NDMaterial *theMaterial = 0;

  int    iData[2];
  double dData[1];
  int numData;
  numData = 2;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tags" << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid E" << endln;
    return 0;	
  }

  UniaxialMaterial* theShearModel;
  int matTag = iData[1];
  theShearModel = OPS_GetUniaxialMaterial(matTag);
  if (theShearModel == 0) {
	 opserr<<"Wrapped Material, shear model: Uniaxial Material " << matTag << " not found." << endln;
			return 0;
   }
 
  // create a new material
  theMaterial = new WrappedMaterial(iData[0], iData[1], theShearModel);

  if (theMaterial == 0) {
    opserr << "WARNING could not create NDMaterial of type WrappedMaterial\n";
    return 0;
  }
  
  // return the material
  return theMaterial;
}

// full constructor
WrappedMaterial::WrappedMaterial(int tag, double _E, UniaxialMaterial* PassedShearModel)
	:NDMaterial(tag, 0), stress(2), stressCommitted(2), u(2), uCommitted(2), Kpen(2,2), K(2,2), KCommitted(2,2), E(_E), theShearMat(0)
{ 
	// create a copy of the shear model
	theShearMat = PassedShearModel->getCopy();

   Kpen(0,0) = _E;
   Kpen(1, 1) = theShearMat->getInitialTangent();

   K = Kpen;
   KCommitted = Kpen;

}

WrappedMaterial::WrappedMaterial()
:NDMaterial(0, 0), stress(2), stressCommitted(2), u(2), uCommitted(2), Kpen(2,2), K(2,2), KCommitted(2,2), E(0.0), theShearMat(0)
{ 
}

WrappedMaterial::~WrappedMaterial() {
}


// standard methods
int 
WrappedMaterial::setTrialStrain(const Vector &strain, const Vector &rate) {
	return this->setTrialStrain(strain);
}

int
WrappedMaterial::setTrialStrain(const Vector &strain) { 
		u = strain;
		theShearMat->setTrialStrain( u(1) );

		// stress update
		stress(0) = Kpen(0,0)*u(0);

		// tangent operator update
		K(0,0) = E;

    return 0;
}

int
WrappedMaterial::setTrialStrainIncr (const Vector &strain)
{
	opserr << "WARNING: WrappedMaterial has no method setTrialStrainIncr implemented" << endln; 
	return 0;
}

int
WrappedMaterial::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  return this->setTrialStrainIncr(strain);
}

const Matrix& 
WrappedMaterial::getCommittedTangent( )  {
	return KCommitted;
}
	
const Matrix& 
WrappedMaterial::getTangent( ) 
{
	double stiff = theShearMat->getTangent();
	double stiff0 = theShearMat->getInitialTangent();

	// use a fake small stiffness if zero stiffness
	if (abs(stiff/stiff0)>0.001)
		K(1,1) = stiff;
	else
		K(1,1) = 0.001*stiff0;

	return K;
}

const Matrix&
WrappedMaterial::getInitialTangent (void)
{ 
	Kpen.Zero();
	Kpen(0,0) = E;
	Kpen(1,1) = theShearMat->getInitialTangent();

	return Kpen;
}

const Vector&
WrappedMaterial::getStress (void)
{
	stress(1) = theShearMat->getStress();
	return stress;
}

const Vector&
WrappedMaterial::getStrain (void)
{
  return u;
}

int
WrappedMaterial::commitState (void)
{
	theShearMat->commitState();

	KCommitted = this->getTangent();
	uCommitted = u;
	stressCommitted = this->getStress();

	return 0;
}

int
WrappedMaterial::revertToLastCommit (void)     
{ 
	theShearMat->revertToLastCommit();

	K = KCommitted;
	u = uCommitted;
	stress = stressCommitted;

  return 0;
}

int
WrappedMaterial::revertToStart (void)
{ 
	theShearMat->revertToStart();

	K.Zero();
	K(0,0) = E;
	K(1,1) = theShearMat->getInitialTangent();
	
	KCommitted = K;

	Kpen = K;

	u.Zero();
	stress.Zero();

  return 0;
}

NDMaterial*
WrappedMaterial::getCopy (void)
{
	WrappedMaterial *theCopy = new WrappedMaterial(this->getTag(), E, theShearMat);
  return theCopy;
}

NDMaterial*
WrappedMaterial::getCopy (const char *type) {
  return this->getCopy();
}

const char*
WrappedMaterial::getType (void) const
{
  return "BeamFiber";
}

int
WrappedMaterial::getOrder (void) const
{
  return 2;
}

//print out material data
void 
WrappedMaterial::Print( OPS_Stream &s, int flag )
{
  s << endln ;
  s << "WrappedMaterial : " ; 
  s << this->getType( ) << endln ;
  s << "UniaxialMaterial tag : " << theShearMat->getTag() << endln;
  s << "UniaxialMaterial class type : " << theShearMat->getClassType() << endln;
  s << endln ;
}

int 
WrappedMaterial::sendSelf(int commitTag, Channel &theChannel) 
{  return -1;  }

int 
WrappedMaterial::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{  return -1;  }


// methods to overwrite the standard implementation in NDMaterial for asking outputs
//----------------------------------------------------------------------------------

Response*
WrappedMaterial::setResponse (const char **argv, int argc, OPS_Stream &output) {

  Response *theResponse = 0;
  output.tag("NdMaterialOutput");
  output.attr("matType",this->getClassType());
  output.attr("matTag",this->getTag());

  if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0) {
    //const Vector &res = this->getStress();
    //int size = res.Size();
    output.tag("ResponseType","sigma");
	output.tag("ResponseType","tau");   
    theResponse =  new MaterialResponse(this, 1, this->getStress());

  } else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0) {
    //const Vector &res = this->getStrain();
    //int size = res.Size();
	output.tag("ResponseType","epsilon");
	output.tag("ResponseType","gamma");
    theResponse =  new MaterialResponse(this, 2, this->getStrain());

  } else if (strcmp(argv[0],"tangent") == 0 || strcmp(argv[0],"stiffness") == 0) {
	//const Matrix &res = this->getCommittedTangent();
    //int size = 4;
	output.tag("ResponseType","K11");
	output.tag("ResponseType","K12");
	output.tag("ResponseType","K21");
	output.tag("ResponseType","K22");
    theResponse =  new MaterialResponse(this, 3, this->getCommittedTangent());
  }

	// in-plane shear state  
    else if (strcmp(argv[0],"shearMaterial") == 0 || strcmp(argv[0],"shearMat") == 0) {
		output.tag("ShearMaterial");
		theResponse = theShearMat->setResponse(&argv[1], argc-1, output);
		output.endTag();
	}
  

  output.endTag(); // NdMaterialOutput

  return theResponse;
}

int 
WrappedMaterial::getResponse (int responseID, Information &matInfo)
{
  Vector ep_dot(2);
  double mode;

  switch (responseID) {
  case 1:
    return matInfo.setVector(this->getStress());
    
  case 2:
    return matInfo.setVector(this->getStrain());

  case 3:
	  return matInfo.setMatrix(this->getCommittedTangent());

  default:
    return -1;
  }
}
