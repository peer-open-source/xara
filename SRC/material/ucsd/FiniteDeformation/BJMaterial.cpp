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
**                                                                    **
** Additions and changes by:                                          **
**   Boris Jeremic (@ucdavis.edu)                                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.26 $                                                              
// $Date: 2010-09-13 21:29:28 $                                                                  
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/BJMaterial.cpp,v $                                                                
                                                                        
// File: ~/material/BJMaterial.C
//
// Written: MHS 
// Created: Feb 2000
// Revision: A
//
// Description: This file contains the class implementation for BJMaterial.
//
// What: "@(#) BJMaterial.C, revA"

#include "BJMaterial.h"
#include <Information.h>
#include <OPS_Globals.h>
#include <Matrix.h>
#include <Vector.h>
#include <stresst.h>
#include <straint.h>
#include <MaterialResponse.h>

#include <PlaneStressMaterial.h>
#include <BeamFiberMaterial.h>
#include <PlateFiberMaterial.h>

Matrix BJMaterial::errMatrix(1,1);
Vector BJMaterial::errVector(1);
Tensor BJMaterial::errTensor(2, def_dim_2, 0.0 );
stresstensor BJMaterial::errstresstensor;
straintensor BJMaterial::errstraintensor;

BJMaterial::BJMaterial(int tag, int classTag)
:NDMaterial(tag,classTag)
{

}

BJMaterial::BJMaterial()
:NDMaterial(0, 0)
{

}

BJMaterial::~BJMaterial()
{

}


NDMaterial*
BJMaterial::getCopy(void) {
  return (NDMaterial*)(this->getCopyBJ());
}

#if 0
NDMaterial*
BJMaterial::getCopy(const char *type)
{

  if (strcmp(type,"PlaneStress") == 0 ||
      strcmp(type,"PlaneStress2D") == 0) {
    BJMaterial *copy = this->getCopyBJ("ThreeDimensional");
    PlaneStressMaterial *clone = new PlaneStressMaterial(this->NDMaterial::getTag(),*copy);
    return clone;
  }
  else if (strcmp(type,"BeamFiber") == 0 ||
	   strcmp(type,"TimoshenkoFiber") == 0) {
    BJMaterial *copy = this->getCopyBJ("ThreeDimensional");
    BeamFiberMaterial *clone = new BeamFiberMaterial(this->NDMaterial::getTag(),*copy);
    return clone;
  }
  else if (strcmp(type,"PlateFiber") == 0) {
    BJMaterial *copy = this->getCopyBJ("ThreeDimensional");
    PlateFiberMaterial *clone = new PlateFiberMaterial(this->NDMaterial::getTag(),*copy);
    return clone;
  }
  else
    return nullptr;
}
#endif


double
BJMaterial::getRho(void)
{
  return 0.0;
}

// methods to set and retrieve state.
int 
BJMaterial::setTrialStrain(const Vector &v)
{
   opserr << "BJMaterial::setTrialStrain -- subclass responsibility\n";
   return -1;    
}

int 
BJMaterial::setTrialStrain(const Vector &v, const Vector &r)
{
   opserr << "BJMaterial::setTrialStrain -- subclass responsibility\n";
   return -1;    
}

int 
BJMaterial::setTrialStrainIncr(const Vector &v)
{
   opserr << "BJMaterial::setTrialStrainIncr -- subclass responsibility\n";
   return -1;    
}

int 
BJMaterial::setTrialStrainIncr(const Vector &v, const Vector &r)
{
   opserr << "BJMaterial::setTrialStrainIncr -- subclass responsibility\n";
   return -1;    
}

const Matrix &
BJMaterial::getTangent(void)
{
   opserr << "BJMaterial::getTangent -- subclass responsibility\n";
   return errMatrix;    
}

const Vector &
BJMaterial::getStress(void)
{
   opserr << "BJMaterial::getStress -- subclass responsibility\n";
   return errVector;    
}

const Vector &
BJMaterial::getStrain(void)
{
   opserr << "BJMaterial::getStrain -- subclass responsibility\n";
   return errVector;    
}

int 
BJMaterial::setTrialStrain(const Tensor &v)
{
   opserr << "BJMaterial::setTrialStrainIncr -- subclass responsibility\n";
   return -1;    
}

int 
BJMaterial::setTrialStrain(const Tensor &v, const Tensor &r)    
{
   opserr << "BJMaterial::setTrialStrainIncr -- subclass responsibility\n";
   return -1;    
}

int 
BJMaterial::setTrialStrainIncr(const Tensor &v)
{
   opserr << "BJMaterial::setTrialStrainIncr -- subclass responsibility\n";
   return -1;    
}

int 
BJMaterial::setTrialStrainIncr(const Tensor &v, const Tensor &r)
{
   opserr << "BJMaterial::setTrialStrainIncr -- subclass responsibility\n";
   return -1;    
}


// added Sept 22 2003 for Large Deformation, F is the Deformation Grandient

int
BJMaterial::setTrialF(const straintensor &f)
{
   opserr << "BJMaterial::setTrialF -- subclass responsibility\n";
   return -1;
}

int
BJMaterial::setTrialFIncr(const straintensor &df)
{
   opserr << "BJMaterial::setTrialF -- subclass responsibility\n";
   return -1;
}

int
BJMaterial::setTrialC(const straintensor &c)
{
   opserr << "BJMaterial::setTrialC -- subclass responsibility\n";
   return -1;
}

int
BJMaterial::setTrialCIncr(const straintensor &c)
{
   opserr << "BJMaterial::setTrialC -- subclass responsibility\n";
   return -1;
}

const stresstensor& BJMaterial::getPK1StressTensor(void)
{
   opserr << "BJMaterial::getPK1StressTensor -- subclass responsibility\n";
   return errstresstensor;    
}

const stresstensor& BJMaterial::getCauchyStressTensor(void)
{
   opserr << "BJMaterial::getCauchyStressTensor -- subclass responsibility\n";
   return errstresstensor;    
}

const straintensor& BJMaterial::getF(void)
{
   opserr << "BJMaterial::getF -- subclass responsibility\n";
   return errstraintensor;    
}

const straintensor& BJMaterial::getC(void)
{
   opserr << "BJMaterial::getF -- subclass responsibility\n";
   return errstraintensor;    
}

const straintensor& BJMaterial::getFp(void)
{
   opserr << "BJMaterial::getFp -- subclass responsibility\n";
   return errstraintensor;    
}
// Only For Large Deformation, END////////////////////////////

const Tensor &
BJMaterial::getTangentTensor(void)
{
   opserr << "BJMaterial::getTangentTensor -- subclass responsibility\n";
   return errTensor;    
}

const stresstensor& BJMaterial::getStressTensor(void)
{
   opserr << "BJMaterial::getStressTensor -- subclass responsibility\n";
   return errstresstensor;    
}

const straintensor& BJMaterial::getStrainTensor(void)
{
   opserr << "BJMaterial::getStrainTensor -- subclass responsibility\n";
   return errstraintensor;    
}

const straintensor& BJMaterial::getPlasticStrainTensor(void)
{
   opserr << "BJMaterial::getPlasticStrainTensor -- subclass responsibility\n";
   return errstraintensor;    
}


//const Tensor &
//BJMaterial::getStrainTensor(void)
//{
//   opserr << "BJMaterial::getStrainTensor -- subclass responsibility\n";
//   return errTensor;    
//}

#if 0
Response*
BJMaterial::setResponse (const char **argv, int argc, 
			 OPS_Stream &output)
{
  Response *theResponse =0;
  const char *matType = this->getType();

  output.tag("NdMaterialOutput");
  output.attr("matType",this->getClassType());
  output.attr("matTag",this->getTag());

  if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0) {
    const Vector &res = this->getStress();
    int size = res.Size();
    
    if ( (strcmp(matType,"PlaneStress") == 0 && size == 3) ||
	 (strcmp(matType,"PlaneStrain") == 0 && size == 3)) {
	output.tag("ResponseType","sigma11");
	output.tag("ResponseType","sigma22");
	output.tag("ResponseType","sigma12");
    } else if (strcmp(matType,"ThreeDimensional") == 0 && size == 6) {
	output.tag("ResponseType","sigma11");
	output.tag("ResponseType","sigma22");
	output.tag("ResponseType","sigma33");
	output.tag("ResponseType","sigma12");
	output.tag("ResponseType","sigma13");
	output.tag("ResponseType","sigma23");
    } else {
      for (int i=0; i<size; i++) 
	output.tag("ResponseType","UnknownStress");
    }
    theResponse =  new MaterialResponse(this, 1, this->getStress());

  } else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0) {
    const Vector &res = this->getStrain();
    int size = res.Size();
    if ( (strcmp(matType,"PlaneStress") == 0 && size == 3) ||
	 (strcmp(matType,"PlaneStrain") == 0 && size == 3)) {
	output.tag("ResponseType","eta11");
	output.tag("ResponseType","eta22");
	output.tag("ResponseType","eta12");
    } else if (strcmp(matType,"ThreeDimensional") == 0 && size == 6) {
	output.tag("ResponseType","eps11");
	output.tag("ResponseType","eps22");
	output.tag("ResponseType","eps33");
	output.tag("ResponseType","eps12");
	output.tag("ResponseType","eps13");
	output.tag("ResponseType","eps23");
    } else {
      for (int i=0; i<size; i++) 
	output.tag("ResponseType","UnknownStrain");
    }      
    theResponse =  new MaterialResponse(this, 2, this->getStress());
  }

  output.endTag(); // NdMaterialOutput

  return theResponse;
}

int 
BJMaterial::getResponse (int responseID, Information &matInfo)
{
  switch (responseID) {
  case 1:
    return matInfo.setVector(this->getStress());
    
  case 2:
    return matInfo.setVector(this->getStrain());
    
  default:
    return -1;
  }
}
#endif



// AddingSensitivity:BEGIN ////////////////////////////////////////
const Vector &
BJMaterial::getStressSensitivity(int gradIndex, bool conditional)
{
	static Vector dummy(1);
	return dummy;
}

const Vector &
BJMaterial::getStrainSensitivity(int gradIndex)
{
	static Vector dummy(1);
	return dummy;
}

double
BJMaterial::getRhoSensitivity(int gradIndex)
{
	return 0.0;
}

const Matrix &
BJMaterial::getDampTangentSensitivity(int gradIndex)
{
	static Matrix dummy(1,1);
	return dummy;
}
const Matrix &
BJMaterial::getTangentSensitivity(int gradIndex)
{
	static Matrix dummy(1,1);
	return dummy;
}
const Matrix &
BJMaterial::getInitialTangentSensitivity(int gradIndex)
{
	static Matrix dummy(1,1);
	return dummy;
}
int
BJMaterial::commitSensitivity(Vector & strainSensitivity, int gradIndex, int numGrads)
{
	return 0;
}
// AddingSensitivity:END //////////////////////////////////////////


