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
// Written: MHS 
// Created: Feb 2000
// Revision: A
//
// Description: This file contains the class implementation for SectionForceDeformation.
//
//
#include <SectionForceDeformation.h>
#include <Information.h>
#include <Matrix.h>
#include <Vector.h>

#include <SensitiveResponse.h>
typedef SensitiveResponse<SectionForceDeformation> SectionResponse;

#include <string.h>

#include <TaggedObject.h>
#include <api/runtimeAPI.h>


SectionForceDeformation::SectionForceDeformation(int tag, int classTag)
  :TaggedObject(tag), MovableObject(classTag), fDefault(0), sDefault(0)
{

}

SectionForceDeformation::SectionForceDeformation()
    : TaggedObject(0),MovableObject(0), fDefault(0), sDefault(0)
{

}

SectionForceDeformation::~SectionForceDeformation()
{
  if (fDefault != 0)
    delete fDefault;
  if (sDefault != 0)
    delete sDefault;
}

const Matrix&
SectionForceDeformation::getSectionFlexibility ()
{
  int order = this->getOrder();
  
  if (fDefault == nullptr) {                
    fDefault = new Matrix(order,order);
  }

  const Matrix &k = this->getSectionTangent();
  
  switch(order) {
  case 1:
    if (k(0,0) != 0.0)
      (*fDefault)(0,0) = 1.0/k(0,0);
    break;
  default:
    k.Invert(*fDefault);
    break;
  }

  return *fDefault;
}

const Matrix&
SectionForceDeformation::getInitialFlexibility ()
{
  int order = this->getOrder();
  
  if (fDefault == 0) {                
    fDefault = new Matrix(order,order);
    if (fDefault == 0) {
      opserr << "SectionForceDeformation::getInitialFlexibility -- failed to allocate flexibility matrix\n";
      exit(-1);
    }
  }
  
  const Matrix &k = this->getInitialTangent();
  
  switch(order) {
  case 1:
    if (k(0,0) != 0.0)
      (*fDefault)(0,0) = 1.0/k(0,0);
    break;
  default:
    k.Invert(*fDefault);
    break;
  }
  
  return *fDefault;
}

double 
SectionForceDeformation::getRho(void) 
{
  return 0.0 ;
}

Response*
SectionForceDeformation::setResponse(const char **argv, int argc,
                                     OPS_Stream &output)
{
  const ID &type = this->getType();
  int typeSize = this->getOrder();
  
  Response *theResponse =0;

  output.tag("SectionOutput");
  output.attr("secType", this->getClassType());
  output.attr("secTag", this->getTag());

  // deformations
  if (strcmp(argv[0],"deformations") == 0 || strcmp(argv[0],"deformation") == 0) {
    for (int i=0; i<typeSize; i++) {
      int code = type(i);
      switch (code){
      case SECTION_RESPONSE_MZ:
        output.tag("ResponseType","kappaZ");
        break;
      case SECTION_RESPONSE_P:
        output.tag("ResponseType","eps");
        break;
      case SECTION_RESPONSE_VY:
        output.tag("ResponseType","gammaY");
        break;
      case SECTION_RESPONSE_MY:
        output.tag("ResponseType","kappaY");
        break;
      case SECTION_RESPONSE_VZ:
        output.tag("ResponseType","gammaZ");
        break;
      case SECTION_RESPONSE_T:
        output.tag("ResponseType","theta");
        break;
      case SECTION_RESPONSE_FXX:
          output.tag("ResponseType", "epsXX");
          break;
      case SECTION_RESPONSE_FYY:
          output.tag("ResponseType", "epsYY");
          break;
      case SECTION_RESPONSE_FXY:
          output.tag("ResponseType", "epsXY");
          break;
      case SECTION_RESPONSE_MXX:
          output.tag("ResponseType", "kappaXX");
          break;
      case SECTION_RESPONSE_MYY:
          output.tag("ResponseType", "kappaYY");
          break;
      case SECTION_RESPONSE_MXY:
          output.tag("ResponseType", "kappaXY");
          break;
      case SECTION_RESPONSE_VXZ:
          output.tag("ResponseType", "gammaXZ");
          break;
      case SECTION_RESPONSE_VYZ:
          output.tag("ResponseType", "gammaYZ");
          break;
      default:
        output.tag("ResponseType","Unknown");
      }
    }
    theResponse =  new SectionResponse(*this, 1, this->getSectionDeformation());
  
  // forces
  } else if (strcmp(argv[0], "resultant")==0 || 
             strcmp(argv[0], "forces") == 0 || 
             strcmp(argv[0], "force") == 0) {
    for (int i=0; i<typeSize; i++) {
      int code = type(i);
      switch (code){
      case SECTION_RESPONSE_MZ:
        output.tag("ResponseType","Mz");
        break;
      case SECTION_RESPONSE_P:
        output.tag("ResponseType","P");
        break;
      case SECTION_RESPONSE_VY:
        output.tag("ResponseType","Vy");
        break;
      case SECTION_RESPONSE_MY:
        output.tag("ResponseType","My");
        break;
      case SECTION_RESPONSE_VZ:
        output.tag("ResponseType","Vz");
        break;
      case SECTION_RESPONSE_T:
        output.tag("ResponseType","T");
        break;
      case SECTION_RESPONSE_FXX:
          output.tag("ResponseType", "Fxx");
          break;
      case SECTION_RESPONSE_FYY:
          output.tag("ResponseType", "Fyy");
          break;
      case SECTION_RESPONSE_FXY:
          output.tag("ResponseType", "Fxy");
          break;
      case SECTION_RESPONSE_MXX:
          output.tag("ResponseType", "Mxx");
          break;
      case SECTION_RESPONSE_MYY:
          output.tag("ResponseType", "Myy");
          break;
      case SECTION_RESPONSE_MXY:
          output.tag("ResponseType", "Mxy");
          break;
      case SECTION_RESPONSE_VXZ:
          output.tag("ResponseType", "Vxz");
          break;
      case SECTION_RESPONSE_VYZ:
          output.tag("ResponseType", "Vyz");
          break;
      default:
        output.tag("ResponseType","Unknown");
      }
    }
    theResponse =  new SectionResponse(*this, 2, this->getStressResultant());
  
  // force and deformation
  } else if (strcmp(argv[0],"forceAndDeformation") == 0) { 
    for (int i=0; i<typeSize; i++) {
      int code = type(i);
      switch (code){
      case SECTION_RESPONSE_MZ:
        output.tag("ResponseType","kappaZ");
        break;
      case SECTION_RESPONSE_P:
        output.tag("ResponseType","eps");
        break;
      case SECTION_RESPONSE_VY:
        output.tag("ResponseType","gammaY");
        break;
      case SECTION_RESPONSE_MY:
        output.tag("ResponseType","kappaY");
        break;
      case SECTION_RESPONSE_VZ:
        output.tag("ResponseType","gammaZ");
        break;
      case SECTION_RESPONSE_T:
        output.tag("ResponseType","theta");
        break;
      case SECTION_RESPONSE_FXX:
          output.tag("ResponseType", "epsXX");
          break;
      case SECTION_RESPONSE_FYY:
          output.tag("ResponseType", "epsYY");
          break;
      case SECTION_RESPONSE_FXY:
          output.tag("ResponseType", "epsXY");
          break;
      case SECTION_RESPONSE_MXX:
          output.tag("ResponseType", "kappaXX");
          break;
      case SECTION_RESPONSE_MYY:
          output.tag("ResponseType", "kappaYY");
          break;
      case SECTION_RESPONSE_MXY:
          output.tag("ResponseType", "kappaXY");
          break;
      case SECTION_RESPONSE_VXZ:
          output.tag("ResponseType", "gammaXZ");
          break;
      case SECTION_RESPONSE_VYZ:
          output.tag("ResponseType", "gammaYZ");
          break;
      default:
        output.tag("ResponseType","Unknown");
      }
    }
    for (int j=0; j<typeSize; j++) {
      int code = type(j);
      switch (code){
      case SECTION_RESPONSE_MZ:
        output.tag("ResponseType","Mz");
        break;
      case SECTION_RESPONSE_P:
        output.tag("ResponseType","P");
        break;
      case SECTION_RESPONSE_VY:
        output.tag("ResponseType","Vy");
        break;
      case SECTION_RESPONSE_MY:
        output.tag("ResponseType","My");
        break;
      case SECTION_RESPONSE_VZ:
        output.tag("ResponseType","Vz");
        break;
      case SECTION_RESPONSE_T:
        output.tag("ResponseType","T");
        break;
      case SECTION_RESPONSE_FXX:
          output.tag("ResponseType", "Fxx");
          break;
      case SECTION_RESPONSE_FYY:
          output.tag("ResponseType", "Fyy");
          break;
      case SECTION_RESPONSE_FXY:
          output.tag("ResponseType", "Fxy");
          break;
      case SECTION_RESPONSE_MXX:
          output.tag("ResponseType", "Mxx");
          break;
      case SECTION_RESPONSE_MYY:
          output.tag("ResponseType", "Myy");
          break;
      case SECTION_RESPONSE_MXY:
          output.tag("ResponseType", "Mxy");
          break;
      case SECTION_RESPONSE_VXZ:
          output.tag("ResponseType", "Vxz");
          break;
      case SECTION_RESPONSE_VYZ:
          output.tag("ResponseType", "Vyz");
          break;
      default:
        output.tag("ResponseType","Unknown");
      }
    }
    
    theResponse =  new SectionResponse(*this, 4, Vector(2*this->getOrder()));
  }
  
  else if (strcmp(argv[0],"stiffness") == 0) {
    theResponse =  new SectionResponse(*this, 12, this->getSectionTangent());
  }

  else if (strcmp(argv[0],"flexibility") == 0) {
    theResponse =  new SectionResponse(*this, 13, this->getSectionFlexibility());
  }
  

  output.endTag(); // SectionOutput
  return theResponse;
}

int 
SectionForceDeformation::getResponse(int responseID, Information &secInfo)
{
  switch (responseID) {
  case 1:
    return secInfo.setVector(this->getSectionDeformation());
    
  case 2:
    return secInfo.setVector(this->getStressResultant());
    
  case 4: {
    Vector &theVec = *(secInfo.theVector);
    const Vector &e = this->getSectionDeformation();
    const Vector &s = this->getStressResultant();
    int order = this->getOrder();
    for (int i = 0; i < order; i++) {
      theVec(i) = e(i);
      theVec(i+order) = s(i);
    }
    
    return secInfo.setVector(theVec);
  }

  case 12:
    return secInfo.setMatrix(this->getSectionTangent());

  case 13:
    return secInfo.setMatrix(this->getSectionFlexibility());

  default:
    return -1;
  }
}

int 
SectionForceDeformation::getResponseSensitivity(int responseID, int gradIndex,
                                                Information &secInfo)
{
  Vector &theVec = *(secInfo.theVector);

  switch (responseID) {
  case 1:
    theVec = this->getSectionDeformationSensitivity(gradIndex);
    return secInfo.setVector(theVec);
    
  case 2: {
    const Matrix &ks = this->getSectionTangent();
    const Vector &dedh = this->getSectionDeformationSensitivity(gradIndex);
    const Vector &dsdh = this->getStressResultantSensitivity(gradIndex, true);
    theVec.addMatrixVector(0.0, ks, dedh, 1.0);
    theVec.addVector(1.0, dsdh, 1.0);
    return secInfo.setVector(theVec);
  }

  default:
    return -1;
  }
}

// AddingSensitivity:BEGIN ////////////////////////////////////////
const Vector &
SectionForceDeformation::getStressResultantSensitivity(int gradIndex, bool conditional)
{
  if (sDefault == 0)
    sDefault = new Vector (this->getOrder());
  return *sDefault;
}

const Vector &
SectionForceDeformation::getSectionDeformationSensitivity(int gradIndex)
{
  if (sDefault == 0)
    sDefault = new Vector (this->getOrder());
  return *sDefault;
}

const Matrix &
SectionForceDeformation::getSectionTangentSensitivity(int gradIndex)
{
  int order = this->getOrder();
  
  if (fDefault == 0) {                
    fDefault = new Matrix(order,order);
    if (fDefault == 0) {
      opserr << "SectionForceDeformation::getSectionTangentSensitivity -- failed to allocate matrix\n";
      exit(-1);
    }
  }

  fDefault->Zero();

  return *fDefault;
}

const Matrix &
SectionForceDeformation::getInitialTangentSensitivity(int gradIndex)
{
  int order = this->getOrder();
  
  if (fDefault == 0) {                
    fDefault = new Matrix(order,order);
    if (fDefault == 0) {
      opserr << "SectionForceDeformation::getInitialTangentSensitivity -- failed to allocate matrix\n";
      exit(-1);
    }
  }

  fDefault->Zero();

  return *fDefault;
}

const Matrix&
SectionForceDeformation::getSectionFlexibilitySensitivity(int gradIndex)
{
  int order = this->getOrder();
  
  if (fDefault == 0) {                
    fDefault = new Matrix(order,order);
    if (fDefault == 0) {
      opserr << "SectionForceDeformation::getSectionFlexibilitySensitivity -- failed to allocate matrix\n";
      exit(-1);
    }
  }

  const Matrix &dksdh = this->getSectionTangentSensitivity(gradIndex);
  
  const Matrix &fs = this->getSectionFlexibility();

  *fDefault = (fs * dksdh * fs) * -1;

  return *fDefault;
}

const Matrix&
SectionForceDeformation::getInitialFlexibilitySensitivity(int gradIndex)
{
  int order = this->getOrder();
  
  if (fDefault == 0) {                
    fDefault = new Matrix(order,order);
    if (fDefault == 0) {
      opserr << "SectionForceDeformation::getInitialFlexibilitySensitivity -- failed to allocate matrix\n";
      exit(-1);
    }
  }
  
  const Matrix &dksdh = this->getInitialTangentSensitivity(gradIndex);
  
  const Matrix &fs = this->getInitialFlexibility();

  *fDefault = (fs * dksdh * fs) * -1;
  
  return *fDefault;
}

double
SectionForceDeformation::getRhoSensitivity(int gradIndex)
{
  return 0.0;
}

int
SectionForceDeformation::commitSensitivity(const Vector& defSens,
                                           int gradIndex, int numGrads)
{
  return -1;
}
// AddingSensitivity:END ///////////////////////////////////////////

//--- Adding Thermal Functions:[BEGIN]   by UoE OpenSees Group ----//
int
SectionForceDeformation::setTrialSectionDeformation (const Vector&) //JZ
{
  opserr << "SectionForceDeformation::setTrialSectionDeformation(strain) - should not be called\n";
  return -1;
}

int
SectionForceDeformation::setTrialSectionDeformation(const Vector& nouse, const Vector &data) //JZ
{
  opserr << "SectionForceDeformation::setTrialSectionDeformationTemperature (strain, tData) - should not be called\n";
  return -1;
}

static Vector errRes(3);

const Vector &
SectionForceDeformation::getTemperatureStress(const Vector &tData) //PK
{
  opserr << "SectionForceDeformation::getTemperatureStress(double *dataMixed) - should not be called\n";
  errRes.resize(this->getStressResultant().Size());
  return errRes;
  //  return this->getStressResultant();
}
//--- Adding Thermal Functions:[END]   by UoE OpenSees Group ----//

const Vector& SectionForceDeformation::getThermalElong(void)
{
  errRes.resize(this->getStressResultant().Size());
  return errRes;
}

