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
// Written: fmk
// Created: Sep 2010
//
#include <math.h>

#include <InitStressMaterial.h>
#include <ID.h>
#include <Channel.h>
#include <Vector.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>

int
InitStressMaterial::findInitialStrain()
{
  // determine the initial strain
  double tol  = fabs(sigInit)*1e-12;
  double dSig = sigInit;
  double tStrain = 0.0, tStress = 0.0;
  int count = 0;

  do {
    count++;
    double K = theMaterial->getTangent();
    double dStrain = dSig/K;
    tStrain += dStrain;
    theMaterial->setTrialStrain(tStrain);
    tStress = theMaterial->getStress();
    dSig = sigInit-tStress;
  } while ((fabs(tStress-sigInit) > tol) && (count <= 100));

  epsInit = tStrain;

  if ((fabs(tStress-sigInit) < tol)) 
    theMaterial->setTrialStrain(epsInit);
  else {
    opserr << "WARNING: InitStressMaterial - could not find initStrain to within tol for material: " << theMaterial->getTag();
    opserr << " wanted sigInit: " << sigInit << " using tStress: " << theMaterial->getStress() << endln;
    return -1;
  }

  return 0;
}

InitStressMaterial::InitStressMaterial(int tag, 
                                        UniaxialMaterial &material,
                                        double sigini)
  :UniaxialMaterial(tag,MAT_TAG_InitStress), theMaterial(0),
   epsInit(0.0), sigInit(sigini)
{
  theMaterial = material.getCopy();

  if (this->findInitialStrain() == 0)
    theMaterial->commitState();
}

InitStressMaterial::InitStressMaterial()
  :UniaxialMaterial(0,MAT_TAG_InitStress), theMaterial(0),
   epsInit(0.0)
{

}

InitStressMaterial::~InitStressMaterial()
{
  if (theMaterial)
    delete theMaterial;
}

int 
InitStressMaterial::setTrialStrain(double strain, double strainRate)
{
  return theMaterial->setTrialStrain(strain+epsInit, strainRate);
}

double 
InitStressMaterial::getStress()
{
  return theMaterial->getStress();
}

double 
InitStressMaterial::getTangent()
{
  return theMaterial->getTangent();  
}

double 
InitStressMaterial::getDampTangent()
{
  return theMaterial->getDampTangent();
}

double 
InitStressMaterial::getStrain()
{
  return theMaterial->getStrain();
}

double 
InitStressMaterial::getStrainRate()
{
  return theMaterial->getStrainRate();
}

int 
InitStressMaterial::commitState()
{	
  return theMaterial->commitState();
}

int 
InitStressMaterial::revertToLastCommit()
{
  return theMaterial->revertToLastCommit();
}

int 
InitStressMaterial::revertToStart()
{
  int res = 0;
  res = theMaterial->revertToStart();
  res += theMaterial->setTrialStrain(epsInit);
  res += theMaterial->commitState();
  return res;
}

UniaxialMaterial *
InitStressMaterial::getCopy()
{
  InitStressMaterial *theCopy = 
    new InitStressMaterial(this->getTag(), *theMaterial, sigInit);
        
  return theCopy;
}

int 
InitStressMaterial::sendSelf(int cTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();

  static ID dataID(3);
  dataID(0) = this->getTag();
  dataID(1) = theMaterial->getClassTag();
  int matDbTag = theMaterial->getDbTag();
  if ( matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    theMaterial->setDbTag(matDbTag);
  }
  dataID(2) = matDbTag;
  if (theChannel.sendID(dbTag, cTag, dataID) < 0) {
    opserr << "InitStressMaterial::sendSelf() - failed to send the ID\n";
    return -1;
  }

  static Vector dataVec(1);
  dataVec(0) = epsInit;

  if (theChannel.sendVector(dbTag, cTag, dataVec) < 0) {
    opserr << "InitStressMaterial::sendSelf() - failed to send the Vector\n";
    return -2;
  }

  if (theMaterial->sendSelf(cTag, theChannel) < 0) {
    opserr << "InitStressMaterial::sendSelf() - failed to send the Material\n";
    return -3;
  }

  return 0;
}

int 
InitStressMaterial::recvSelf(int cTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static ID dataID(3);
  if (theChannel.recvID(dbTag, cTag, dataID) < 0) {
    opserr << "InitStressMaterial::recvSelf() - failed to get the ID\n";
    return -1;
  }
  this->setTag(int(dataID(0)));

  // as no way to change material, don't have to check classTag of the material 
  if (theMaterial == 0) {
    int matClassTag = int(dataID(1));
    theMaterial = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << "InitStressMaterial::recvSelf() - failed to create Material with classTag " 
	   << dataID(0) << endln;
      return -2;
    }
  }
  theMaterial->setDbTag(dataID(2));

  static Vector dataVec(1);
  if (theChannel.recvVector(dbTag, cTag, dataVec) < 0) {
    opserr << "InitStressMaterial::recvSelf() - failed to get the Vector\n";
    return -3;
  }

  epsInit = dataVec(0);
  
  if (theMaterial->recvSelf(cTag, theChannel, theBroker) < 0) {
    opserr << "InitStressMaterial::recvSelf() - failed to get the Material\n";
    return -4;
  }
  return 0;
}

void 
InitStressMaterial::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"initial_stress\": " << sigInit;
    s << "}";
    return;
  }
  else {
    s << "InitStressMaterial tag: " << this->getTag() << endln;
    s << "\tMaterial: " << theMaterial->getTag() << endln;
    s << "\tinitital strain: " << epsInit << endln;
  }
}

int 
InitStressMaterial::setParameter(const char **argv, int argc, Parameter &param)
{
  if (strcmp(argv[0],"sig0") == 0 || strcmp(argv[0],"f0") == 0 || strcmp(argv[0],"F0") == 0) {
    param.setValue(sigInit);
    return param.addObject(1, this);
  }
  return theMaterial->setParameter(argv, argc, param);
}

int
InitStressMaterial::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
  case -1:
    return -1;
  case 1:
    this->sigInit = info.theDouble;
    this->findInitialStrain();
    break;
  }

  return 0;
}

double
InitStressMaterial::getStressSensitivity(int gradIndex, bool conditional)
{
  return theMaterial->getStressSensitivity(gradIndex, conditional);
}

double
InitStressMaterial::getInitialTangentSensitivity(int gradIndex)
{
  return theMaterial->getInitialTangentSensitivity(gradIndex);
}

int
InitStressMaterial::commitSensitivity(double strainGradient, 
                                      int gradIndex, int numGrads)
{
  return theMaterial->commitSensitivity(strainGradient, gradIndex, numGrads);
}
