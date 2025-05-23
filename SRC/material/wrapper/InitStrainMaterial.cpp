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
                                                                        
// $Revision: 1.1 $
// $Date: 2010-09-11 00:45:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/InitStrainMaterial.cpp,v $

// Written: fmk
// Created: Sep 2010
//
#include <InitStrainMaterial.h>
#include <ID.h>
#include <Channel.h>
#include <Vector.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>


InitStrainMaterial::InitStrainMaterial(int tag, 
                                       UniaxialMaterial &material,
                                       double epsini)
  :UniaxialMaterial(tag,MAT_TAG_InitStrain), theMaterial(0),
   epsInit(epsini), localStrain(0.0)
{
  theMaterial = material.getCopy();
  theMaterial->setTrialStrain(epsInit);
  theMaterial->commitState();
}

InitStrainMaterial::InitStrainMaterial()
  :UniaxialMaterial(0,MAT_TAG_InitStrain), theMaterial(0),
   epsInit(0.0), localStrain(0.0)
{

}

InitStrainMaterial::~InitStrainMaterial()
{
  if (theMaterial)
    delete theMaterial;
}

int 
InitStrainMaterial::setTrialStrain(double strain, double strainRate)
{
  localStrain = strain;

  if (theMaterial)
    return theMaterial->setTrialStrain(strain+epsInit, strainRate);
  else
    return -1;
}

double 
InitStrainMaterial::getStress()
{
  if (theMaterial)
    return theMaterial->getStress();
  else
    return 0.0;
}

double 
InitStrainMaterial::getTangent()
{
  if (theMaterial)
    return theMaterial->getTangent();
  else
    return 0.0;
}

double 
InitStrainMaterial::getDampTangent()
{
  if (theMaterial)
    return theMaterial->getDampTangent();
  else
    return 0.0;
}

double 
InitStrainMaterial::getStrain()
{
  if (theMaterial)
    return theMaterial->getStrain();
  else
    return 0.0;
}

double 
InitStrainMaterial::getStrainRate()
{
  if (theMaterial)
    return theMaterial->getStrainRate();
  else
    return 0.0;
}

int
InitStrainMaterial::commitState()
{
  if (theMaterial)
    return theMaterial->commitState();
  else
    return -1;
}

int 
InitStrainMaterial::revertToLastCommit()
{
  if (theMaterial)
    return theMaterial->revertToLastCommit();
  else
    return -1;
}

int 
InitStrainMaterial::revertToStart()
{
  int res = 0;
  if (theMaterial) {
    res = theMaterial->revertToStart();
    res += theMaterial->setTrialStrain(epsInit);
    res += theMaterial->commitState();
    return res;
  } else
    return -1;
}

UniaxialMaterial *
InitStrainMaterial::getCopy(void)
{
  InitStrainMaterial *theCopy = 0;
  if (theMaterial)
    theCopy = new InitStrainMaterial(this->getTag(), *theMaterial, epsInit);
        
  return theCopy;
}

int 
InitStrainMaterial::sendSelf(int cTag, Channel &theChannel)
{
  if (theMaterial == 0) {
    opserr << "InitStrainMaterial::sendSelf() - theMaterial is null, nothing to send\n";
    return -1;
  }
  
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
    opserr << "InitStrainMaterial::sendSelf() - failed to send the ID\n";
    return -1;
  }

  static Vector dataVec(2);
  dataVec(0) = epsInit;
  dataVec(1) = localStrain;
  
  if (theChannel.sendVector(dbTag, cTag, dataVec) < 0) {
    opserr << "InitStrainMaterial::sendSelf() - failed to send the Vector\n";
    return -2;
  }

  if (theMaterial->sendSelf(cTag, theChannel) < 0) {
    opserr << "InitStrainMaterial::sendSelf() - failed to send the Material\n";
    return -3;
  }

  return 0;
}

int 
InitStrainMaterial::recvSelf(int cTag, Channel &theChannel, 
                         FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static ID dataID(3);
  if (theChannel.recvID(dbTag, cTag, dataID) < 0) {
    opserr << "InitStrainMaterial::recvSelf() - failed to get the ID\n";
    return -1;
  }
  this->setTag(int(dataID(0)));

  // as no way to change material, don't have to check classTag of the material 
  if (theMaterial == 0) {
    int matClassTag = int(dataID(1));
    theMaterial = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << "InitStrainMaterial::recvSelf() - failed to create Material with classTag " 
           << dataID(0) << endln;
      return -2;
    }
  }
  theMaterial->setDbTag(dataID(2));

  static Vector dataVec(2);
  if (theChannel.recvVector(dbTag, cTag, dataVec) < 0) {
    opserr << "InitStrainMaterial::recvSelf() - failed to get the Vector\n";
    return -3;
  }

  epsInit = dataVec(0);
  localStrain = dataVec(1);
  
  if (theMaterial->recvSelf(cTag, theChannel, theBroker) < 0) {
    opserr << "InitStrainMaterial::recvSelf() - failed to get the Material\n";
    return -4;
  }
  return 0;
}

void 
InitStrainMaterial::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"InitStrainMaterial\", ";
    if (theMaterial)
      s << "\"Material\": " << theMaterial->getTag() << ", ";
    else
      s << "\"Material\": " << "null" << ", ";
    s << "\"initial_strain\": " << epsInit;
    s <<  "}";
    return;
  }
  else {
    s << "InitStrainMaterial tag: " << this->getTag() << endln;
    if (theMaterial)
      s << "\tMaterial: " << theMaterial->getTag() << endln;
    else
      s << "\tMaterial is NULL" << endln;
    s << "\tinitital strain: " << epsInit << endln;
  }
}

int 
InitStrainMaterial::setParameter(const char **argv, int argc, Parameter &param)
{
  if (strcmp(argv[0], "initial_strain")==0 || strcmp(argv[0],"epsInit") == 0) {
    param.setValue(epsInit);
    return param.addObject(1, this);
  }

  // Otherwise, pass it on to the wrapped material
  if (theMaterial)
    return theMaterial->setParameter(argv, argc, param);
  else
    return -1;
}

int
InitStrainMaterial::updateParameter(int parameterID, Information &info)
{
  if (parameterID == 1) {
    this->epsInit = info.theDouble;
    if (theMaterial) {
      theMaterial->setTrialStrain(localStrain+epsInit);
      theMaterial->commitState();
    } else
      return -1;
  }

  return 0;
}

double
InitStrainMaterial::getStressSensitivity(int gradIndex, bool conditional)
{
  if (theMaterial)
    return theMaterial->getStressSensitivity(gradIndex, conditional);
  else
    return 0.0;
}

double
InitStrainMaterial::getInitialTangentSensitivity(int gradIndex)
{
  if (theMaterial)
    return theMaterial->getInitialTangentSensitivity(gradIndex);
  else
    return 0.0;
}

int
InitStrainMaterial::commitSensitivity(double strainGradient, 
                                      int gradIndex, int numGrads)
{
  if (theMaterial)
    return theMaterial->commitSensitivity(strainGradient, gradIndex, numGrads);
  else
    return -1;
}
