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
// Created: 07/98
// Revision: A
//
// Description: This file contains the class implementation for 
// ElasticMaterial. 
//
#include <ElasticMaterial.h>
#include <Vector.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <string.h>


ElasticMaterial::ElasticMaterial(int tag, double E)
: UniaxialMaterial(tag,MAT_TAG_ElasticMaterial),
  trialStrain(0.0),  
  trialStrainRate(0.0),
  committedStrain(0.0), 
  committedStrainRate(0.0),
  Epos(E), Eneg(E), eta(0.0), density(0.0),
  parameterID(0)
{

}


ElasticMaterial::ElasticMaterial(int tag, double Epos, double eta, double Eneg, double density)
: UniaxialMaterial(tag,MAT_TAG_ElasticMaterial),
  trialStrain(0.0),  trialStrainRate(0.0),
  committedStrain(0.0),  committedStrainRate(0.0),
  Epos(Epos), Eneg(Eneg), 
  eta(eta), 
  density(density),
  parameterID(0)
{

}


ElasticMaterial::~ElasticMaterial()
{
  // does nothing
}


int 
ElasticMaterial::setTrialStrain(double strain, double strainRate)
{
    trialStrain     = strain;
    trialStrainRate = strainRate;
    return 0;
}


int 
ElasticMaterial::setTrial(double strain, double &stress, double &tangent, double strainRate)
{
    trialStrain     = strain;
    trialStrainRate = strainRate;

    if (trialStrain >= 0.0) {
        stress = Epos*trialStrain + eta*trialStrainRate;
        tangent = Epos;
    } else {
        stress = Eneg*trialStrain + eta*trialStrainRate;
        tangent = Eneg;
    }

    return 0;
}


double 
ElasticMaterial::getStress(void)
{
    if (trialStrain >= 0.0)
        return Epos*trialStrain + eta*trialStrainRate;
    else
        return Eneg*trialStrain + eta*trialStrainRate;
}


double 
ElasticMaterial::getTangent()
{
    if (trialStrain > 0.0)
        return Epos;
    else if (trialStrain < 0.0)
        return Eneg;
    else
        return (Epos > Eneg) ? Epos : Eneg;
}


double 
ElasticMaterial::getInitialTangent()
{
    return (Epos > Eneg) ? Epos : Eneg;
}


double 
ElasticMaterial::getRho()
{
  return density;
}

int 
ElasticMaterial::commitState()
{
  committedStrain = trialStrain;
  committedStrainRate = trialStrainRate;
  return 0;
}


int 
ElasticMaterial::revertToLastCommit()
{
  trialStrain = committedStrain;
  trialStrainRate = committedStrainRate;
  return 0;
}


int 
ElasticMaterial::revertToStart()
{
  trialStrain      = 0.0;
  trialStrainRate  = 0.0;
  return 0;
}


UniaxialMaterial *
ElasticMaterial::getCopy()
{
  ElasticMaterial *theCopy = new ElasticMaterial(this->getTag(),Epos,eta,Eneg,density);
  theCopy->trialStrain     = trialStrain;
  theCopy->trialStrainRate = trialStrainRate;
  theCopy->committedStrain     = committedStrain;
  theCopy->committedStrainRate = committedStrainRate;
  return theCopy;
}




void 
ElasticMaterial::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
    s << "ElasticMaterial tag: " << this->getTag() << "\n";
    s << "  Epos: " << Epos << " Eneg: " << Eneg << " eta: " << eta << "\n";
  } 
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"ElasticMaterial\", ";
    s << "\"Epos\": " << Epos << ", ";
    s << "\"Eneg\": " << Eneg << ", ";
    s << "\"density\": " << density << ", ";
    s << "\"eta\": " << eta << "}";
  }
}


int
ElasticMaterial::setParameter(const char **argv, int argc, Parameter &param)
{

  if (strcmp(argv[0],"E") == 0) {
    param.setValue(Epos);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"Epos") == 0) {
    param.setValue(Epos);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"Eneg") == 0) {
    param.setValue(Eneg);
    return param.addObject(3, this);
  }
  if (strcmp(argv[0],"eta") == 0) {
    param.setValue(eta);
    return param.addObject(4, this);
  }
  return -1;
}


int 
ElasticMaterial::updateParameter(int parameterID, Information &info)
{
  switch(parameterID) {
  case 1:
    Epos = info.theDouble;
    Eneg = info.theDouble;
    return 0;
  case 2:
    Epos = info.theDouble;
    return 0;
  case 3:
    Eneg = info.theDouble;
    return 0;
  case 4:
    eta = info.theDouble;
    return 0;
  default:
    return -1;
  }
}


int
ElasticMaterial::activateParameter(int paramID)
{
  parameterID = paramID;
  return 0;
}


double
ElasticMaterial::getStressSensitivity(int gradIndex, bool conditional)
{
  if (parameterID == 1)
    return trialStrain;
  if (parameterID == 2 && trialStrain > 0.0)
    return trialStrain;
  if (parameterID == 3 && trialStrain < 0.0)
    return trialStrain;
  if (parameterID == 4)
    return trialStrainRate;

  return 0.0;
}


double
ElasticMaterial::getTangentSensitivity(int gradIndex)
{
  if (parameterID == 1)
    return 1.0;
  if (parameterID == 2 && trialStrain >= 0.0)
    return 1.0;
  if (parameterID == 3 && trialStrain <= 0.0)
    return 1.0;

  return 0.0;
}


double
ElasticMaterial::getInitialTangentSensitivity(int gradIndex)
{
  if (parameterID == 1)
    return 1.0;
  if (parameterID == 2)
    return 1.0;
  if (parameterID == 3)
    return 1.0;

  return 0.0;
}


int
ElasticMaterial::commitSensitivity(double strainGradient,
				   int gradIndex, int numGrads)
{
  // Nothing to commit ... path independent
  return 0;
}
