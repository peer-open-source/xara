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

// $Revision: 1.0 $
// $Date: 2012-05-26 22:03:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/PlateRebarMaterial.cpp,v $

//
// Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)
//
// Generic Plate Rebar Material
//
// Ref: Lu X, Lu XZ, Guan H, Ye LP, Collapse simulation of reinforced 
//      concrete high-rise building induced by extreme earthquakes, 
//      Earthquake Engineering & Structural Dynamics, 2013, 42(5): 705-723
//
#include <PlateRebarMaterial.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MaterialResponse.h>
#include <math.h>

Vector PlateRebarMaterial::stress(5);
Matrix PlateRebarMaterial::tangent(5, 5);

PlateRebarMaterial::PlateRebarMaterial() : NDMaterial(0, ND_TAG_PlateRebarMaterial), strain(5) {}


PlateRebarMaterial::PlateRebarMaterial(int tag, UniaxialMaterial& uniMat, double ang)
 : NDMaterial(tag, ND_TAG_PlateRebarMaterial), strain(5), angle(ang)
{
  theMat      = uniMat.getCopy();
  double rang = ang * 4.0 * asin(1.0) / 360.0;
  c           = cos(rang);
  s           = sin(rang);
}


PlateRebarMaterial::~PlateRebarMaterial()
{
  if (theMat != nullptr)
    delete theMat;
}


NDMaterial*
PlateRebarMaterial::getCopy()
{
  return new PlateRebarMaterial(this->getTag(), *theMat,angle);
}


NDMaterial*
PlateRebarMaterial::getCopy(const char* type)
{
  if (strcmp(type, this->getType()) == 0)
    return this->getCopy();

  return this->NDMaterial::getCopy(type);
}


int
PlateRebarMaterial::getOrder() const
{
  return 5;
}


const char*
PlateRebarMaterial::getType() const
{
  return "PlateFiber";
}


int
PlateRebarMaterial::commitState()
{
  return theMat->commitState();
}


int
PlateRebarMaterial::revertToLastCommit()
{
  return theMat->revertToLastCommit();
}


int
PlateRebarMaterial::revertToStart()
{
  strain.Zero();
  return theMat->revertToStart();
}


//mass per unit volume
double
PlateRebarMaterial::getRho()
{
  return theMat->getRho();
}


//receive the strain
int
PlateRebarMaterial::setTrialStrain(const Vector& strainFromElement)
{
  strain(0) = strainFromElement(0);
  strain(1) = strainFromElement(1);
  strain(2) = strainFromElement(2);
  strain(3) = strainFromElement(3);
  strain(4) = strainFromElement(4);

  if (angle == 0)
    return theMat->setTrialStrain(strain(0));
  else if (angle == 90)
    return theMat->setTrialStrain(strain(1));

  return theMat->setTrialStrain(strain(0) * c * c + strain(1) * s * s + strain(2) * c * s, 0);
}


const Vector&
PlateRebarMaterial::getStrain()
{
  return strain;
}


const Vector&
PlateRebarMaterial::getStress()
{
  double sig = theMat->getStress();

  stress.Zero();
  if (angle == 0)
    stress(0) = sig;

  else if (angle == 90)
    stress(1) = sig;

  else {
    stress(0) = sig * c * c;
    stress(1) = sig * s * s;
    stress(2) = sig * c * s;
  }

  return stress;
}


const Matrix&
PlateRebarMaterial::getTangent()
{
  double tan = theMat->getTangent();

  tangent.Zero();
  if (angle == 0)
    tangent(0, 0) = tan;
  else if (angle == 90)
    tangent(1, 1) = tan;
  else {
    tangent(0, 0) = tan * c * c * c * c;
    tangent(0, 1) = tan * c * c * c * s;
    tangent(0, 2) = tan * c * c * s * s;
    tangent(1, 0) = tangent(0, 1);
    tangent(1, 1) = tangent(0, 2);
    tangent(1, 2) = tan * c * s * s * s;
    tangent(2, 0) = tangent(0, 2);
    tangent(2, 1) = tangent(1, 2);
    tangent(2, 2) = tan * s * s * s * s;
  }

  return tangent;
}

const Matrix&
PlateRebarMaterial::getInitialTangent()
{
  double tan = theMat->getInitialTangent();
  tangent.Zero();
  if (angle == 0)
    tangent(0, 0) = tan;
  else if (angle == 90)
    tangent(1, 1) = tan;
  else {
    tangent(0, 0) = tan * c * c * c * c;
    tangent(0, 1) = tan * c * c * c * s;
    tangent(0, 2) = tan * c * c * s * s;
    tangent(1, 0) = tangent(0, 1);
    tangent(1, 1) = tangent(0, 2);
    tangent(1, 2) = tan * c * s * s * s;
    tangent(2, 0) = tangent(0, 2);
    tangent(2, 1) = tangent(1, 2);
    tangent(2, 2) = tan * s * s * s * s;
  }

  return tangent;
}


void
PlateRebarMaterial::Print(OPS_Stream& s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
    s << "PlateRebar Material tag: " << this->getTag() << endln;
    s << "angle: " << angle << endln;
    s << "using uniaxial material: " << endln;
    theMat->Print(s, flag);
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"PlateRebarMaterial\", ";
    s << "\"angle\": " << angle << ", ";
    s << "\"material\": \"" << theMat->getTag() << "\"}";
  }
}


int
PlateRebarMaterial::sendSelf(int commitTag, Channel& theChannel)
{
  int res = 0;

  int dataTag = this->getDbTag();

  int matDbTag;

  static ID idData(3);
  idData(0) = dataTag;
  idData(1) = theMat->getClassTag();
  matDbTag  = theMat->getDbTag();
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    theMat->setDbTag(matDbTag);
  }
  idData(2) = matDbTag;

  res = theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "PlateRebarMaterial::sendSelf() - failed to send data" << endln;
    return res;
  }

  static Vector vecData(1);
  vecData(0) = angle;

  res = theChannel.sendVector(dataTag, commitTag, vecData);
  if (res < 0) {
    opserr << "PlateRebarMaterial::sendSelf() - failed to send data" << endln;
    return res;
  }

  // now send the materials data
  res += theMat->sendSelf(commitTag, theChannel);
  if (res < 0)
    opserr << "PlateRebarMaterial::sendSelf() - failed to send material1" << endln;

  return res;
}

int
PlateRebarMaterial::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
  int res = 0;

  int dataTag = this->getDbTag();

  // recv an id containing the tag and associated materials class and db tags
  static ID idData(3);
  res = theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "PlateRebarMaterial::sendSelf() - failed to receive id data" << endln;
    return res;
  }

  this->setTag(idData(0));
  int matClassTag = idData(1);
  if (theMat->getClassTag() != matClassTag) {
    if (theMat != 0)
      delete theMat;
    theMat = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theMat == 0) {
      opserr << "PlateRebarMaterial::recvSelf() - failed to get a material of type: " << matClassTag
             << endln;
      return -1;
    }
  }
  theMat->setDbTag(idData(2));

  static Vector vecData(1);
  res = theChannel.recvVector(dataTag, commitTag, vecData);
  if (res < 0) {
    opserr << "PlateRebarMaterial::sendSelf() - failed to receive vector data" << endln;
    return res;
  }
  angle       = vecData(0);
  double rang = angle * 0.0174532925;
  c           = cos(rang);
  s           = sin(rang);

  // now receive the materials data
  res = theMat->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0)
    opserr << "PlateRebarMaterial::sendSelf() - failed to receive material1" << endln;

  return res;
}
