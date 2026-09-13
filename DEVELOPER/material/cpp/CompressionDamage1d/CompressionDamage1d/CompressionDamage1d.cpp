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
                                                                        
// $Revision: 1.7 $
// $Date: 2009/03/23 23:17:04 $
// $Source: /usr/local/cvs/OpenSees/PACKAGES/NewMaterial/cpp/ElasticPPcpp.cpp,v $
                                                                        
// Written: fmk 
//
// Description: This file contains the class implementation for 
// ElasticMaterial. 
//
// What: "@(#) ElasticPPcpp.C, revA"

#include <elementAPI.h>
#include "CompressionDamage1d.h"

#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>


#ifdef _USRDLL
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif

static int numCompressionDamage1d = 0;

OPS_Export void *
OPS_CompressionDamage1d()
{
  // print out some KUDO's
  if (numCompressionDamage1d == 0) {
    opserr << "CompressionDamage1d uniaxial material - Written by Francesco Vanin 2018 - Loaded from external library. \n";
    numCompressionDamage1d =1;
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  //
  // parse the input line for the material parameters
  //

  int    iData[1];
  double dData[2];
  int numData;
  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial CompressionDamage1d tag" << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid E, fc. \n";
    return 0;	
  }

  
  double maxDuctility = -1.0;
  double r = 1;
  bool triangular = false;
  const char* type;
  
  while (OPS_GetNumRemainingInputArgs() > 0) {
	  type = OPS_GetString();

	  if (strcmp(type, "-maxDuctility") == 0) {
		  numData = 2;
		  double dData2[2];
		  if (OPS_GetDoubleInput(&numData, dData2) < 0) {
			  opserr << "WARNING invalid maxDuctility (" << maxDuctility << ", input ignored \n";
			  maxDuctility = -1;
			  r = 1;
		  }
		  else {
			  maxDuctility = dData2[1];
			  r = dData2[0];
		  }
	  }
	  if (strcmp(type, "-triangular") == 0) {
		  triangular = true;
	  }
  }

  // 
  // create a new material
  //


  theMaterial = new CompressionDamage1d(iData[0], dData[0], dData[1], maxDuctility, r, triangular);       

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type CompressionDamage1d\n";
    return 0;
  }

  // return the material
  return theMaterial;
}




CompressionDamage1d::CompressionDamage1d(int _tag, double _E, double _fc, double _maxDuctility, double r, bool triangular)
:UniaxialMaterial(_tag, 0),
E(_E), fc(_fc), D(0.0), Dcommit(0.0), maxDuctility(_maxDuctility), r(r), triangular(triangular),
 trialStrain(0.0), trialStress(0.0), trialTangent(_E),
 commitStrain(0.0), commitStress(0.0), commitTangent(_E)
{ ey = fc/E;
}

CompressionDamage1d::CompressionDamage1d()
:UniaxialMaterial(0, 0),
E(0.0), fc(0.0), ey(0.0), D(0.0), Dcommit(0.0), maxDuctility(0.0), r(0.0), triangular(false),
 trialStrain(0.0), trialStress(0.0), trialTangent(0.0),
 commitStrain(0.0), commitStress(0.0), commitTangent(0.0)
{ }

CompressionDamage1d::~CompressionDamage1d()
{
  // does nothing
}

int 
CompressionDamage1d::setTrialStrain(double strain, double strainRate)
{
	
    if (fabs(trialStrain - strain) < DBL_EPSILON){
		//opserr << "not Needed : " << strain  << ", stress: " << trialStress << ", tangent: " << trialTangent << endln;
      return 0;
	}
	
    trialStrain = strain;
	
	if (trialStrain>=0.0) {
		trialStress = 0.0;
		trialTangent = 0.0001*E;
	} 	
	else {
		// ductility demand
		double mu = -strain/ey;

		// damage (perfect plasticity)
		if (mu<1.0) {
			D = Dcommit;
		} else {
			if (mu>maxDuctility && maxDuctility>1) {
				//D = 1 - r*(1 - D);
				D = 1 - r / mu;

			} else {
				if (triangular && maxDuctility>1) {
					D = ((mu-1)*(maxDuctility -r))/(mu*(maxDuctility-1));
				}
				else {
					D = 1. - 1. / mu;
				}
			}
		}
		
		if (D<Dcommit)  D=Dcommit;

		trialStress = (1.0-D)*E*trialStrain;

		if (D==Dcommit) {
			trialTangent = (1.0-D)*E;
		} else {
			trialTangent = 0.0001*E;
		}

	}
    return 0;
}

double 
CompressionDamage1d::getStrain(void)
{
  return trialStrain;
}

double 
CompressionDamage1d::getStress(void)
{
  return trialStress;
}


double 
CompressionDamage1d::getTangent(void)
{
  return trialTangent;
}

double 
CompressionDamage1d::getInitialTangent(void) {
   return E;
}

int 
CompressionDamage1d::commitState(void)
{
    commitStrain  = trialStrain;
    commitTangent = trialTangent;
    commitStress  = trialStress;

	Dcommit = D;

    return 0;
}	


int 
CompressionDamage1d::revertToLastCommit(void)
{
  trialStrain = commitStrain;
  trialTangent = commitTangent;
  trialStress = commitStress;

  D = Dcommit;

  return 0;
}


int 
CompressionDamage1d::revertToStart(void)
{
  trialStrain = commitStrain = 0.0;
  trialTangent = commitTangent = E;
  trialStress = commitStress = 0.0;

  D = Dcommit = 0.0;

  return 0;
}


UniaxialMaterial *
CompressionDamage1d::getCopy(void)
{
  CompressionDamage1d* theCopy = new CompressionDamage1d(this->getTag(),E,fc,maxDuctility,r,triangular);
  theCopy->Dcommit = this->Dcommit;
  
  return theCopy;
}


int 
CompressionDamage1d::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(8);
  data(0) = this->getTag();
  data(1) = Dcommit;
  data(2) = E;
  data(3) = fc;
  data(4) = ey;
  data(5) = commitStrain;
  data(6) = commitStress;
  data(7) = commitTangent;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ElasticPPcpp::sendSelf() - failed to send data\n";

  return res;
}

int 
CompressionDamage1d::recvSelf(int cTag, Channel &theChannel, 
				 FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(9);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ElasticPPcpp::recvSelf() - failed to recv data\n";
  else {
    this->setTag(data(0));
    Dcommit    = data(1);
    E     = data(2);
    fc = data(3);
    ey   = data(4);  
    commitStrain=data(5);
    commitStress=data(6);
    commitTangent=data(7);

    trialStrain = commitStrain;
    trialTangent = commitTangent;
    trialStress = commitStress;
  }

  return res;
}

void 
CompressionDamage1d::Print(OPS_Stream &s, int flag)
{
  s << "ElasticPPcpp tag: " << this->getTag() << endln;
  s << "  E: " << E << endln;
  s << "  D: " << D << endln;
  s << "  stress: " << trialStress << " tangent: " << trialTangent << endln;
}


