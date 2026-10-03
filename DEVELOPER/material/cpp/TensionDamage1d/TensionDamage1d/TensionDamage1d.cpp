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
#include "TensionDamage1d.h"

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

static int numTensionDamage1d = 0;

OPS_Export void *
OPS_TensionDamage1d()
{
  // print out some KUDO's
  if (numTensionDamage1d == 0) {
    opserr << "TensionDamage1d unaxial material - Written by Francesco Vanin 2018 - Loaded from external library. \n";
    numTensionDamage1d =1;
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  //
  // parse the input line for the material parameters
  //

  int    iData[1];
  double dData[3];
  int numData;
  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial TensionDamage1d tag" << endln;
    return 0;
  }
 
  numData = 3;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid E, ft, and GfI \n";
    return 0;	
  }

    double secantTangent = 0.0;
    numData = 1;
    while (OPS_GetNumRemainingInputArgs() > 0) {
		OPS_GetDoubleInput(&numData,&secantTangent);
		opserr << "Read Secant tangent" << secantTangent << endln;
	}

	

  // 
  // create a new material
  //

	theMaterial = new TensionDamage1d(iData[0], dData[0], dData[1], dData[2], secantTangent);       

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type TensionDamage1d\n";
    return 0;
  }

  // return the material
  return theMaterial;
}




TensionDamage1d::TensionDamage1d(int _tag, double _E, double _ft, double _GfI, double secantTangent)
:UniaxialMaterial(_tag, 0),
 E(_E), ft(_ft), GfI(_GfI), D(0.0), Dcommit(0.0),
 trialStrain(0.0), trialStress(0.0), trialTangent(_E),
 commitStrain(0.0), commitStress(0.0), commitTangent(_E), secantTangent(secantTangent)
{ ey = ft/E;
}

TensionDamage1d::TensionDamage1d()
:UniaxialMaterial(0, 0),
 E(0.0), ft(0.0), GfI(0.0), ey(0.0), D(0.0), Dcommit(0.0),
 trialStrain(0.0), trialStress(0.0), trialTangent(0.0),
 commitStrain(0.0), commitStress(0.0), commitTangent(0.0)
{ }

TensionDamage1d::~TensionDamage1d()
{
  // does nothing
}

int 
TensionDamage1d::setTrialStrain(double strain, double strainRate)
{
	/*
    if (fabs(trialStrain - strain) < DBL_EPSILON){
		//opserr << "not Needed : " << strain  << ", stress: " << trialStress << ", tangent: " << trialTangent << endln;
      return 0;
	}
	*/
    trialStrain = strain;
	
	if (trialStrain<=0.0) {
		trialStress = E*trialStrain;
		trialTangent = E;
	} 
	
	else {
		// update damage  variable
		double k = trialStrain - ey;

		if (k>0.0) {
			D = 1.0 - (exp(-ft/GfI*k)*ey) / (ey+k);  // exponential softening
		} else {
			D = Dcommit;
		}
	
		if (D<Dcommit)  D=Dcommit;

		trialStress = (1.0-D)*E*trialStrain;

		if (D==Dcommit) {
			trialTangent = (1.0-D)*E;
		} else {
			if (Dcommit>DBL_EPSILON) {
				trialTangent = -exp(-ft/GfI*k) *ft*ft/GfI; // exponential softening

				//trialTangent = (trialStress-commitStress)/(trialStrain-commitStrain);

			} else {
				trialTangent = -exp(-ft/GfI*k) *ft*ft/GfI; // exponential softening
				trialTangent = (1.0-secantTangent)*trialTangent + secantTangent*(1.0-D)*E;  // secant/tangent stiffness

				trialTangent = (trialStress-commitStress)/(trialStrain-commitStrain);
			}

			/*
			if (abs(commitStrain-trialStrain)>DBL_EPSILON) {
				trialTangent = (commitStress-trialStress) / (commitStrain-trialStrain);
			}
			*/

			//opserr << "called set strain: " << strain  << ", stress: " << trialStress << ", tangent: " << trialTangent << endln;
		}

	}

	//opserr << "called set strain: " << strain  << ", stress: " << trialStress << ", tangent: " << trialTangent << endln;

    return 0;
}

double 
TensionDamage1d::getStrain(void)
{
  return trialStrain;
}

double 
TensionDamage1d::getStress(void)
{
  return trialStress;
}


double 
TensionDamage1d::getTangent(void)
{
	//opserr << "called get tangent: " << trialTangent << endln;
  return trialTangent;
}

int 
TensionDamage1d::commitState(void)
{
    commitStrain  = trialStrain;
    commitTangent = trialTangent;
    commitStress  = trialStress;

	Dcommit = D;

    return 0;
}	


int 
TensionDamage1d::revertToLastCommit(void)
{
  trialStrain = commitStrain;
  trialTangent = commitTangent;
  trialStress = commitStress;

  D = Dcommit;

  return 0;
}


int 
TensionDamage1d::revertToStart(void)
{
  trialStrain = commitStrain = 0.0;
  trialTangent = commitTangent = E;
  trialStress = commitStress = 0.0;

  D = Dcommit = 0.0;

  return 0;
}


UniaxialMaterial *
TensionDamage1d::getCopy(void)
{
  TensionDamage1d* theCopy = new TensionDamage1d(this->getTag(),E,ft,GfI,secantTangent);
  theCopy->Dcommit = this->Dcommit;
  
  return theCopy;
}


int 
TensionDamage1d::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(9);
  data(0) = this->getTag();
  data(1) = Dcommit;
  data(2) = E;
  data(3) = ft;
  data(4) = GfI;
  data(5) = ey;
  data(6) = commitStrain;
  data(7) = commitStress;
  data(8) = commitTangent;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ElasticPPcpp::sendSelf() - failed to send data\n";

  return res;
}

int 
TensionDamage1d::recvSelf(int cTag, Channel &theChannel, 
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
    ft = data(3);
    GfI   = data(4);
    ey   = data(5);  
    commitStrain=data(6);
    commitStress=data(7);
    commitTangent=data(8);
    trialStrain = commitStrain;
    trialTangent = commitTangent;
    trialStress = commitStress;
  }

  return res;
}

void 
TensionDamage1d::Print(OPS_Stream &s, int flag)
{
  s << "ElasticPPcpp tag: " << this->getTag() << endln;
  s << "  E: " << E << endln;
  s << "  D: " << D << endln;
  s << "  stress: " << trialStress << " tangent: " << trialTangent << endln;
}


