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
#include "BeamFrictionSupport.h"

#include <Channel.h>
#include <math.h>
#include <float.h>
#include <Information.h>
#include <Channel.h>
#include <NDMaterial.h>
#include <MaterialResponse.h>


#ifdef _USRDLL
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif

static int numBeamFrictionSupport = 0;

OPS_Export void *
OPS_BeamFrictionSupport()
{
  // print out some KUDO's
  if (numBeamFrictionSupport == 0) {
    opserr << "BeamFrictionSupport, NdMaterial - Written by Francesco Vanin - EPFL, 2018 - v 1.1 \n";
    numBeamFrictionSupport =1;
  }

  // Pointer to a NDmaterial that will be returned
  NDMaterial *theMaterial = 0;

   // parse the input line for the material parameters
  //
  int    iData[1];
  double dData[3];
  int numData;
  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid NDMaterial BeamFrictionSupport tag" << endln;
    return 0;
  }

  numData = 3;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING BeamFrictionSupport::input structure has to be $E, $G, $mu <-anchoring $dx> <-gap $gapX, $gapY>. \n";
    return 0;	
  }

  double E  = dData[0];
  double G  = dData[1];
  double mu = dData[2];

  double dx, gapX, gapY;
  dx = gapX = gapY = -1.;

  // options for input
  while(OPS_GetNumRemainingInputArgs() > 0) {
	const char* type = OPS_GetString();
    if (strcmp(type,"-anchoring") == 0) {
	    if (OPS_GetNumRemainingInputArgs() > 0) {
			numData = 1;
		    if (OPS_GetDoubleInput(&numData,&dx) < 0) {
		        opserr<<"WARNING BeamFrictionSupport:: invalid anchoring length, assumes infinite sliding\n";
		        dx = -1.;
		    }
	    }
	} else if (strcmp(type,"-gap") == 0) {
             numData = 2;
			 double dData2[2];
			 if (OPS_GetDoubleInput(&numData,&dData2[0]) < 0) {
		        opserr<<"WARNING BeamFrictionSupport:: invalid gap, assumes infinite sliding\n";
		        gapX = gapY = -1.;
			 } else {
				 gapX = dData2[0];
				 gapY = dData2[1];
			 }
		   }
   }
   
  // 
  // create the material
  //

  theMaterial = new BeamFrictionSupport(iData[0], E, G, mu, dx, gapX, gapY);       

  if (theMaterial == 0) {
    opserr << "WARNING could not create NDMaterial of type BeamFrictionSupport\n";
    return 0;
  }

  // return the material
  return theMaterial;
}


// full constructor
BeamFrictionSupport::BeamFrictionSupport(int _tag, double _E, double _G, double _mu, double _anchoringLength, double _gapX, double _gapY)
:NDMaterial(_tag, 0), 
mu(_mu), dx(_anchoringLength), gapX(_gapX), gapY(_gapY), failed(false), failedTrial(false),
slidingTrial(3), sliding(3), e(3), eTrial(3), s(3), sTrial(3), D(3,3), Dtrial(3,3), Dcommitted(3,3) 
{
   D(0,0) = _E;
   D(1,1) = D(2,2) = _G;
   Dtrial = D;
   Dcommitted = D;

   fakeStiffnessFactor =  1.0e-6; 
}

// null constructor
BeamFrictionSupport::BeamFrictionSupport()
:NDMaterial(0, 0), 
mu(0.0), dx(0.0), gapX(0.0), gapY(0.0), failed(false),  failedTrial(false),
slidingTrial(3), sliding(3), e(3), eTrial(3), s(3), sTrial(3), D(3,3), Dtrial(3,3), Dcommitted(3,3)   {
	fakeStiffnessFactor = 0.0;
}

// destructor
BeamFrictionSupport::~BeamFrictionSupport()
{
  // does nothing, no object to delete was created
}

//set strain methods
int 
BeamFrictionSupport::setTrialStrain(const Vector& strain, const Vector& strainRate)
{
	return this->setTrialStrain(strain);
}

int 
BeamFrictionSupport::setTrialStrain(const Vector& strain)
{
	Vector increment(3);
	increment = eTrial - strain;
	
    if (fabs(increment.Norm()) < DBL_EPSILON) {
      return 0;
	}

	//opserr << "entered. \n ";
    eTrial = strain;
	slidingTrial = sliding;
	increment = eTrial - e;

	if (failed) {
		sTrial.Zero();
		Dtrial.Zero();
		return 0;
	}
	
	// check for eccessive sliding in the beam direction
	if ((eTrial(1)>dx) && (dx>0.0)) {
		sTrial.Zero();
		Dtrial.Zero();
		Dtrial(0,0) = 0.0*D(0,0);

		Dtrial(1,1) = 0.0*D(1,1);
		Dtrial(2,2) = 0.0*D(2,2);

		failedTrial = true;
		return 0;
	}
	
	// vertical stress
	sTrial(0) = D(0,0)*eTrial(0);

	// separate simple 2d case from 3d
	if (gapY == 0.0 && gapX == 0.0) {
		
		sTrial(2) = D(2, 2)*eTrial(2);
		Dtrial(2, 2) = D(2, 2);

		//elastic in y direction , possible sliding only in positive x direction
		if (sTrial(0)>DBL_EPSILON) {
			sTrial.Zero();
			Dtrial = fakeStiffnessFactor*D;

			sTrial(2) = D(2, 2)*eTrial(2);
		    Dtrial(2, 2) = D(2, 2);

			if (eTrial(1)<0.0) { 
				sTrial(1) = D(1, 1)*eTrial(1);
				Dtrial(1, 1) = D(1, 1);
			}
			else {
				slidingTrial(1) += increment(1);
				//slidingTrial(1) = eTrial(1);
			}
			return 0;
		}

		// else in compression: check sliding condition
		double tauTrial = D(1, 1)* (eTrial(1) - slidingTrial(1));
		double f2 = mu*sTrial(0) + abs(tauTrial);
		bool correctedForNegative = false;
		

		double signTau = 1.0;
		if (tauTrial < 0.0) { signTau = -1.0; }
		
		double deltaLambda = 0.;

		if (f2 <= DBL_EPSILON) {
			// no sliding occurs
			sTrial = D*(eTrial - sliding);
			Dtrial = D;

		}
		else {
			
			deltaLambda = f2 / D(1, 1);
			slidingTrial(1) += deltaLambda * signTau;

			if (slidingTrial(1)<0.0) {
				slidingTrial(1) = 0.0;
				correctedForNegative = true;
			}

			sTrial = D*(eTrial - slidingTrial);
			
			Dtrial.Zero();
			Dtrial(0, 0) = D(0, 0);
			Dtrial(1, 0) = -D(0, 0)*mu*signTau;
			Dtrial(1, 1) = D(1, 1)*fakeStiffnessFactor;
			Dtrial(2, 2) = D(2, 2);

			if (correctedForNegative) {
				Dtrial(1, 0) = 0.0;
				Dtrial(1, 1) = D(1, 1);
			}
		}

		return 0;
	}

	if (sTrial(0)>DBL_EPSILON) {
		sTrial.Zero();

		slidingTrial(1) += increment(1);
		slidingTrial(2) += increment(2);

		// no stiffness provided at the interface level
		Dtrial = fakeStiffnessFactor*D;
		return 0;
	}

	// free sliding
	double tau1 = D(1,1)* ( eTrial(1) - slidingTrial(1) );
	double tau2 = D(2,2)* ( eTrial(2) - slidingTrial(2) );
	double tauTrial = sqrt(tau1*tau1 + tau2*tau2);
	
	double f2 = mu*sTrial(0) + tauTrial;
	double deltaLambda = 0.;

	if (f2<=DBL_EPSILON) {
		// no sliding occurs
		sTrial = D*(eTrial-sliding);
		Dtrial = D;

	} else {
		deltaLambda = f2 / D(1,1);
		slidingTrial(1) += deltaLambda * tau1 / tauTrial;
		slidingTrial(2) += deltaLambda * tau2 / tauTrial;

		sTrial = D*(eTrial-slidingTrial);
		Dtrial.Zero();

		double fy = -mu*sTrial(0);

		Dtrial(0,0) = D(0,0);
		Dtrial(1,1) = (fy * tau2*tau2 / pow(tauTrial,3.)  + fakeStiffnessFactor) *D(1,1);
		Dtrial(2,2) = (fy * tau1*tau1 / pow(tauTrial,3.)  + fakeStiffnessFactor) *D(2,2);
		Dtrial(1,2) = Dtrial(2,1) = -fy * tau1*tau2 /  pow(tauTrial,3.) *D(1,1);

	}

    // check for closing gaps
	if ((-eTrial(1)>gapX) && (gapX>=0.0)) {
		sTrial(1) -= (abs(eTrial(1)) - gapX) *D(1,1);
		Dtrial(1,1) = D(1,1);
	}

	if ((abs(eTrial(2))>gapY) && (gapY>=0.0)) {
		if (-eTrial(2) > gapY) {
			sTrial(2) -= (abs(eTrial(2)) - gapY) *D(2,2);
		} else {
			sTrial(2) += (abs(eTrial(2)) - gapY) *D(2,2);
		}

		Dtrial(2,2) = D(2,2);
	}

    return 0;
}


int
BeamFrictionSupport::setTrialStrainIncr (const Vector &strain) {

	return this->setTrialStrain(e + strain);
}


int
BeamFrictionSupport::setTrialStrainIncr (const Vector &strain, const Vector &rate) {  
	return this->setTrialStrainIncr(strain);
}



// get strain, stress and tangents methods
const Matrix& 
BeamFrictionSupport::getTangent(void) {
	return Dtrial;
}

const Matrix& 
BeamFrictionSupport::getInitialTangent(void) {
	return D;
 }

const Matrix& 
BeamFrictionSupport::getCommittedTangent(void) {
	return Dcommitted;
 }


const Vector& 
BeamFrictionSupport::getStress(void) {
	 return sTrial;
 }

const Vector& 
BeamFrictionSupport::getStrain(void) {
	return eTrial;
}

// commit variables methods
int 
BeamFrictionSupport::commitState(void) {

	// update state variables
	e = eTrial;
	s = sTrial;
	Dcommitted = Dtrial;
	sliding = slidingTrial;
	if ((!failed) && failedTrial)
		opserr << "WARNING: instance of material " << this->getTag() << " failed. \n";

	if (sTrial(0)>DBL_EPSILON) 
		opserr << "committing an open interface\n";

	failed = failedTrial;
	
	return 0;
}

int 
BeamFrictionSupport::revertToLastCommit(void) {
	eTrial = e;
	sTrial = s;
	Dtrial = Dcommitted;
	slidingTrial = sliding;
	failedTrial = failed;
	
	return 0;
}

int 
BeamFrictionSupport::revertToStart(void) {
	e.Zero();
	s.Zero();
	Dcommitted = D;
	sliding.Zero();

    eTrial = e;
	sTrial = s;
	Dtrial = Dcommitted;
	slidingTrial = sliding;

	failed = false;
	failedTrial = false;
	
	return 0;
}



// get copies methods
NDMaterial* 
BeamFrictionSupport::getCopy(void) {
	BeamFrictionSupport* theCopy;
	theCopy = new BeamFrictionSupport(this->getTag(), D(0,0), D(1,1), mu, dx, gapX, gapY);
	return theCopy;
};


NDMaterial* 
BeamFrictionSupport::getCopy(const char *code) {
	// meaningful only in order 3 version
	return this->getCopy();
}


const char *
BeamFrictionSupport::getType(void) const {
	return "BeamFiber";
}

int 
BeamFrictionSupport::getOrder(void) const {
	return 3;
}



// response methods
Response *
BeamFrictionSupport::setResponse(const char **argv, int argc, OPS_Stream &output) {
     Response *theResponse = 0;
     output.tag("NdMaterialOutput");
     output.attr("matType",this->getClassType());
     output.attr("matTag",this->getTag());

  if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0) {
    const Vector &res = this->getStress();
    int size = res.Size();
    output.tag("ResponseType","tauX");
	output.tag("ResponseType","tauY");   
	output.tag("ResponseType","sigmaZ"); 
    theResponse =  new MaterialResponse(this, 1, this->getStress());

  } else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0) {
    const Vector &res = this->getStrain();
    int size = res.Size();
	output.tag("ResponseType","gammaX");
	output.tag("ResponseType","gammaY");
	output.tag("ResponseType","epsilonZ");
    theResponse =  new MaterialResponse(this, 2, this->getStrain());

  } else if (strcmp(argv[0],"tangent") == 0 || strcmp(argv[0],"stiffness") == 0) {
	const Matrix &res = this->getCommittedTangent();
    int size = 9;
	output.tag("ResponseType","K11");
	output.tag("ResponseType","K12");
	output.tag("ResponseType","K13");
	output.tag("ResponseType","K21");
	output.tag("ResponseType","K22");
	output.tag("ResponseType","K23");
	output.tag("ResponseType","K31");
	output.tag("ResponseType","K32");
	output.tag("ResponseType","K33");
    theResponse =  new MaterialResponse(this, 3, this->getCommittedTangent());
  
  } else if (strcmp(argv[0],"slip") == 0 || strcmp(argv[0],"sliding") == 0) {
	Vector &res = sliding;
    int size = 3;
	output.tag("ResponseType","slipX");
	output.tag("ResponseType","slipY");
	output.tag("ResponseType","zero");
    theResponse =  new MaterialResponse(this, 4, sliding);
  }

  output.endTag(); // NdMaterialOutput

  return theResponse;
}
 
int 
BeamFrictionSupport::getResponse(int responseID, Information &matInfo) {

  switch (responseID) {
  case 1:
    return matInfo.setVector(this->getStress());
    
  case 2:
    return matInfo.setVector(this->getStrain());

  case 3:
	  return matInfo.setMatrix(this->getCommittedTangent());

  case 4:
	  return matInfo.setVector(sliding);

  default:
    return -1;
  }
}



// communication methods
void 
BeamFrictionSupport::Print( OPS_Stream &s, int flag ) {
  s << endln ;
  s << "BeamFrictionSupport : " ; 
  s << "tag: " << this->getTag() << endln ;
  s << "E:   " << D(0,0) << endln ;
  s << "G:   " << D(1,1) << endln ;
  s << "mu:  " << mu     << endln ;
  s << "dx:  " << dx     << endln ;
  s << "gapX:" << gapX   << endln ;
  s << "gapY:" << gapY   << endln ;
  if (failed)
	  s << "failed" << endln ;
  else
	  s << "not failed" << endln ;
  s << endln ;

}

int  
BeamFrictionSupport::sendSelf(int commitTag, Channel &theChannel) 	{ 
	opserr << "BeamFrictionSupport::sendSelf not yet implemented. \n";
	return -1;  
}

int  
BeamFrictionSupport::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) { 
	opserr << "BeamFrictionSupport::recvSelf not yet implemented. \n";
	return -1;  
}