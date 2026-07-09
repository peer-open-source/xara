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
Source: /DEVELOPER/element/cpp/Macroelement3d/Macroelement3d/FrictionCohesion3d.cpp
Written by Francesco Vanin (francesco.vanin@epfl.ch)
Ecole Polytechnique Federale de Lausanne, Switzerland,
Earthquake Engineering and Structural Dynamics laboratory, 2019
*/


#include <elementAPI.h>
#include "FrictionCohesion3d.h"

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

static int numFrictionCohesion3d = 0;

OPS_Export void *
OPS_FrictionCohesion3d()
{
  // print out some KUDO's
  if (numFrictionCohesion3d == 0) {
    opserr << "FrictionCohesion3d - Loaded from external library. Written by Igor Tomic, EPFL, 2019.\n";
    numFrictionCohesion3d =1;
  }

  // Pointer to the nD material that will be returned
  NDMaterial *theMaterial = 0;
  bool imposedSigma = false;
  int numData = 0;
  double sigma0 = 0.0;


  int iData[1];
  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[0]) < 0) {
	  opserr << "WARNING: invalid tag. Required structure: ndMaterial FrictionCohesion3d $tag $G $mu $c $GfII <-flags> " << endln;
	  return 0;
  }


  // inputs: 
  double dData[4];
  numData = 4;
  if (OPS_GetDoubleInput(&numData, &dData[0]) < 0) {
	  opserr << "WARNING: invalid double inputs. Required structure: ndMaterial FrictionCohesion3d $tag $G $mu $c $GfII<-flags> " << endln;
	  return 0;
  }
  

  const char* inputStructure = OPS_GetString();

  if (strcmp(inputStructure, "-imposedSigma") == 0) {
	  double dData2[1];
	  numData = 1;
	  if (OPS_GetDoubleInput(&numData, &dData2[0]) < 0) {
		  opserr << "WARNING: invalid double inputs. Required structure: ndMaterial FrictionCohesion3d $tag $G $mu $c $GfII -imposedSigma $sigma0 " << endln;
		  return 0;
	  }
	  imposedSigma = true;
	  sigma0 = dData2[0];

	  theMaterial = new FrictionCohesion3d(iData[0], dData[0], dData[1], dData[2], dData[3], sigma0);
	  return theMaterial;
  }
  else {
	  if (strcmp(inputStructure, "-axialMaterial") == 0) {
		  int iData2[1];
		  numData = 1;
		  if (OPS_GetIntInput(&numData, &iData2[0]) < 0) {
			  opserr << "WARNING: invalid integer inputs. Required structure: ndMaterial FrictionCohesion3d $tag $G $mu $c $GfII -axialMaterial $axialMaterialTag " << endln;
			  return 0;
		  }

		  UniaxialMaterial* theAxialModel;
		  theAxialModel = OPS_GetUniaxialMaterial(iData2[0]);
		  if (theAxialModel == 0) {
			  opserr << "WARNING: Material " << iData2[0] <<  " not found." << endln;
			  return 0;
		  }

		  theMaterial = new FrictionCohesion3d(iData[0], theAxialModel, dData[0], dData[1], dData[2], dData[3]);
		  return theMaterial;
	  }
  }
 
  if (theMaterial == 0) {
    opserr << "WARNING could not create NDMaterial of type FrictionCohesion3d\n";
    return 0;
  }
  
}




// full constructors
FrictionCohesion3d::FrictionCohesion3d(int tag, UniaxialMaterial* comprMat, double G, double mu, double c, double GfII)
	:NDMaterial(tag, 0), stress(3), stressCommitted(3), u(3), uCommitted(3), up(3), upCommitted(3), Kpen(3,3), K(3,3), KCommitted(3,3), theAxialMat(0), mu(mu), c(c), G(G), GfII(GfII),
	k(0.0), kCommitted(0.0), deltaLambda(0.0), imposedSigma(false), sigma0(0.0)
{ 
	// create a copy of the axial model
	theAxialMat = comprMat->getCopy();

	this->E = theAxialMat->getInitialTangent();

    Kpen(0,0) = E;
    Kpen(1,1) = Kpen(2,2) = G;

    K = Kpen;
    KCommitted = Kpen;

	double GfII_lim = c*c / G;
	if (this->GfII < 1.1*GfII_lim) {
		this->GfII = 1.1*GfII_lim;
		opserr << "FrictionCohesion3d::tag:" << tag << ": corrected GfII (fracture energy in modeII) to " << this->GfII << " to avoid snap-back at material level." << endln;
	}

}


FrictionCohesion3d::FrictionCohesion3d(int tag, double G, double mu, double c, double GfII, double sigma0)
	:NDMaterial(tag, 0), stress(3), stressCommitted(3), u(3), uCommitted(3), up(3), upCommitted(3), Kpen(3, 3), K(3, 3), KCommitted(3, 3), theAxialMat(0), mu(mu), c(c), G(G), GfII(GfII),
	k(0.0), kCommitted(0.0),deltaLambda(0.0), imposedSigma(true), sigma0(sigma0)
{
	
	Kpen(0, 0) = 0.0;
	Kpen(1, 1) = Kpen(2, 2) = G;

	K = Kpen;
	KCommitted = Kpen;

	if (sigma0 > 0.0) {
		sigma0 = -sigma0;
		opserr << "FrictionCohesion3d::tag:" << tag << ": the constant vertical stress sigma0 is assumed to act in compression." << endln;
	}

	double GfII_lim = c*c / G;
	if (this->GfII < 1.1*GfII_lim) {
		this->GfII = 1.1*GfII_lim;
		opserr << "FrictionCohesion3d::tag:" << tag << ": corrected GfII (fracture energy in modeII) to " << this->GfII << " to avoid snap-back at material level." << endln;
	}

}


FrictionCohesion3d::FrictionCohesion3d()
:NDMaterial(0, 0), stress(3), stressCommitted(3), u(3), uCommitted(3), up(3), upCommitted(3), Kpen(3, 3), K(3, 3), KCommitted(3, 3), theAxialMat(0), mu(mu), c(0.0), G(0.0), GfII(0.0),
k(0.0), kCommitted(0.0), deltaLambda(0.0), imposedSigma(false), sigma0(0.0)
{ 
}


FrictionCohesion3d::~FrictionCohesion3d() {
	if (theAxialMat != 0)
		delete theAxialMat;
}


double 
FrictionCohesion3d::hardeningLaw(double kk) {
	return c*exp(-c/GfII*kk);
}

double
FrictionCohesion3d::der_hardeningLaw(double kk) {
	return -1.0*c*c/GfII *exp(-c / GfII*kk);
}

// standard methods
int 
FrictionCohesion3d::setTrialStrain(const Vector &strain, const Vector &rate) {
	return this->setTrialStrain(strain);
}

int
FrictionCohesion3d::setTrialStrain(const Vector &strain) { 
	
	Vector increment(3);
	increment = strain - u;
	if (increment.Norm() < DBL_EPSILON)
		return 0;
	
	u = strain;
	double axialStiffness = 0.;
	
	if (!imposedSigma) {
		theAxialMat->setTrialStrain(u(0));
		sigma0         = theAxialMat->getStress();
		axialStiffness = theAxialMat->getTangent();
	}
	else {
		axialStiffness = 0.0;
	}



	// shear problem
	// --------------------------------------------------------------------
	Vector sigmaTrial(3);
	double f;

	// elastic update
	sigmaTrial = Kpen*(u-upCommitted);
	f = sqrt(pow(sigmaTrial(1), 2.) + pow(sigmaTrial(2), 2.)) + mu*sigma0 - this->hardeningLaw(kCommitted);
	
	double tol = G*DBL_EPSILON;


	if (f <= tol) {
		// elastic step
		stress = sigmaTrial;
		K.Zero();
		K(0, 0) = axialStiffness;
		K(1, 1) = K(2, 2) = G;
	}
	else {
		int maxIter = 50;
		
		int iter = 0;
		deltaLambda = 0.0;
		k = kCommitted;
		stress = sigmaTrial;

		double ay = sigmaTrial(1) / sqrt(pow(sigmaTrial(1), 2.) + pow(sigmaTrial(2), 2.));
		double az = sigmaTrial(2) / sqrt(pow(sigmaTrial(1), 2.) + pow(sigmaTrial(2), 2.));
		double J = 0.0;

		//opserr << "RETURN MAPPING --------------" << endln;
		//opserr << "iter\tdeltaLambda \tf \t\stress \tcohesion \tup \tk" << endln;

		while (abs(f) > tol && iter < maxIter) {
			iter += 1;
			
			/*
			opserr << iter << "\t";
			opserr << deltaLambda << "\t";
			opserr << f << "\t";
			opserr << stress(1) << "," << stress(2) << "\t";
			opserr << this->hardeningLaw(k) << "\t";
			opserr << up(1) << "," << up(2) << "\t";
			opserr << k << "\n";
			*/

			J = -(ay*ay+az*az)*G	- this->der_hardeningLaw(k);

			deltaLambda -= f / J;

			k = kCommitted + deltaLambda;

			stress(1) = sigmaTrial(1) - G *deltaLambda*ay;
			stress(2) = sigmaTrial(2) - G *deltaLambda*az;

			up(1) = upCommitted(1) + deltaLambda*ay;
			up(2) = upCommitted(2) + deltaLambda*az;

			f = sqrt(pow(stress(1), 2.) + pow(stress(2), 2.)) + mu*sigma0 - this->hardeningLaw(k);

		}

		if (iter >= maxIter) {
			opserr << "FrictionCohesion3d: did not converge in the return mapping." << endln;
			return -1;
		}

		// tangent operator
		K = Kpen;

		K(1, 0) += mu*axialStiffness*G*ay / J;
		K(2, 0) += mu*axialStiffness*G*az / J;

		K(1, 1) += G*G*ay*ay / J;
		K(2, 1) += G*G*ay*az / J;

		K(1, 2) += G*G*ay*az / J;
		K(2, 2) += G*G*az*az / J;
	}


	if (imposedSigma) {
		stress(0) = 0.0;
	}

    return 0;
}


int
FrictionCohesion3d::setTrialStrainIncr (const Vector &strain)
{
	opserr << "WARNING: FrictionCohesion3d has no method setTrialStrainIncr implemented" << endln; 
	return 0;
}

int
FrictionCohesion3d::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  return this->setTrialStrainIncr(strain);
}

const Matrix& 
FrictionCohesion3d::getCommittedTangent( )  {
	return KCommitted;
}
	
const Matrix& 
FrictionCohesion3d::getTangent( ) 
{
	// use a fake small stiffness if zero stiffness
	double factor = 0.0000000001;

	if (abs(K(1, 1) / G) < factor) {
		K(1, 1) = factor * G;
	}

	if (abs(K(2, 2) / G) < factor) {
		K(2, 2) = factor * G;
	}
	

	return K;
}

const Matrix&
FrictionCohesion3d::getInitialTangent (void)
{ 
	return Kpen;
}

const Vector&
FrictionCohesion3d::getStress (void)
{
	return stress;
}

const Vector&
FrictionCohesion3d::getStrain (void)
{
  return u;
}

int
FrictionCohesion3d::commitState (void)
{
	if (theAxialMat !=0)
		theAxialMat->commitState();

	KCommitted = this->getTangent();
	uCommitted = u;
	upCommitted = up;
	kCommitted = k;
	stressCommitted = this->getStress();

	return 0;
}

int
FrictionCohesion3d::revertToLastCommit (void)     
{ 
	if (theAxialMat != 0)
		theAxialMat->revertToLastCommit();

	K = KCommitted;
	u = uCommitted;
	up = upCommitted;
	k = kCommitted;
	stress = stressCommitted;

  return 0;
}

int
FrictionCohesion3d::revertToStart (void)
{ 
	if (theAxialMat != 0)
		theAxialMat->revertToStart();

	K =  Kpen;
	KCommitted = K;
	
	u.Zero();
	up.Zero();
	uCommitted.Zero();
	upCommitted.Zero();
	k = kCommitted = 0.0;
	stress.Zero();
	stressCommitted.Zero();

  return 0;
}

NDMaterial*
FrictionCohesion3d::getCopy (void)
{
	FrictionCohesion3d *theCopy = 0;

	if (!imposedSigma)
		theCopy = new FrictionCohesion3d(this->getTag(), theAxialMat, G, mu, c, GfII);
	else
		theCopy = new FrictionCohesion3d(this->getTag(), G, mu, c, GfII, sigma0);

  return theCopy;
}

NDMaterial*
FrictionCohesion3d::getCopy (const char *type) {
  return this->getCopy();
}

const char*
FrictionCohesion3d::getType (void) const
{
  return "Interface3D";
}

int
FrictionCohesion3d::getOrder (void) const
{
  return 3;
}

//print out material data
void 
FrictionCohesion3d::Print( OPS_Stream &s, int flag )
{
  s << endln ;
  s << "FrictionCohesion3d : " ; 
  s << this->getType( ) << endln ;
  s << endln ;
}

int 
FrictionCohesion3d::sendSelf(int commitTag, Channel &theChannel) 
{  return -1;  }

int 
FrictionCohesion3d::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{  return -1;  }


// methods to overwrite the standard implementation in NDMaterial for asking outputs
//----------------------------------------------------------------------------------

Response*
FrictionCohesion3d::setResponse (const char **argv, int argc, OPS_Stream &output) {

  Response *theResponse = 0;
  output.tag("NdMaterialOutput");
  output.attr("matType",this->getClassType());
  output.attr("matTag",this->getTag());

  if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0) {
    const Vector &res = this->getStress();
    int size = res.Size();
    output.tag("ResponseType","sigma");
	output.tag("ResponseType","tauY"); 
	output.tag("ResponseType", "tauZ");
    theResponse =  new MaterialResponse(this, 1, this->getStress());

  } else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0) {
    const Vector &res = this->getStrain();
    int size = res.Size();
	output.tag("ResponseType","epsilon");
	output.tag("ResponseType","gammaY");
	output.tag("ResponseType", "gammaZ");
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
  }


  output.endTag(); // NdMaterialOutput

  return theResponse;
}

int 
FrictionCohesion3d::getResponse (int responseID, Information &matInfo)
{

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