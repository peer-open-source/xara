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
Source: /DEVELOPER/element/cpp/Macroelement3d/Macroelement3d/GambarottaLagomarsinoModel.cpp
Written by Francesco Vanin (francesco.vanin@epfl.ch)
Ecole Polytechnique Federale de Lausanne, Switzerland,
Earthquake Engineering and Structural Dynamics laboratory, 2019

Reference: Vanin F., Penna A., Beyer K.;"A three dimensional macro-element
for modelling of the in-plane and out-of-plane response of masonry walls",
submitted to Earthquake Engineering and Structural Dynamics (2019)

Last edit: 27 Feb 2019
*/

#include <elementAPI.h>
#include "GambarottaLagomarsinoModel.h"


#include <Channel.h>
#include <math.h>
#include <float.h>
#include <Information.h>
#include <Channel.h>
#include <Message.h>
#include <NDMaterial.h>
#include <MaterialResponse.h>



#ifdef _USRDLL
#define OPS_Export extern "C" _declspec(dllexport)
#elif _MACOSX
#define OPS_Export extern "C" __attribute__((visibility("default")))
#else
#define OPS_Export extern "C"
#endif

static int numGambarottaLagomarsinoModel = 0;
static double tol = 1e-7;
static int maxIter = 100;

OPS_Export void *
OPS_GambarottaLagomarsinoModel()
{
  // print out some KUDO's
  if (numGambarottaLagomarsinoModel == 0) {
    opserr << "GambarottaLagomarsinoModel loaded from external library. Written by Francesco Vanin - EPFL - 2019.\n";
    numGambarottaLagomarsinoModel =1;
  }

  // Pointer to a uniaxial material that will be returned
  NDMaterial *theMaterial = 0;

  //
  // parse the input line for the material parameters
  //
  int    iData[1];
  double dData[9];
  int numData;
  numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag" << endln;
    return 0;
  }

  numData = 9;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid input parameters" << endln;
    return 0;	
  }

  bool _elastic= false;
  while (OPS_GetNumRemainingInputArgs() > 0) {
	  const char* type = OPS_GetString();
	  if (strcmp(type, "-elastic") == 0 || strcmp(type, "-Elastic") == 0) {
		  _elastic = true;
	  }
  }
  
  double h = dData[0];
  double b = dData[1];
  double t = dData[2];
  double E_ = dData[3];
  double G = dData[4];
  double c = dData[5];
  double mu = dData[6];
  double Gc = 1.0 + dData[7];
  double beta = dData[8];

  theMaterial = new GambarottaLagomarsinoModel(0, E_, G, c, mu, Gc / G, beta, b, t, h, _elastic);

  if (theMaterial == 0) {
    opserr << "WARNING could not create NDMaterial of type GambarottaLagomarsinoModel\n";
    return 0;
  }

  // return the material
  return theMaterial;
}

// full constructor
GambarottaLagomarsinoModel::GambarottaLagomarsinoModel(int tag, double _E, double _G, double _c, double _mu, double _ct, double _beta,  double _L, double _t, double _h, bool _elasticSolution)
						  :NDMaterial(tag, 0), stress(2), stressCommitted(2), u(2), uCommitted(2), Kpen(2,2), K(2,2), KCommitted(2,2), elasticSolution(_elasticSolution),
						  c(_c), mu(_mu), ct(_ct), beta(_beta), alpha(0.0), alphaCommitted(0.0), deltaLambda(0.0), deltaAlpha(0.0), s(0.0), sCommitted(0.0),
						  L(_L), t(_t), h(_h)
{ 
   Kpen(0,0) = _E*L*t/h;
   Kpen(1,1) = _G*L*t/h;

   K = Kpen;
   KCommitted = Kpen;
}

GambarottaLagomarsinoModel::GambarottaLagomarsinoModel()
:NDMaterial(0, 0), stress(2), stressCommitted(2), u(2), uCommitted(2), Kpen(2,2), K(2,2), KCommitted(2,2), elasticSolution(false),
			c(0.0), mu(0.0), ct(0.0), beta(0.0), alpha(0.0), alphaCommitted(0.0), deltaLambda(0.0), deltaAlpha(0.0), s(0.0), sCommitted(0.0),
			L(0.0), t(0.0), h(0.0)
{ 
}

GambarottaLagomarsinoModel::~GambarottaLagomarsinoModel() {
}


//internal methods for toughness function
double 
GambarottaLagomarsinoModel::toughnessFunction(double _alpha) {
	
	if (_alpha<=1.0) {
		return 0.5*c*c*ct*h*L*t*_alpha;
	} else {
		return 0.5*c*c*ct*h*L*t*pow(_alpha, -beta);
	}
}

double 
GambarottaLagomarsinoModel::derToughnessFunction(double _alpha) {

	if (_alpha<=1.0) {
		return 0.5*c*c*ct*h*L*t;
	} else {
		return 0.5*c*c*ct*h*L*t*(-beta)*pow(_alpha, -beta-1.0);
	}
}

double 
GambarottaLagomarsinoModel::sign(double _V) {

	if (_V>=0.0) {
		return 1.0;
	} else {
		return -1.0;
	}
}

// standard methods
int 
GambarottaLagomarsinoModel::setTrialStrain(const Vector &strain, const Vector &rate) {
	return this->setTrialStrain(strain);
}

int
GambarottaLagomarsinoModel::setTrialStrain(const Vector &strain) { 

	Vector increment = strain - u;
	if (increment.Norm() > DBL_EPSILON) {
		u = strain;
		stress(0) = Kpen(0,0)*u(0);
		stress(1) = Kpen(1,1)*(u(1) - sCommitted);
		Vector trialStress = stress;
		K = Kpen;

		if (elasticSolution)
			return 0;
		
		// initialise values of state variables to last committed
		alpha = alphaCommitted;
		s = sCommitted;
		deltaLambda = 0.0;
		deltaAlpha = 0.0;

		double csi = 0.;
		if (alphaCommitted>0.0) csi = sCommitted / alphaCommitted;

		double N = stress(0);
		double Vf;
		// kill it here if in tension
		if (N>-DBL_EPSILON) {
			stress(1) = 0.0;
			K.Zero();
			return 0;
		}

		Vector residuals(2);
		// damage function	
		residuals(0) = 0.5*L*t / ( ct*h)*pow(csi,2) - this->toughnessFunction(alpha);
		Vf = Kpen(1,1)*(u(1) - csi*alpha)  - L*t/(ct*h)*csi ;
		residuals(1) = this->sign(Vf)* (Vf) + mu*N; 
		
		if (residuals(0)<=tol && residuals(1)<=tol ) {
			// no evolution of damage parameters, elastic step 
			
		} else {
			// evolution of damage parameters
			if (residuals(0)<=-tol && alpha>DBL_EPSILON) {   // -tol because otherwise we enter here even if there should be evolution of damage: for the first iteration the damage function is equal to the last converged value
				
				// try first sliding without damage evolution (unloading/reloading)			
				deltaLambda = residuals(1) / (Kpen(1,1) + L*t/(ct*h*alpha));
				s += deltaLambda * this->sign(Vf);
								
				// stress update
				stress(1) = Kpen(1,1)*(u(1) - s);

				// algorithmic operator
				double df_dLambda = -L*t/(ct*h*alpha);

				Vector df_dsigma(2);
				df_dsigma(0) = mu;
				df_dsigma(1) = this->sign(Vf);

				Vector dg_dsigma(2);
				dg_dsigma(1) = this->sign(Vf);

				K.Zero();
				K.addMatrixTripleProduct(0.0, Kpen, (dg_dsigma % df_dsigma), Kpen, -1.0);
				K /= (df_dsigma^(Kpen*dg_dsigma)) - df_dLambda;
				K += Kpen;
			} 
			else {
				// sliding with evolution of damage (new loading)
				// Newton Raphson scheme to solve the nonlinear system			
				int iter = 0;
			
				Matrix inverseDerivative(2,2);
				Matrix M(2,2);
				Matrix preM(2,2);
				preM(0,0) = 1.;
				preM(1,1) = this->sign(Vf);

				Vector updates(2);
				double G = Kpen(1,1) / (L*t/h);

				while ( (residuals.Norm() >= tol || (alpha<DBL_EPSILON && iter<=1)) && iter<=maxIter) {
					iter++;
					
					M(0,0) = this->derToughnessFunction(alpha);
					M(0,1) = -L*t/(ct*h)*csi;
					M(1,0) = G*L*t/h*csi * this->sign(Vf);
					M(1,1) = this->sign(Vf)*(Kpen(1,1)*alpha + L*t/(ct*h));

					inverseDerivative(0,0) = M(1,1);
					inverseDerivative(0,1) = -M(0,1);
					inverseDerivative(1,0) = -M(1,0);
					inverseDerivative(1,1) = M(0,0);

					inverseDerivative /= M(0,0)*M(1,1) - M(0,1)*M(1,0);
					updates = (inverseDerivative)*residuals;
				
					deltaAlpha  += updates(0);
					if (deltaAlpha<DBL_EPSILON)      deltaAlpha=0.0;
					alpha += updates(0);
					csi += updates(1);

					deltaLambda = (csi*alpha - sCommitted) / this->sign(Vf);
					s = csi*alpha;

					residuals(0) = 0.5*L*t / ( ct*h)*pow(csi,2) - this->toughnessFunction(alpha);
					residuals(1) = this->sign(Vf)* ( Kpen(1,1)*(u(1) - csi*alpha)  - L*t/(ct*h)*csi ) + mu*N;

					
					if (alpha<alphaCommitted)   {  // if we enter here: there should have been the update only for deltaLambda and we are here because the damage function is close to zero

						alpha=alphaCommitted;
						deltaAlpha = 0.0;
						csi = sCommitted / alphaCommitted;					
						s = sCommitted;

						residuals(1) = this->sign(Vf)* ( Kpen(1,1)*(u(1) - csi*alpha)  - L*t/(ct*h)*csi ) + mu*N;

						csi +=  residuals(1) / ( this->sign(Vf)*(Kpen(1,1)*alpha + L*t/(ct*h)) );
						deltaLambda = (csi*alpha - sCommitted) / this->sign(Vf);
						s = csi*alpha;

						residuals(1) = this->sign(Vf)* ( Kpen(1,1)*(u(1) - csi*alpha)  - L*t/(ct*h)*csi ) + mu*N;

						residuals(0) = 0.5*L*t / ( ct*h*pow(alpha,2) )*pow(s,2) - this->toughnessFunction(alpha);
						if (residuals(0)<0.0)
							residuals(0) = 0.0;
						
					}
					
				}  // end while (NR iterations)

				// stress update
				stress(1) = Kpen(1,1)*(u(1) - s);

				if (iter>maxIter) {
					// maybe insert here warning or error
					K=Kpen*0;
					stress(1) = -mu*N* this->sign(Vf);
					s = sCommitted;
					alpha = alphaCommitted;
					deltaLambda = 0.0;
					deltaAlpha = 0.0;

					return 0;

				}	else {

				// tangent operator update
				double hh;				
				hh = -(L*t*alpha*s*this->sign(Vf)) / (L*t*s*s + ct*h*pow(alpha,3)* this->derToughnessFunction(alpha));

				if (deltaAlpha<= DBL_EPSILON)   hh = 0.0;

				double df_dalpha = (L*t*s*this->sign(Vf)) / (ct*h*pow(alpha,2));
				double df_dLambda = -L*t/(ct*h*alpha);

				Vector df_dsigma(2);
				df_dsigma(0) = mu;
				df_dsigma(1) = this->sign(Vf);

				Vector dg_dsigma(2);
				dg_dsigma(1) = this->sign(Vf);

				K.Zero();
				K.addMatrixTripleProduct(0.0, Kpen, (dg_dsigma % df_dsigma), Kpen, -1.0);
				K /= df_dalpha*hh + (df_dsigma^(Kpen*dg_dsigma)) - df_dLambda;
				K += Kpen;
				}
			}  // end second case
		}      // end evolution of damage parameters (2 cases)
	}          // end set new strain
    
	return 0;
}

int
GambarottaLagomarsinoModel::setTrialStrainIncr (const Vector &strain)
{
	opserr << "WARNING: GambarottaLagomarsinoModel has no method setTrialStrainIncr implemented" << endln; 
	return 0;
}

int
GambarottaLagomarsinoModel::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  return this->setTrialStrainIncr(strain);
}

const Matrix& 
GambarottaLagomarsinoModel::getCommittedTangent( )  {
	return KCommitted;
}
	
const Matrix& 
GambarottaLagomarsinoModel::getTangent( ) 
{
	if (alpha>1.0 && alphaCommitted<1.0) {
		return KCommitted;
	} else {
		if (abs(K(0,1))>-DBL_EPSILON && abs(K(1,1))>-DBL_EPSILON)		
			return K;
		else
			return Kpen;
	}
}

const Matrix&
GambarottaLagomarsinoModel::getInitialTangent (void)
{ 
	return Kpen;
}

const Vector&
GambarottaLagomarsinoModel::getStress (void)
{
	return stress;
}

const Vector&
GambarottaLagomarsinoModel::getStrain (void)
{
  return u;
}

int
GambarottaLagomarsinoModel::commitState (void)
{
	alphaCommitted = alpha;
	sCommitted = s;
	KCommitted = K;

	uCommitted = u;
	stressCommitted = stress;

	return 0;
}

int
GambarottaLagomarsinoModel::revertToLastCommit (void)     
{ 
	alpha = alphaCommitted;
	s = sCommitted;
	K = KCommitted;

	u = uCommitted;
	stress = stressCommitted;

  return 0;
}

int
GambarottaLagomarsinoModel::revertToStart (void)
{ 
	alpha = 0;
	s = 0;
	K = Kpen;

	alphaCommitted = 0;
	sCommitted = 0;
	KCommitted = Kpen;

	u.Zero();
	stress.Zero();

  return 0;
}

double 
GambarottaLagomarsinoModel::getMode(void) {
	return alphaCommitted;
}

NDMaterial*
GambarottaLagomarsinoModel::getCopy (void)
{
	GambarottaLagomarsinoModel *theCopy = new GambarottaLagomarsinoModel(this->getTag(), Kpen(0,0)/(L*t/h), Kpen(1,1)/(L*t/h), c, mu, ct, beta, L, t, h, elasticSolution);

  return theCopy;
}

NDMaterial*
GambarottaLagomarsinoModel::getCopy (const char *type) {
  return this->getCopy();
}

const char*
GambarottaLagomarsinoModel::getType (void) const
{
  return "BeamFiber";
}

int
GambarottaLagomarsinoModel::getOrder (void) const
{
  return 2;
}

//print out material data
void 
GambarottaLagomarsinoModel::Print( OPS_Stream &s, int flag )
{
  s << endln ;
  s << "GambarottaLagomarsinoModel : " ; 
  s << this->getType( ) << endln ;
  s << endln ;
}

int 
GambarottaLagomarsinoModel::sendSelf(int commitTag, Channel &theChannel) 
{  return -1;  }

int 
GambarottaLagomarsinoModel::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{  return -1;  }


// methods to overwrite the standard implementation in NDMaterial for asking outputs
//----------------------------------------------------------------------------------

Response*
GambarottaLagomarsinoModel::setResponse (const char **argv, int argc, OPS_Stream &output) {

  Response *theResponse = 0;
  output.tag("NdMaterialOutput");
  output.attr("matType",this->getClassType());
  output.attr("matTag",this->getTag());

  if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0) {
    const Vector &res = this->getStress();
    int size = res.Size();
    output.tag("ResponseType","sigma");
	output.tag("ResponseType","tau");   
    theResponse =  new MaterialResponse(this, 1, this->getStress());

  } else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0) {
    const Vector &res = this->getStrain();
    int size = res.Size();
	output.tag("ResponseType","epsilon");
	output.tag("ResponseType","gamma");
    theResponse =  new MaterialResponse(this, 2, this->getStrain());

  } else if (strcmp(argv[0],"tangent") == 0 || strcmp(argv[0],"stiffness") == 0) {
	const Matrix &res = this->getCommittedTangent();
    int size = 4;
	output.tag("ResponseType","K11");
	output.tag("ResponseType","K12");
	output.tag("ResponseType","K21");
	output.tag("ResponseType","K22");
    theResponse =  new MaterialResponse(this, 3, this->getCommittedTangent());
  }

    else if (strcmp(argv[0], "stateVariables") == 0 || strcmp(argv[0], "state") == 0) {
	  double res = this->getMode();
	  int size = 1;
	  output.tag("ResponseType", "mode");
	  theResponse = new MaterialResponse(this, 4, this->getMode());
  }
  

  output.endTag(); // NdMaterialOutput

  return theResponse;
}


int 
GambarottaLagomarsinoModel::getResponse (int responseID, Information &matInfo)
{
  Vector ep_dot(2);
  double mode;

  switch (responseID) {
  case 1:
    return matInfo.setVector(this->getStress());
    
  case 2:
    return matInfo.setVector(this->getStrain());

  case 3:
	  return matInfo.setMatrix(this->getCommittedTangent());

  case 4:
	  return matInfo.setDouble(this->getMode());

  default:
    return -1;
  }
}