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
Source: /DEVELOPER/element/cpp/Macroelement3d/Macroelement3d/CohesiveSurface.cpp
Written by Francesco Vanin (francesco.vanin@epfl.ch)
Ecole Polytechnique Federale de Lausanne, Switzerland,
Earthquake Engineering and Structural Dynamics laboratory, 2019

Reference: Vanin F., Penna A., Beyer K.;"A three dimensional macro-element
for modelling of the in-plane and out-of-plane response of masonry walls",
submitted to Earthquake Engineering and Structural Dynamics (2019)

Last edit: 27 Feb 2019
*/

/* -------------------------------------------------------------------------- */
#include "CohesiveSurface.h"
/* -------------------------------------------------------------------------- */
#include <OPS_Globals.h>
#include <OPS_Stream.h>
#include <TaggedObject.h>
#include <MapOfTaggedObjects.h>
/* -------------------------------------------------------------------------- */
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
//#include <string.h>
#include <cmath>
/* -------------------------------------------------------------------------- */


#ifndef DBL_EPSILON
#define DBL_EPSILON (std::numeric_limits<double>::epsilon())
#endif

CohesiveSurface::CohesiveSurface(void)
	:Kalg(2, 2), Kpen(2, 2), sigma(2), s_di(2), s_di_nplus1(2) {
}

CohesiveSurface::CohesiveSurface(Matrix& Kpen, double c, double mu0, double muR, double Gc, double dropDrift, bool elasticSolution)
	: Kalg(Kpen), Kpen(Kpen), KCommitted(Kpen), sigma(2), s_di(2), dD_ds(2), mu0(mu0), muR(muR), c(c), Gc(Gc), dropDrift(dropDrift), flagAlg(false), deltaLambda(0.0), s_di_nplus1(2), 
	  elasticSolution(elasticSolution), k(0.), x(0.), xCommitted(0.) {
		
		D = 0.;
		Dnplus1 = D;
}

CohesiveSurface::~CohesiveSurface(void) {
}

double 
CohesiveSurface::signTau(double tau) {
	if (tau == 0.0)
		return 0.0;
	else if (tau > 0.0)
		return 1.0;
	else
		return -1.0;
}

double
CohesiveSurface::heaviside(double arg) {
	if (arg > 0.0)
		return 1.0;
	else
		return 0.0;
}

double 
CohesiveSurface::macauley(double arg) {
	if (arg > 0.0)
		return arg;
	else
		return 0.0;
}

double 
CohesiveSurface::getMu(double kk) {
      return muR; 
}

double 
CohesiveSurface::getdMu_dk(double kk) {
    return 0; 
}

void 
CohesiveSurface::updatePlasticStrain(Vector& sigmaTrial) {

	// return mapping for constant muR, no iterations needed
	double yieldFunction = muR*sigmaTrial(0) + abs(sigmaTrial(1));
	deltaLambda = 0.0;
	s_di_nplus1 = s_di;

	if (yieldFunction > DBL_EPSILON) {
		deltaLambda = 1.0 / (Kpen(1, 1)) * yieldFunction;
		s_di_nplus1(1) += deltaLambda * this->signTau(sigmaTrial(1));
		double sign = this->signTau(sigmaTrial(1));
		sigmaTrial(1) -= deltaLambda * Kpen(1, 1) * sign;
	}

	return;
}

double 
CohesiveSurface::evolveDamage(const Vector& s) {
	
	double Dtrial;
	double n = min(0.0,Kpen(0,0)*s(0));
	
	// force-capacity model 
	double VMax = c - n*mu0;          // mohr-coulomb law
	double dVMax_dn = -1.0*mu0;       // derivative of VMax with respect to N

	double s0 = -n*muR/Kpen(1,1);
	double sMax = -muR*n/Kpen(1,1) + Gc*(VMax - (-n)*muR)/Kpen(1,1);
	double sDrop = dropDrift;         // constant drift model

	double zeta = 0.20;               // capacity drop percentage
	double slope = zeta *(c-mu0*n)/ (std::max(sDrop, 1.05*sMax) - sMax);

	if (abs(s(1))<s0) {
		x = 0.0;
		Dtrial = 0.0;
	} else {
		if (abs(s(1))<sMax) {
			x = (abs(s(1)) - (-muR*n/Kpen(1,1)))/(sMax - (-muR*n/Kpen(1,1)));
			Dtrial = (1.0 - 1.0/Gc)*pow(x, 1.0/(Gc - 1.0));
		} else {
			x = (abs(s(1)) - (-muR*n/Kpen(1,1)))/(sMax - (-muR*n/Kpen(1,1)));  // not used but calculated for output
			Dtrial = (Kpen(1,1)*Kpen(1,1)*abs(s(1)) + muR*n*slope + Kpen(1,1)*abs(s(1))*slope - muR*n*Gc*slope - Kpen(1,1)*VMax - Gc*slope*VMax) / (Kpen(1,1)*(muR*n + Kpen(1,1)*abs(s(1))));
		}
	}

	if (Dtrial < 0.0)   Dtrial = 0.0;
	if (Dtrial > 1.0)   Dtrial = 1.0;
	
	if (Dtrial == 1.0 || D > Dtrial) {
		flagAlg = false;
	    dD_ds.Zero();

	} else {
		flagAlg = true;
		dD_ds.Zero();

		double sign;
		if (s(1)>=0.0)
			sign = 1.0;
		else
			sign = -1.0;

	    if (abs(s(1))<=sMax && Dtrial>0.0) {

			dD_ds(0)  = -((pow((muR*n + Kpen(1,1)*abs(s(1)))/(Gc*muR*n + Gc*VMax),pow(-1 + Gc,-1))*(k*muR*abs(s(1)) - muR*VMax + (muR*n + Kpen(1,1)*abs(s(1)))* dVMax_dn))
				         /(Gc*(muR*n + Kpen(1,1)*abs(s(1)))*(muR*n + VMax))); 
			
			dD_ds(0)*= Kpen(0,0);

	    	dD_ds(1)  = sign*  (Kpen(1,1)*pow((muR*n + Kpen(1,1)*abs(s(1)))/(Gc*muR*n + Gc*VMax),pow(-1 + Gc,-1)))/(Gc*(muR*n + Kpen(1,1)*abs(s(1))));
	    } else {

			double dSlope_dN = (Kpen(1,1)*zeta*(muR*(-1 + Gc)*VMax + (Kpen(1,1)*sDrop + muR*(n - n*Gc))*dVMax_dn))/pow(-(Kpen(1,1)*sDrop) + muR*n*(-1 + Gc) + Gc*VMax, 2.0);

			if (sDrop<0.0) {

				dD_ds(1)  = -((mu0*n - muR*n - c)*(Kpen(1,1) + 2.0*slope)* sign) / pow(muR*n + Kpen(1,1)*abs(s(1)), 2);

				dD_ds(0)  = ((c*muR + Kpen(1,1)*(mu0-muR)*abs(s(1)))*(Kpen(1,1) + 2*slope) + (muR*n + Kpen(1,1)*abs(s(1)))*(2*mu0*n - muR*n - 
					          2*c + Kpen(1,1)*abs(s(1)))*dSlope_dN)/(Kpen(1,1)*pow(muR*n + Kpen(1,1)*abs(s(1)),2));
			    
				dD_ds(0)*= Kpen(0,0);

			}   else {

				dD_ds(0)  = (-((Kpen(1,1) + Gc*slope)*(Kpen(1,1)*muR*abs(s(1)) - muR*VMax + (muR*n + Kpen(1,1)*abs(s(1)))* dVMax_dn)) + 
                            (muR*n + Kpen(1,1)*abs(s(1)))*(Kpen(1,1)*abs(s(1)) + muR*(n - n*Gc) - Gc*VMax)* dSlope_dN)/(Kpen(1,1)*pow(muR*n + Kpen(1,1)*abs(s(1)),2.0));
			    
				dD_ds(0)*= Kpen(0,0);
				
								
				dD_ds(1)  = sign*((Kpen(1,1) + Gc*slope)*(muR*n + VMax))/pow(muR*n + Kpen(1,1)*abs(s(1)), 2.0);
	    }  
		}		
	}

	if (x<xCommitted)
		x=xCommitted;

	if (Dtrial > D)
		return Dtrial;
	else
		return D;

}

int 
CohesiveSurface::setTrialDisplacement(const Vector& s) {

	if (elasticSolution) {

		sigma = Kpen * s;
		Kalg = Kpen;

	} else {

    x = xCommitted;
	Dnplus1 = this->evolveDamage(s);
	Matrix H = Kpen;
	H(0,0) *= 1.0 - this->heaviside(s(0));
	H(1,1) *= 1.000;
	Vector sigmaTrial = H*(s-s_di);
	
	this->updatePlasticStrain(sigmaTrial);

	sigma = ((1.0 - Dnplus1) * (Kpen*s)) + (Dnplus1 * sigmaTrial);	

	// update tangent matrix
	Matrix Ht(2, 2);
	if (s(0) <= -1.0*DBL_EPSILON) {
		if (deltaLambda > DBL_EPSILON) {
			  Ht(0, 0) = Kpen(0, 0);
			  Ht(1, 0) = -1.0*muR*Kpen(0,0) *this->signTau(sigma(1));
		}
		else {
			Ht = Kpen;
			Ht(1,1) *= 1.000;
		}
	}

	Kalg = (1.0 - Dnplus1)*Kpen + Dnplus1*Ht;
	
	if (flagAlg) {
	  if (Dnplus1>D) {
		Kalg -= Kpen* (s % dD_ds);
		Kalg += sigmaTrial % dD_ds;
	  }

	}
	}

	return 0;

}

double 
CohesiveSurface::getMode(Vector& s) {
	return xCommitted;
}

Vector 
CohesiveSurface::getSigma() {
	return sigma;
}

Matrix 
CohesiveSurface::getElasticTangent() {
	return Kpen;
}

Matrix 
CohesiveSurface::getAlgorithmicTangent() {
	if (elasticSolution) {
		return Kpen;
	} else {
		return Kalg;
	}
}

Matrix
CohesiveSurface::getCommittedTangent() {
	return KCommitted;
}
	
void
CohesiveSurface::commit() {
	if (!elasticSolution) {
		D = Dnplus1;
		xCommitted = x;
		s_di = s_di_nplus1;
		KCommitted = Kalg;
		k += deltaLambda;
	}
	return;
}

CohesiveSurface* 
CohesiveSurface::getCopy() {
	CohesiveSurface* theCopy;
	theCopy = new CohesiveSurface(Kpen, c, mu0, muR, Gc, dropDrift, elasticSolution);	
	return theCopy;
}
