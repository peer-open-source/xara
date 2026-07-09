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
Source: /DEVELOPER/element/cpp/Macroelement3d/Macroelement3d/CohesiveSurface.h
Written by Francesco Vanin (francesco.vanin@epfl.ch)
Ecole Polytechnique Federale de Lausanne, Switzerland,
Earthquake Engineering and Structural Dynamics laboratory, 2019

Reference: Vanin F., Penna A., Beyer K.;"A three dimensional macro-element
for modelling of the in-plane and out-of-plane response of masonry walls",
submitted to Earthquake Engineering and Structural Dynamics (2019)

Last edit: 27 Feb 2019
*/


#pragma once

#include "Vector.h"
#include "Matrix.h"
#include <string>


class CohesiveSurface
{
public:
	CohesiveSurface(void);
	CohesiveSurface(Matrix& Kpen, double c, double mu0, double muR, double beta, double dropDrift, bool elasticSolution);
	~CohesiveSurface(void);

public:
	double evolveDamage(const Vector& s);              // call damage evolution law
	void updatePlasticStrain(Vector& sigmaTrial);      // call return mapping for platicity

	Matrix getAlgorithmicTangent();
	Matrix getElasticTangent();
	Matrix getCommittedTangent();
	Vector getSigma();

	int setTrialDisplacement(const Vector& s);

	double signTau(double tau);                         // sign function
	double heaviside(double arg);                       // heaviside fuction
	double macauley(double arg);                        // macauley function

	// methods for possible implementation of a hardening/softening plasticity model for friction
	double getMu(double kk);                            
	double getdMu_dk(double kk);

	void commit();
	CohesiveSurface* getCopy();
	double getMode(Vector& epsTrial);
	

protected:


	double mu0;                  // friction coeffcient defining the maximum force capacity 
	double muR;                  // residual friction coeffcient
	double c;                    // cohesion
	double Gc;                   // parameter for pre-peak softening reponse
	double dropDrift;            // parameter defining the post-peak linear force drop

	
	Vector s_di, s_di_nplus1;    // inelastic displacement (committed, trial)
	double D, Dnplus1;           // damage (committed, trial)
	double k;                    // hardening variable (committed)
	double x, xCommitted;        // damage variable x (trial, committed) 

	Vector sigma;                // current stress
	Matrix Kpen;                 // initial stiffness matrix
	Matrix Kalg;                 // algorithmic tangent matrix
	Matrix KCommitted;           // committed stiffness matrix
	Vector dD_ds;                // derivative of D with respect to the displacements s
	double deltaLambda;          // plastic multiplier increase
	bool flagAlg;                // option: use algorithmic tangent or initial
	bool elasticSolution;        // option: define a linear elastic response
};

