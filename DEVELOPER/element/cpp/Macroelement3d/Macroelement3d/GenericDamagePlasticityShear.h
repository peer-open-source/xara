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

#ifndef GenericDamagePlasticityShear_h
#define GenericDamagePlasticityShear_h

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include <string>


class GenericDamagePlasticityShear : public NDMaterial
{
public:
	GenericDamagePlasticityShear(void);
	GenericDamagePlasticityShear(int tag, double G, int modelType, Vector matProp, double Gc, double dropDrift, double alpha, bool elasticSolution, double initialDamage=0.0);
	~GenericDamagePlasticityShear(void);

public:
	double evolveDamage(const Vector& s);              // call damage evolution law
	void updatePlasticStrain(Vector& sigmaTrial, double maxShearCapacity);      // call return mapping for platicity

	const Matrix &getTangent();
	const Matrix &getInitialTangent(void);
	const Vector &getStress(void);
	const Vector &getStrain(void);
	const Matrix& getCommittedTangent(void);

	int setTrialStrain(const Vector &s);
	int setTrialStrain(const Vector &s, const Vector &rate);
	int setTrialStrainExtraParams(const Vector &s, Vector params);

	double signTau(double tau);                         // sign function
	double heaviside(double arg);                       // heaviside fuction
	double macauley(double arg);                        // macauley function

	double Vmax(const Vector& s);
	double dVmax_dN(const Vector& s);


	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);


	NDMaterial *getCopy(void);
	NDMaterial *getCopy(const char *type);
	const char *getType(void) const;
	int getOrder(void) const;
	double getMode();


	void Print(OPS_Stream &s, int flag);
	int  sendSelf(int commitTag, Channel &theChannel);
	int  recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

	

protected:

	int modelType;
	Vector props;              // vector with the generic material properties list of the mode
	double Gc;                   // parameter for pre-peak softening reponse
	double dropDrift;            // parameter defining the post-peak linear force drop
	double alpha;                // scale factor between the yield surface and the full strength domain
	Vector params;               // additional parameters passed by the macroelement (H0, axial load ratio)

	
	Vector s_di, s_di_nplus1;    // inelastic displacement (committed, trial)
	double D, Dnplus1;           // damage (committed, trial)
	double k;                    // hardening variable (committed)
	double x, xCommitted;        // damage variable x (trial, committed) 

	Vector eps;                  // last trial displacements
	Vector sigma;                // current stress
	Matrix Kpen;                 // initial stiffness matrix
	Matrix Kalg;                 // algorithmic tangent matrix
	Matrix KCommitted;           // committed stiffness matrix
	Vector dD_ds;                // derivative of D with respect to the displacements s
	double deltaLambda;          // plastic multiplier increase
	bool flagAlg;                // option: use algorithmic tangent or initial
	bool elasticSolution;        // option: define a linear elastic response
};

#endif