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
Source: /DEVELOPER/element/cpp/Macroelement3d/Macroelement3d/CohesionFriction3d.h
Written by Igor Tomić (igor.tomic@epfl.ch)
Ecole Polytechnique Federale de Lausanne, Switzerland,
Earthquake Engineering and Structural Dynamics laboratory, 2019

Last edit: 27 October 2019
*/
                                                                        
#ifndef CohesionFriction3d_h
#define CohesionFriction3d_h

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>


class Information;
class Response;
class Channel;
class NDMaterial;
class UniaxialMaterial;


class CohesionFriction3d : public NDMaterial
{
  public:
	// full constructor
    CohesionFriction3d(int tag, UniaxialMaterial* comprMat, double G, double mu, double c, double GfII);
	CohesionFriction3d(int tag, double G, double mu, double c, double GfII, double sigma0);

	//null constructor
    CohesionFriction3d();    

	//destructor
    ~CohesionFriction3d();

	//methods
    int setTrialStrain (const Vector &strain);
    int setTrialStrain (const Vector &strain, const Vector &rate);
    int setTrialStrainIncr (const Vector &strain);
    int setTrialStrainIncr (const Vector &strain, const Vector &rate);

	const Matrix &getTangent(void) ;
    const Matrix &getInitialTangent (void);
    const Vector &getStress (void);
    const Vector &getStrain (void);
	const Matrix& getCommittedTangent(void);

	int commitState (void);
    int revertToLastCommit (void);
    int revertToStart (void);
    
    NDMaterial *getCopy (void);
	NDMaterial *getCopy (const char *type);
    const char *getType (void) const;
    int getOrder (void) const;

	void Print( OPS_Stream &s, int flag );
    int  sendSelf(int commitTag, Channel &theChannel);
	int  recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

	Response* setResponse (const char **argv, int argc, OPS_Stream &output);
	int getResponse (int responseID, Information &matInfo);

	double hardeningLaw(double k);
	double der_hardeningLaw(double k);
  
//--------------------------------------------------------------------------------
  private:	

    Matrix Kpen;                        // elastic stiffness matrix
	double E;                           // stiffness in axial direction
	double G;                           // stiffness in shear direction
	UniaxialMaterial*  theAxialMat;     // uniaxial material model defining the axial response

	Vector stress, stressCommitted;     // stress vector (trial, committed)
	Vector u, uCommitted;               // deformation (trial, committed)
	Vector up, upCommitted;             // deformation (trial, committed)
	Matrix K, KCommitted;		        // stiffness matrix (trial, committed)

	// parameters
	double c, mu, GfII;                 // cohesion, friction coeffcient, fracture energy mode II
	double k, deltaLambda, kCommitted;  // hardening parameter, plastic increment
	double sigma0;                      // imposed value of axial stress, if any
	bool imposedSigma;                  // impose, or not, one axial load.

};


#endif



