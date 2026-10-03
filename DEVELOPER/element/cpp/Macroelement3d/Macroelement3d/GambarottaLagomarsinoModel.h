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
Source: /DEVELOPER/element/cpp/Macroelement3d/Macroelement3d/GambarottaLagomarsinoModel.h
Written by Francesco Vanin (francesco.vanin@epfl.ch)
Ecole Polytechnique Federale de Lausanne, Switzerland,
Earthquake Engineering and Structural Dynamics laboratory, 2019

Reference: Vanin F., Penna A., Beyer K.;"A three dimensional macro-element
for modelling of the in-plane and out-of-plane response of masonry walls",
submitted to Earthquake Engineering and Structural Dynamics (2019)

Last edit: 27 Feb 2019
*/
                                                                        
#ifndef GambarottaLagomarsinoModel_h
#define GambarottaLagomarsinoModel_h

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>

class Information;
class Response;
class Channel;
class NDMaterial;

class GambarottaLagomarsinoModel : public NDMaterial
{
  public:
	// full constructor
    GambarottaLagomarsinoModel(int tag, double _E, double _G, double _c, double _mu, double _ct, double _beta, double _L, double _t, double _h, bool _elasticSolution=false);

	//null constructor
    GambarottaLagomarsinoModel();    

	//destructor
    ~GambarottaLagomarsinoModel();

	//methods
    int setTrialStrain (const Vector &strain);
    int setTrialStrain (const Vector &strain, const Vector &rate);
    int setTrialStrainIncr (const Vector &strain);
    int setTrialStrainIncr (const Vector &strain, const Vector &rate);

	const Matrix &getTangent( ) ;
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
	
	double getMode();

	void Print( OPS_Stream &s, int flag );
    int  sendSelf(int commitTag, Channel &theChannel);
	int  recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

	Response* setResponse (const char **argv, int argc, OPS_Stream &output);
	int getResponse (int responseID, Information &matInfo);

	// internal methods
	double toughnessFunction(double alpha);
	double derToughnessFunction(double alpha);
	double sign(double V);

  
//--------------------------------------------------------------------------------
  private:	

    Matrix Kpen;                           // initial stiffness matrix

	double L,t,h;                          // dimension of the rectangular section and height of the panel

	double mu;                             // friction coefficient
	double c;                              // cohesion
	double ct;                             // parameter defining the pre-peak response
	double beta;                           // parameter defining the post-peak response

	double alpha, alphaCommitted;          // damage variable (trial and committed)
	double s, sCommitted;                  // section inelastic slip (trial and committed)
	double deltaLambda;                    // increase of the plastic multiplier
	double deltaAlpha;                     // increase of the damage variable alpha


	Vector stress, stressCommitted;        // section forces (trial, committed)
	Vector u, uCommitted;                  // section displacements (trial, committed)
	Matrix K, KCommitted;                  // stiffness matrix (trial, committed)

	bool elasticSolution;                  // option: linear elastic response (default: false)
	
};

#endif



