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
Source: /DEVELOPER/element/cpp/Macroelement3d/Macroelement3d/DamageShearInterface.h
Written by Francesco Vanin (francesco.vanin@epfl.ch)
Ecole Polytechnique Federale de Lausanne, Switzerland,
Earthquake Engineering and Structural Dynamics laboratory, 2019

Reference: Vanin F., Penna A., Beyer K.;"A three dimensional macro-element
for modelling of the in-plane and out-of-plane response of masonry walls",
submitted to Earthquake Engineering and Structural Dynamics (2019)

Last edit: 27 Feb 2019
*/

                                                                        
#ifndef DamageShearInterface_h
#define DamageShearInterface_h

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>
#include "CohesiveSurface.h"

class Information;
class Response;
class Channel;
class NDMaterial;

class DamageShearInterface : public NDMaterial
{
  public:
	// full constructor
    DamageShearInterface(int tag, double _E, double _G, double c, double mu, double muR, double beta, double dropDrift, bool elasticSolution);

	// constructor for copies
    DamageShearInterface(int tag, CohesiveSurface* _copyCohesiveSurface);

	//null constructor
    DamageShearInterface();    

	//destructor
    ~DamageShearInterface();

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

  
//--------------------------------------------------------------------------------
  private:	

	  CohesiveSurface* cohesiveSurface;         // cohesive surface. Can be updated using a set of cohesive surfaces with different orientations
	  Vector CommittedStrain;                   // committed strain, to avoid useless iterations
	  Vector epsTrial;                          // trial strain

	  Vector sigma;                             // trial stress
	  Matrix K;                                 // trial stiffness matrix
	  bool elasticSolution;                     // option: provide a linar elastic solution (default: false)	
};

#endif



