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
Source: /DEVELOPER/element/cpp/Macroelement3d/Macroelement3d/WrappedMaterial.h
Written by Francesco Vanin (francesco.vanin@epfl.ch)
Ecole Polytechnique Federale de Lausanne, Switzerland,
Earthquake Engineering and Structural Dynamics laboratory, 2019

Builds a wrapped nDMaterial (dim: 2) that uses a linear elastic model for the
first (axial) response and a generic material model for the second (shear) response.
No coupling between the two directions.
Scope: defining the shear response of the macroelement through a uniaxial material
model, and apply it to a nDMaterial to be attached to the macroelement.

Last edit: 27 Feb 2019
*/
                                                                        
#ifndef WrappedMaterial_h
#define WrappedMaterial_h

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>


class Information;
class Response;
class Channel;
class NDMaterial;
class UniaxialMaterial;


class WrappedMaterial : public NDMaterial
{
  public:
	// full constructor
    WrappedMaterial(int tag, double E,  UniaxialMaterial* shearMat);

	//null constructor
    WrappedMaterial();    

	//destructor
    ~WrappedMaterial();

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

	void Print( OPS_Stream &s, int flag );
    int  sendSelf(int commitTag, Channel &theChannel);
	int  recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

	Response* setResponse (const char **argv, int argc, OPS_Stream &output);
	int getResponse (int responseID, Information &matInfo);
  
//--------------------------------------------------------------------------------
  private:	

    Matrix Kpen;                        // elastic stiffness matrix
	double E;                           // stiffness in axial direction
	UniaxialMaterial*  theShearMat;     // uniaxial material model defining the shear response

	Vector stress, stressCommitted;     // stress vector (trial, committed)
	Vector u, uCommitted;               // deformation (trial, committed)
	Matrix K, KCommitted;		        // stiffness matrix (trial, committed)

};


#endif



