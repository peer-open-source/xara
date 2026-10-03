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
                                                                        
// $Revision: 1.1 $
// $Date: 2008/12/09 20:00:16 $
// $Source: /usr/local/cvs/OpenSees/PACKAGES/NewMaterial/cpp/ElasticPPcpp.h,v $
                                                                        
#ifndef ElasticPPcpp_h
#define ElasticPPcpp_h

// Written: fmk 
//
// Description: This file contains the class definition for 
// ElasticPPcpp. ElasticPPcpp provides the abstraction
// of an elastic perfectly plastic uniaxial material, 
//
// What: "@(#) ElasticPPcpp.h, revA"

#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>

class Information;
class Response;
class Channel;
class NDMaterial;

class BeamFrictionSupport : public NDMaterial
{
  public:
    BeamFrictionSupport(int tag, double E, double G, double mu, double anchoringLength, double gapX, double gapY);    
    BeamFrictionSupport();    

    ~BeamFrictionSupport();

	int setTrialStrain(const Vector &v);
    int setTrialStrain(const Vector &v, const Vector &r);
    int setTrialStrainIncr(const Vector &v);
    int setTrialStrainIncr(const Vector &v, const Vector &r);
    
	const Matrix& getTangent(void);
    const Matrix& getInitialTangent(void);
	const Matrix& getCommittedTangent(void);

    const Vector& getStress(void);
    const Vector& getStrain(void);

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    NDMaterial *getCopy(void);
    NDMaterial *getCopy(const char *code);

    const char *getType(void) const;
    int getOrder(void) const;  

    Response *setResponse (const char **argv, int argc, OPS_Stream &s);
    int getResponse (int responseID, Information &matInformation);

	void Print( OPS_Stream &s, int flag );
    int  sendSelf(int commitTag, Channel &theChannel);
	int  recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
  protected:
    
  private:

	  double mu;          // friction coefficient
	  double dx;          // anchorage length (positive value, negative value for infinite sliding)
	  double gapX, gapY;  // gap of the embedded beam (positive value, negative value for infinite sliding)

	  Vector slidingTrial, sliding;    // plastic sliding, trial and committed
	  Vector e, eTrial, s, sTrial;     // strain and stress vectors (committed and trial) 

	  Matrix D, Dtrial, Dcommitted;   // stiffness matrix

	  bool failed, failedTrial; 
	  double fakeStiffnessFactor;
};


#endif



