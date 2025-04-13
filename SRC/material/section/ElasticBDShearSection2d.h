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
//
// NOTE[cmp] : The only difference between this and ElasticShearSection2d
// is that it is constructed from dimensions b and d???
//
#ifndef ElasticBDShearSection2d_h
#define ElasticBDShearSection2d_h

#include <SectionForceDeformation.h>
#include <Matrix.h>
#include <Vector.h>

class Channel;
class FEM_ObjectBroker;
class Information;

class ElasticBDShearSection2d: public SectionForceDeformation
{
 public:
  ElasticBDShearSection2d(int tag, double E, double b, double d,
			double G, double alpha);
  ElasticBDShearSection2d();    
  ~ElasticBDShearSection2d();
  
  int commitState();
  int revertToLastCommit();
  int revertToStart();
  
  const char *getClassType() const {return "ElasticBDShearSection2d";};
  
  int setTrialSectionDeformation(const Vector&);
  const Vector &getSectionDeformation();
  
  const Vector &getStressResultant();
  const Matrix &getSectionTangent();
  const Matrix &getInitialTangent();
  const Matrix &getSectionFlexibility();
  const Matrix &getInitialFlexibility();
  
  SectionForceDeformation *getCopy();
  const ID &getType();
  int getOrder() const;
  
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel,
	       FEM_ObjectBroker &theBroker);
  
  void Print(OPS_Stream &s, int flag =0);
  
  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int parameterID, Information &info);
  int activateParameter(int parameterID);
  const Vector& getStressResultantSensitivity(int gradIndex,
					      bool conditional);
  const Matrix& getInitialTangentSensitivity(int gradIndex);
  
 protected:
  
 private:
  
  double E, b, d, G, alpha;
  
  Vector e;			// section trial deformations
  
  static Vector s;
  static Matrix ks;
  static ID code;
  
  int parameterID;
};

#endif
