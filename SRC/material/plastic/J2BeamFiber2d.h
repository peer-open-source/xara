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
// Written: MHS
// Created: Aug 2001
//
// Description: Elastic isotropic model where stress components 22, 33, and 23
// are condensed out.
//
#ifndef J2BeamFiber2d_h
#define J2BeamFiber2d_h

#define ND_TAG_J2BeamFiber2d 91625

#include <NDMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Information.h>
#include <Parameter.h>

class J2BeamFiber2d : public NDMaterial {
public:
  J2BeamFiber2d(int tag, double E, double G, double sigY, double Hi, double Hk);
  J2BeamFiber2d();
  ~J2BeamFiber2d();

  int setTrialStrain(const Vector& v);
  int setTrialStrain(const Vector& v, const Vector& r);
  int setTrialStrainIncr(const Vector& v);
  int setTrialStrainIncr(const Vector& v, const Vector& r);
  const Matrix& getTangent(void);
  const Matrix& getInitialTangent(void);
  const Vector& getStress(void);
  const Vector& getStrain(void);

  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);

  NDMaterial* getCopy(void);
  NDMaterial* getCopy(const char* type);
  const char* getType(void) const;
  int getOrder(void) const;

  int sendSelf(int commitTag, Channel& theChannel);
  int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

  void Print(OPS_Stream& s, int flag = 0);

  int setParameter(const char** argv, int argc, Parameter& param);
  int updateParameter(int parameterID, Information& info);
  int activateParameter(int paramID);

  const Vector& getStressSensitivity(int gradIndex, bool conditional);
  int commitSensitivity(const Vector& depsdh, int gradIndex, int numGrads);

protected:
private:
  double E;
  double nu;
  double sigmaY;
  double Hiso;
  double Hkin;

  int parameterID;
  Matrix* SHVs;

  static Vector sigma; // Stress vector ... class-wide for returns
  static Matrix D;     // Elastic constants
  Vector Tepsilon;     // Trial strains

  double alphan;
  double alphan1;

  double epsPn[2];
  double epsPn1[2];

  double dg_n1;
};


#endif
