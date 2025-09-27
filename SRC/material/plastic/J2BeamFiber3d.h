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


// $Revision: 1.3 $
// $Date: 2002-12-05 22:49:09 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/J2BeamFiber3d.h,v $

// Written: MHS
// Created: Aug 2001
//
// Description: Elastic isotropic model where stress components 22, 33, and 23
// are condensed out.

#ifndef J2BeamFiber3d_h
#define J2BeamFiber3d_h

#define ND_TAG_J2BeamFiber3d 92516

#include <NDMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Information.h>
#include <Parameter.h>
#include <Matrix3D.h>
#include <VectorND.h>
#include <Vector3D.h>


class J2BeamFiber3d : public NDMaterial {
public:
  J2BeamFiber3d(int tag, double E, double nu, double sigY, double Hi, double Hk);
  J2BeamFiber3d();
  ~J2BeamFiber3d();

  const char* getClassType() const override { return "J2BeamFiber3d"; }

  int setTrialStrain(const Vector&) override;
  int setTrialStrainIncr(const Vector& v) override;
  const Matrix& getTangent() override;
  const Matrix& getInitialTangent() override;
  const Vector& getStress() override;
  const Vector& getStrain() override;
  Vector3D getFrameStress();
  OpenSees::Matrix3D getFrameTangent();

  int commitState() override;
  int revertToLastCommit() override;
  int revertToStart() override;

  NDMaterial* getCopy() override;
  NDMaterial* getCopy(const char* type) override;
  const char* getType() const override;
  int getOrder() const;

  int sendSelf(int commitTag, Channel& ) override;
  int recvSelf(int commitTag, Channel& , FEM_ObjectBroker&) override;

  void Print(OPS_Stream& s, int flag) override;

  int setParameter(const char** argv, int argc, Parameter&) override;
  int updateParameter(int parameterID, Information&) override;
  int activateParameter(int paramID);

  const Vector& getStressSensitivity(int gradIndex, bool conditional);
  int commitSensitivity(const Vector& depsdh, int gradIndex, int numGrads);


private:
  double E;
  double nu;
  double sigmaY;
  double Hiso;
  double Hkin;

  int parameterID;
  Matrix* SHVs;

  static Matrix D;     // Elastic constants
  Vector Tepsilon;     // Trial strains

  double alphan;
  double alphan1;

  double epsPn[3];
  double epsPn1[3];
  // Vector3D xsi;

  double dg_n1;
};


#endif
