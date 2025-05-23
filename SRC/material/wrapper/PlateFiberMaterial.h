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

// $Revision: 1.6 $
// $Date: 2006-08-04 18:18:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/PlateFiberMaterial.h,v $

// Ed "C++" Love
//
// Generic Plate Fiber Material
//
#ifndef PlateFiberMaterial_h
#define PlateFiberMaterial_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <NDMaterial.h>

class PlateFiberMaterial : public NDMaterial {

public:
  PlateFiberMaterial();

  PlateFiberMaterial(int tag, NDMaterial& the3DMaterial);

  virtual ~PlateFiberMaterial();

  virtual const char*
  getClassType(void) const
  {
    return "PlateFiberMaterial";
  }

  // make a clone of this material
  NDMaterial* getCopy();
  NDMaterial* getCopy(const char* type);

  int getOrder() const;

  const char* getType() const;

  int commitState();
  int revertToLastCommit();
  int revertToStart();
  int setTrialStrain(const Vector& strainFromElement);

  const Vector& getStrain();
  const Vector& getStress();
  const Matrix& getTangent();
  const Matrix& getInitialTangent();

  //density
  double getRho();

  void Print(OPS_Stream& s, int flag);

  int sendSelf(int commitTag, Channel& theChannel);
  int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

  int setParameter(const char** argv, int argc, Parameter& param);
  Response* setResponse(const char** argv, int argc, OPS_Stream& s);

  const Vector& getStressSensitivity(int gradIndex, bool conditional);

private:
  //out of plane strain
  double Tstrain22;
  double Cstrain22;

  NDMaterial* theMaterial; //pointer to three dimensional material

  Vector strain;

  static Vector stress;

  static Matrix tangent;
}; //end of PlateFiberMaterial declarations


#endif
