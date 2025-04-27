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

// $Revision: 1.0 $
// $Date: 2012-05-27 21:11:45 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/PlateFromPlaneStressMaterial.h,v $

// Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)
//
// Generic Plate Material from Plane Stress Material
//
/* Ref: Lu X, Lu XZ, Guan H, Ye LP, Collapse simulation of reinforced 
concrete high-rise building induced by extreme earthquakes, 
Earthquake Engineering & Structural Dynamics, 2013, 42(5): 705-723*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <NDMaterial.h>


class PlateFromPlaneStressMaterial : public NDMaterial {
public:
  PlateFromPlaneStressMaterial();
  PlateFromPlaneStressMaterial(int tag, NDMaterial&, double G);

  virtual ~PlateFromPlaneStressMaterial();

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
  const Matrix& getInitialTangent(); // AV Not Sure if it works

  //density
  double getRho();

  void Print(OPS_Stream& s, int flag);

  int sendSelf(int commitTag, Channel& theChannel);
  int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

  Response* setResponse(const char** argv, int argc, OPS_Stream& s);

private:
  NDMaterial* theMat; //pointer to three dimensional material
  double gmod;

  Vector strain;
  static Vector stress;
  static Matrix tangent;
};
