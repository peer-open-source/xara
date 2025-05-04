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
// $Date: 2012-05-23 21:11:45 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/PlateRebarMaterialThermal.h,v $

// Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)
//
// Generic Plate Rebar Material
//
/* Ref: Lu X, Lu XZ, Guan H, Ye LP, Collapse simulation of reinforced 
concrete high-rise building induced by extreme earthquakes, 
Earthquake Engineering & Structural Dynamics, 2013, 42(5): 705-723*/

#ifndef PlateRebarMaterialThermal_h
#define PlateRebarMaterialThermal_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <UniaxialMaterial.h>
#include <NDMaterial.h>


class PlateRebarMaterialThermal : public NDMaterial {
public:
  PlateRebarMaterialThermal();
  PlateRebarMaterialThermal(int tag, UniaxialMaterial& uniMat, double ang);

  virtual ~PlateRebarMaterialThermal();

  NDMaterial* getCopy();
  NDMaterial* getCopy(const char* type);

  int getOrder() const;

  const char* getType() const;

  //swap history variables
  int commitState();

  int revertToLastCommit();

  int revertToStart();

  //get the strain
  int setTrialStrain(const Vector& strainFromElement);

  const Vector& getStrain();

  const Vector& getStress();

  const Matrix& getTangent();

  Response* setResponse(const char** argv, int argc, OPS_Stream& output);
  int getResponse(int responseID, Information& matInformation);

  const Matrix& getInitialTangent(); // AV Not Sure if it works
  double getThermalTangentAndElongation(double& TempT, double& ET, double& Elong);

  //density
  double getRho();

  void Print(OPS_Stream& s, int flag);

  int sendSelf(int commitTag, Channel& theChannel);
  int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

private:
  UniaxialMaterial* theMat;
  double angle, c, s;
  double temperature;

  Vector strain;
  static Vector stress;
  static Matrix tangent;
};


#endif
