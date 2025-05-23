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
// $Date: 2012-05-21 23:49:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/LayeredShellFiberSectionThermal.h,v $

// Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)
//
// Layered Shell Section
//
/* Ref: Lu X, Lu XZ, Guan H, Ye LP, Collapse simulation of reinforced 
concrete high-rise building induced by extreme earthquakes, 
Earthquake Engineering & Structural Dynamics, 2013, 42(5): 705-723*/
// Modified for SIF modelling by Liming Jiang [http://openseesforfire.github.io]


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <NDMaterial.h>

#include <SectionForceDeformation.h>


class LayeredShellFiberSectionThermal : public SectionForceDeformation {

  //-------------------Declarations-------------------------------

public:
  //null constructor
  LayeredShellFiberSectionThermal();

  //full constructor
  LayeredShellFiberSectionThermal(int tag, int iLayers, double* thickness, NDMaterial** fibers);

  const char*
  getClassType(void) const
  {
    return "LayeredShellFiberSectionThermal";
  }


  virtual ~LayeredShellFiberSectionThermal();

  SectionForceDeformation* getCopy();

  //mass per unit area
  double getRho();

    int getOrder() const;

  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information& info);

  const ID& getType();

  int commitState();

  int revertToLastCommit();

  int revertToStart();

  //get the strain and integrate plasticity equations
  int setTrialSectionDeformation(const Vector& strain_from_element);

  const Vector& getSectionDeformation();

  const Vector& getStressResultant();

  const Matrix& getSectionTangent();

  //send back the initial tangent
  const Matrix&
  getInitialTangent()
  {
    return this->getSectionTangent();
  }

  void Print(OPS_Stream& s, int flag);

  int sendSelf(int commitTag, Channel& theChannel);
  int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

  const Vector& getTemperatureStress(const Vector&);                          //Added by LMJ
  double determineFiberTemperature(const Vector& DataMixed, double fiberLoc); //Added by LMJ


private:
  int nLayers;
  //quadrature data
  double* sg;
  double* wg;

  double h; //plate thickness

  NDMaterial** theFibers; //pointers to the materials (fibers)

  Vector strainResultant;

  static Vector stressResultant;

  static Matrix tangent;
  static const double root56; // =sqrt(5/6)

  static ID array;

  Vector* sT;                //  Pointer to sTData
  double* ThermalElongation; // thermal elongation
  int countnGauss;
  double AverageThermalForceP;
  double AverageThermalMomentP;

}; //end of LayeredShellFiberSectionThermal declarations
