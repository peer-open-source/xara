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

// $Revision: 1.7 $
// $Date: 2006-08-03 23:49:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/MembranePlateFiberSection.h,v $

// Ed "C++" Love
//
// Generic Plate Section with membrane
//

#ifndef MembranePlateFiberSection_h
#define MembranePlateFiberSection_h


#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <NDMaterial.h>

#include <SectionForceDeformation.h>


class MembranePlateFiberSection : public SectionForceDeformation {

public:
  MembranePlateFiberSection();
  MembranePlateFiberSection(int tag, double thickness, NDMaterial&);
  // destructor
  virtual ~MembranePlateFiberSection();


  const char*
  getClassType() const
  {
    return "MembranePlateFiberSection";
  }


  // make a clone of this material
  SectionForceDeformation* getCopy();

  // mass per unit area
  double getRho();

  int getOrder() const;

  const ID& getType();

  // history variables
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

  int sendSelf(int commitTag, Channel&);
  int recvSelf(int commitTag, Channel&, FEM_ObjectBroker&);

  Response* setResponse(const char** argv, int argc, OPS_Stream& s);
  int getResponse(int responseID, Information& info);

private:
  constexpr static int numFibers = 5;

  //quadrature data
  static const double sg[numFibers];
  static const double wg[numFibers];

  double h; //plate thickness

  NDMaterial* theFibers[5]; //pointers to five materials (fibers)

  static const double root56; // =sqrt(5/6)

  Vector strainResultant;

  static Vector stressResultant;

  static Matrix tangent;

  static ID array;

};

#endif
