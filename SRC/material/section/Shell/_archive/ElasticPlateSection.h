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
//  Elastic Plate Section
//
// Ed "C++" Love
//
#include <math.h>

#include <Vector.h>
#include <Matrix.h>
#include <ID.h>

#include <SectionForceDeformation.h>


class ElasticPlateSection : public SectionForceDeformation {

public:
  ElasticPlateSection();
  ElasticPlateSection(int tag, double E, double nu, double h = 1.0);
  ~ElasticPlateSection();

  SectionForceDeformation* getCopy();


  int getOrder() const;

  const ID& getType();

  int commitState();
  int revertToLastCommit();
  int revertToStart();
  int setTrialSectionDeformation(const Vector& strain_from_element);

  const Vector& getSectionDeformation();
  const Vector& getStressResultant();
  const Matrix& getSectionTangent();
  const Matrix& getInitialTangent();

  void Print(OPS_Stream& s, int flag);

  int sendSelf(int commitTag, Channel& theChannel);
  int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

private:
  double E;  // elastic modulus
  double nu; // poisson ratio
  double h;  // plate thickness

  static constexpr double shear_factor = 5.0/6.0; // shear correction factor

  Vector strain;
  static Vector stress;
  static Matrix tangent;
  static ID array;
};
