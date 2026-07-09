//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//

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
// Description: J2 plasticity material where stress
// components 22, 33, 13, and 23 are constrained to 0.
//
// Scott, Michael H. 
//   “Direct Differentiation of the J2 Plasticity Model under Assumed Stress States.” n.d.
//
// Written: CMP
// Created: Feb 2026
//
#include <J2BeamThread3d.h>
#include <Channel.h>
#include <string.h>
#include <cmath>
#include <float.h>
#include <MatrixND.h>
#include <VectorND.h>
#include <MaterialResponse.h>
#include <NonlinearJ2.h>
#include <domain/DomainStatus.h>

using namespace OpenSees;

static constexpr double one3   = 1.0 / 3.0;
static constexpr double two3   = 2.0 / 3.0;
static const     double root23 = std::sqrt(2.0 / 3.0);


//
//  Constructors / Destructor
//

J2BeamThread3d::J2BeamThread3d(int tag,
          double e, double g, 
          double sy, double hi, double hk,
          double rho
        )
 : NDMaterial(tag, ND_TAG_J2BeamThread3d),
   E(e),
   nu(g),
   sigmaY(sy),
   Hiso(hi),
   Hkin(hk),
   density(rho),
   parameterID(0),
   SHVs(nullptr),
   Tepsilon{},
   sigma{},
   D_wrapper(D),
   s_wrapper(sigma),
   e_wrapper(Tepsilon),
   alphan(0.0),
   alphan1(0.0),
   dg_n1(0.0),
   xsi_n1{},
   q_n1(0.0)
{
  epsPn[0]  = 0.0;  epsPn[1]  = 0.0;  epsPn[2]  = 0.0;
  epsPn1[0] = 0.0;  epsPn1[1] = 0.0;  epsPn1[2] = 0.0;
}

J2BeamThread3d::J2BeamThread3d()
 : NDMaterial(0, ND_TAG_J2BeamThread3d),
   E(0.0),
   nu(0.0),
   sigmaY(0.0),
   Hiso(0.0),
   Hkin(0.0),
   parameterID(0),
   SHVs(nullptr),
   Tepsilon{},
   sigma{},
   D_wrapper(D),
   s_wrapper(sigma),
   e_wrapper(Tepsilon),
   alphan(0.0),
   alphan1(0.0),
   dg_n1(0.0),
   xsi_n1{},
   q_n1(0.0)
{
  epsPn[0]  = 0.0;  epsPn[1]  = 0.0;  epsPn[2]  = 0.0;
  epsPn1[0] = 0.0;  epsPn1[1] = 0.0;  epsPn1[2] = 0.0;
}

J2BeamThread3d::~J2BeamThread3d()
{
  if (SHVs != nullptr)
    delete SHVs;
}

//
//  Strain interface
//

int
J2BeamThread3d::setTrialStrain(const Vector& strain)
{
  Tepsilon = strain;
  return static_cast<int>(returnMap());
}


int
J2BeamThread3d::setTrialStrainIncr(const Vector& strain)
{
  assert(false);
  return 0;
}

//
//  Tangent
//

const Matrix&
J2BeamThread3d::getTangent()
{
  const double twoG     = E / (1.0 + nu);
  const double G        = 0.5 * twoG;
  const double two3Hkin = two3 * Hkin;

  // Elastic tangent
  if (dg_n1 == 0.0) {
    D(0,0) = E;    D(0,1) = 0.0;  D(0,2) = 0.0;
    D(1,0) = 0.0;  D(1,1) = G;    D(1,2) = 0.0;
    D(2,0) = 0.0;  D(2,1) = 0.0;  D(2,2) = G;
    return D_wrapper;
  }

  // Consistent tangent for plastic step
  //
  // Build the condensed 4x4 system whose top-left 3x3 block, after
  // inversion and scaling, gives D_{ep}.
  const double dg = dg_n1;
  const double q  = q_n1;
  const double den = 1.0 + dg * two3Hkin;

  MatrixND<4,4> Jt{};

  Jt(0,0) = 1.0 + dg * two3 * E / den;
  Jt(1,1) = 1.0 + dg * twoG / den;
  Jt(2,2) = 1.0 + dg * twoG / den;

  Jt(0,3) = (two3 * E  - dg * two3 * E  / den * two3Hkin) * xsi_n1[0];
  Jt(1,3) = (twoG      - dg * twoG      / den * two3Hkin) * xsi_n1[1];
  Jt(2,3) = (twoG      - dg * twoG      / den * two3Hkin) * xsi_n1[2];

  double c = (q > 0.0) ? (1.0 - two3 * Hiso * dg) / (q * den) : 0.0;
  Jt(3,0) = c * xsi_n1[0] * two3;
  Jt(3,1) = c * xsi_n1[1] * 2.0;
  Jt(3,2) = c * xsi_n1[2] * 2.0;

  Jt(3,3) = -q * two3Hkin/den - two3 * Hiso * q;

  MatrixND<4,4> invJ;
  Jt.invert(invJ);

  D(0,0) = invJ(0,0) * E;   D(0,1) = invJ(0,1) * G;   D(0,2) = invJ(0,2) * G;
  D(1,0) = invJ(1,0) * E;   D(1,1) = invJ(1,1) * G;   D(1,2) = invJ(1,2) * G;
  D(2,0) = invJ(2,0) * E;   D(2,1) = invJ(2,1) * G;   D(2,2) = invJ(2,2) * G;

  return D_wrapper;
}

const Matrix&
J2BeamThread3d::getInitialTangent()
{
  double G = 0.5 * E / (1.0 + nu);

  D(0,0) = E;    D(0,1) = 0.0;  D(0,2) = 0.0;
  D(1,0) = 0.0;  D(1,1) = G;    D(1,2) = 0.0;
  D(2,0) = 0.0;  D(2,1) = 0.0;  D(2,2) = G;

  return D_wrapper;
}

//
//  Stress / Strain accessors
//

const Vector&
J2BeamThread3d::getStress()
{
  // sigma was already computed in returnMap (called from setTrialStrain)
  return s_wrapper;
}

const Vector&
J2BeamThread3d::getStrain()
{
  return e_wrapper;
}

//
//  State management
//

int
J2BeamThread3d::commitState()
{
  epsPn[0] = epsPn1[0];
  epsPn[1] = epsPn1[1];
  epsPn[2] = epsPn1[2];
  alphan   = alphan1;
  return 0;
}

int
J2BeamThread3d::revertToLastCommit()
{
  epsPn1[0] = epsPn[0];
  epsPn1[1] = epsPn[1];
  epsPn1[2] = epsPn[2];
  alphan1   = alphan;
  return 0;
}

int
J2BeamThread3d::revertToStart()
{
  Tepsilon.zero();

  epsPn[0]  = 0.0;  epsPn[1]  = 0.0;  epsPn[2]  = 0.0;
  epsPn1[0] = 0.0;  epsPn1[1] = 0.0;  epsPn1[2] = 0.0;

  alphan  = 0.0;
  alphan1 = 0.0;
  dg_n1   = 0.0;
  q_n1    = 0.0;
  xsi_n1  = {};

  if (SHVs != nullptr)
    SHVs->Zero();

  return 0;
}

//
//  Copy
//

NDMaterial*
J2BeamThread3d::getCopy()
{
  return new J2BeamThread3d(this->getTag(), 
                            E, nu, sigmaY, Hiso, Hkin, 
                            density);
}


NDMaterial*
J2BeamThread3d::getCopy(const char* type)
{
  if (strcmp(type, this->getType()) == 0)
    return this->getCopy();
  
  std::vector<double> Ck{Hkin};
  std::vector<double> gammak{0.0};

  NonlinearJ2 *parent = new NonlinearJ2(this->getTag(), E, nu, sigmaY, density, Hiso,
                                        0, 0, 0, 0, Ck, gammak);
  NDMaterial *copy = parent->getCopy(type);
  delete parent;
  return copy;
}

const char*
J2BeamThread3d::getType() const
{
  return "BeamFiber";
}

int
J2BeamThread3d::getOrder() const
{
  return 3;
}


#if 0
Response*
J2BeamThread3d::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse =0;
  const char *matType = this->getType();

  output.tag("NdMaterialOutput");
  output.attr("matType",this->getClassType());
  output.attr("matTag",this->getTag());

  if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0) {
    const Vector &res = this->getStress();
    int size = res.Size();
    
    if ((strcmp(matType,"PlaneStress") == 0 && size == 3) ||
        (strcmp(matType,"PlaneStrain") == 0 && size == 3)) {
      output.tag("ResponseType","sigma11");
      output.tag("ResponseType","sigma22");
      output.tag("ResponseType","sigma12");
        } else if (strcmp(matType,"ThreeDimensional") == 0 && size == 6) {
      output.tag("ResponseType","sigma11");
      output.tag("ResponseType","sigma22");
      output.tag("ResponseType","sigma33");
      output.tag("ResponseType","sigma12");
      output.tag("ResponseType","sigma23");
      output.tag("ResponseType","sigma13");
    } else {
      for (int i=0; i<size; i++) 
      output.tag("ResponseType","UnknownStress");
    }
    theResponse =  new MaterialResponse(this, 1, this->getStress());

  } else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0) {
    const Vector &res = this->getStrain();
    int size = res.Size();
    if ((strcmp(matType,"PlaneStress") == 0 && size == 3) ||
        (strcmp(matType,"PlaneStrain") == 0 && size == 3)) {
      output.tag("ResponseType","eta11");
      output.tag("ResponseType","eta22");
      output.tag("ResponseType","eta12");
    } else if (strcmp(matType,"ThreeDimensional") == 0 && size == 6) {
      output.tag("ResponseType","eps11");
      output.tag("ResponseType","eps22");
      output.tag("ResponseType","eps33");
      output.tag("ResponseType","eps12");
      output.tag("ResponseType","eps23");
      output.tag("ResponseType","eps13");
    } else {
      for (int i=0; i<size; i++) 
        output.tag("ResponseType","UnknownStrain");
    }      
    theResponse =  new MaterialResponse(this, 2, this->getStrain());
  }
  //Adding temperature and thermal expansion output,L.Jiang [SIF]
  else if (strcmp(argv[0], "TempAndElong") == 0 || strcmp(argv[0], "TempAndElong") == 0) {
	  const Vector &res = this->getTempAndElong();
	  int size = res.Size();
	  if (size == 2) {
		  output.tag("ResponseType", "Temp");
		  output.tag("ResponseType", "Elong");
	  }
	  //opserr<<"tempElong "<<this->getTempAndElong()<<endln;
	  theResponse = new MaterialResponse(this, 3, this->getTempAndElong());
  }
  //end of adding output request,L.Jiang [SIF]
  else if (strcmp(argv[0], "Tangent") == 0 || strcmp(argv[0], "tangent") == 0) {
	  const Matrix &res = this->getTangent();
	  theResponse = new MaterialResponse(this, 4, this->getTangent());
  }

  // Massimo Petracca - 28/12/2021:
  // this should be handled by the PlaneStressUserMaterial... not here!
  //default damage output - added by V.K. Papanikolaou [AUTh] - start
  //else if (strcmp(argv[0], "Damage") == 0 || strcmp(argv[0], "damage") == 0) {
  //    static Vector vec = Vector(3);
  //    for (int i = 0; i < 3; i++) vec[i] = 0;
  //    theResponse = new MaterialResponse(this, 5, vec);  // zero vector
  //}
  //default damage output - added by V.K. Papanikolaou [AUTh] - end 

  output.endTag(); // NdMaterialOutput

  return theResponse;
}

int 
J2BeamThread3d::getResponse(int id, Information &info)
{
  switch (id) {
  case int(Request::PlasticStrain):
    return info.setVector(this->getStress());
  default:
	  return NDMaterial::getResponse(id, info);
  }
}
#endif

//
//  Sensitivity
//
const Vector&
J2BeamThread3d::getStressSensitivity(int gradIndex, bool conditional)
{
  // This returns a reference to a local static — matches the API pattern
  static Vector dsigmadh(3);

  dsigmadh(0) = 0.0;
  dsigmadh(1) = 0.0;
  dsigmadh(2) = 0.0;

  double dEdh      = 0.0;
  double dsigmaYdh = 0.0;
  double dHkindh   = 0.0;
  double dHisodh   = 0.0;
  double dGdh      = 0.0;

  if (parameterID == 1) { dEdh = 1.0;  dGdh = 0.5 / (1.0 + nu); }
  if (parameterID == 2) { dGdh = -0.5 * E / ((1.0 + nu) * (1.0 + nu)); }
  if (parameterID == 5) { dsigmaYdh = 1.0; }
  if (parameterID == 6) { dHkindh   = 1.0; }
  if (parameterID == 7) { dHisodh   = 1.0; }

  const double G = 0.5 * E / (1.0 + nu);

  double depsPdh[3] = {0.0, 0.0, 0.0};
  double dalphadh   = 0.0;
  if (SHVs != nullptr) {
    depsPdh[0] = (*SHVs)(0, gradIndex);
    depsPdh[1] = (*SHVs)(1, gradIndex);
    depsPdh[2] = (*SHVs)(2, gradIndex);
    dalphadh   = (*SHVs)(3, gradIndex);
  }

  // Recompute xsi from epsPn1 (same as xsi_n1 stored from return map)
  double xsi[3];
  xsi[0] = E * (Tepsilon(0) - epsPn1[0]) -        Hkin * epsPn1[0];
  xsi[1] = G * (Tepsilon(1) - epsPn1[1]) - one3 * Hkin * epsPn1[1];
  xsi[2] = G * (Tepsilon(2) - epsPn1[2]) - one3 * Hkin * epsPn1[2];

  double q = std::sqrt(two3 * xsi[0]*xsi[0] + 2.0*xsi[1]*xsi[1] + 2.0*xsi[2]*xsi[2]);
  double F = q - root23 * (sigmaY + Hiso * alphan1);

  if (F <= -100.0 * DBL_EPSILON) {
    // Elastic
    dsigmadh(0) = dEdh * (Tepsilon(0) - epsPn1[0]) - E * depsPdh[0];
    dsigmadh(1) = dGdh * (Tepsilon(1) - epsPn1[1]) - G * depsPdh[1];
    dsigmadh(2) = dGdh * (Tepsilon(2) - epsPn1[2]) - G * depsPdh[2];
  } else {
    // Plastic
    const double dg       = dg_n1;
    const double twoG     = E / (1.0 + nu);
    const double two3Hkin = two3 * Hkin;

    MatrixND<4,4> J{};
    VectorND<4>   b{};
    VectorND<4>   dx{};

    J(0,0) = 1.0 + dg * two3 * (E + Hkin);
    J(1,1) = 1.0 + dg * (twoG + two3Hkin);
    J(2,2) = 1.0 + dg * (twoG + two3Hkin);

    J(0,3) = two3 * (E + Hkin)  * xsi[0];
    J(1,3) = (twoG + two3Hkin)  * xsi[1];
    J(2,3) = (twoG + two3Hkin)  * xsi[2];

    double c = (q > 0.0) ? (1.0 - two3 * Hiso * dg) / q : 0.0;
    J(3,0) = c * xsi[0] * two3;
    J(3,1) = c * xsi[1] * 2.0;
    J(3,2) = c * xsi[2] * 2.0;
    J(3,3) = -two3 * Hiso * q;

    b(0) = dEdh * Tepsilon(0) - (E + Hkin)         * depsPdh[0] - (dEdh + dHkindh)         * epsPn1[0];
    b(1) = dGdh * Tepsilon(1) - (G + one3 * Hkin)  * depsPdh[1] - (dGdh + one3 * dHkindh)  * epsPn1[1];
    b(2) = dGdh * Tepsilon(2) - (G + one3 * Hkin)  * depsPdh[2] - (dGdh + one3 * dHkindh)  * epsPn1[2];
    b(3) = root23 * (dsigmaYdh + dHisodh * alphan1 + Hiso * dalphadh);

    J.rsolve(b, dx);

    depsPdh[0] += dx(3) * two3 * xsi[0] + dg * two3 * dx(0);
    depsPdh[1] += dx(3) * 2.0  * xsi[1] + dg * 2.0  * dx(1);
    depsPdh[2] += dx(3) * 2.0  * xsi[2] + dg * 2.0  * dx(2);

    dsigmadh(0) = dx(0) +        Hkin * depsPdh[0] +        dHkindh * epsPn1[0];
    dsigmadh(1) = dx(1) + one3 * Hkin * depsPdh[1] + one3 * dHkindh * epsPn1[1];
    dsigmadh(2) = dx(2) + one3 * Hkin * depsPdh[2] + one3 * dHkindh * epsPn1[2];
  }

  return dsigmadh;
}

int
J2BeamThread3d::commitSensitivity(const Vector& depsdh, int gradIndex, int numGrads)
{
  if (SHVs == nullptr)
    SHVs = new Matrix(4, numGrads);

  if (gradIndex >= SHVs->noCols())
    return 0;

  double dEdh      = 0.0;
  double dsigmaYdh = 0.0;
  double dHkindh   = 0.0;
  double dHisodh   = 0.0;
  double dGdh      = 0.0;

  if (parameterID == 1) { dEdh = 1.0;  dGdh = 0.5 / (1.0 + nu); }
  if (parameterID == 2) { dGdh = -0.5 * E / ((1.0 + nu) * (1.0 + nu)); }
  if (parameterID == 5) { dsigmaYdh = 1.0; }
  if (parameterID == 6) { dHkindh   = 1.0; }
  if (parameterID == 7) { dHisodh   = 1.0; }

  const double G = 0.5 * E / (1.0 + nu);

  double depsPdh[3] = {0.0, 0.0, 0.0};
  double dalphadh   = 0.0;
  if (SHVs != nullptr) {
    depsPdh[0] = (*SHVs)(0, gradIndex);
    depsPdh[1] = (*SHVs)(1, gradIndex);
    depsPdh[2] = (*SHVs)(2, gradIndex);
    dalphadh   = (*SHVs)(3, gradIndex);
  }

  double xsi[3];
  xsi[0] = E * (Tepsilon(0) - epsPn1[0]) -        Hkin * epsPn1[0];
  xsi[1] = G * (Tepsilon(1) - epsPn1[1]) - one3 * Hkin * epsPn1[1];
  xsi[2] = G * (Tepsilon(2) - epsPn1[2]) - one3 * Hkin * epsPn1[2];

  double q = std::sqrt(two3 * xsi[0]*xsi[0] + 2.0*xsi[1]*xsi[1] + 2.0*xsi[2]*xsi[2]);
  double F = q - root23 * (sigmaY + Hiso * alphan1);

  if (F <= -100.0 * DBL_EPSILON) {
    // Elastic, no updates
  } else {
    // Plastic
    const double dg       = dg_n1;
    const double twoG     = E / (1.0 + nu);
    const double two3Hkin = two3 * Hkin;

    MatrixND<4,4> J{};
    VectorND<4>   b{};
    VectorND<4>   dx{};

    J(0,0) = 1.0 + dg * two3 * (E + Hkin);
    J(1,1) = 1.0 + dg * (twoG + two3Hkin);
    J(2,2) = 1.0 + dg * (twoG + two3Hkin);

    J(0,3) = two3 * (E + Hkin)  * xsi[0];
    J(1,3) = (twoG + two3Hkin)  * xsi[1];
    J(2,3) = (twoG + two3Hkin)  * xsi[2];

    double c = (q > 0.0) ? (1.0 - two3 * Hiso * dg) / q : 0.0;
    J(3,0) = c * xsi[0] * two3;
    J(3,1) = c * xsi[1] * 2.0;
    J(3,2) = c * xsi[2] * 2.0;
    J(3,3) = -two3 * Hiso * q;

    b(0) = E * depsdh(0) + dEdh * Tepsilon(0) - (E + Hkin)        * depsPdh[0] - (dEdh + dHkindh)        * epsPn1[0];
    b(1) = G * depsdh(1) + dGdh * Tepsilon(1) - (G + one3 * Hkin) * depsPdh[1] - (dGdh + one3 * dHkindh) * epsPn1[1];
    b(2) = G * depsdh(2) + dGdh * Tepsilon(2) - (G + one3 * Hkin) * depsPdh[2] - (dGdh + one3 * dHkindh) * epsPn1[2];
    b(3) = root23 * (dsigmaYdh + dHisodh * alphan1 + Hiso * dalphadh);

    J.rsolve(b, dx);

    dalphadh += dx(3) * root23 * q
              + dg * root23 * (xsi[0]*two3*dx(0) + xsi[1]*2.0*dx(1) + xsi[2]*2.0*dx(2)) / q;

    depsPdh[0] += dx(3) * two3 * xsi[0] + dg * two3 * dx(0);
    depsPdh[1] += dx(3) * 2.0  * xsi[1] + dg * 2.0  * dx(1);
    depsPdh[2] += dx(3) * 2.0  * xsi[2] + dg * 2.0  * dx(2);

    (*SHVs)(0, gradIndex) = depsPdh[0];
    (*SHVs)(1, gradIndex) = depsPdh[1];
    (*SHVs)(2, gradIndex) = depsPdh[2];
    (*SHVs)(3, gradIndex) = dalphadh;
  }

  return 0;
}

//
//  Serialization
//

int
J2BeamThread3d::sendSelf(int commitTag, Channel& theChannel)
{
  static Vector data(6);

  data(0) = this->getTag();
  data(1) = E;
  data(2) = nu;
  data(3) = sigmaY;
  data(4) = Hiso;
  data(5) = Hkin;

  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0)
    opserr << "J2BeamThread3d::sendSelf -- could not send Vector\n";

  return res;
}

int
J2BeamThread3d::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
  static Vector data(6);

  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "J2BeamThread3d::recvSelf -- could not recv Vector\n";
    return res;
  }

  this->setTag((int)data(0));
  E      = data(1);
  nu     = data(2);
  sigmaY = data(3);
  Hiso   = data(4);
  Hkin   = data(5);

  return res;
}

//
//  Print
//

void
J2BeamThread3d::Print(OPS_Stream& s, int flag)
{
  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "J2 Beam Fiber Material Model" << "\n";
    s << "\tE:  " << E << "\n";
    s << "\tnu:  " << nu << "\n";
    s << "\tsigmaY:  " << sigmaY << "\n";
    s << "\tHiso:  " << Hiso << "\n";
    s << "\tHkin:  " << Hkin << "\n";
  }
  else if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"J2BeamFiber\", ";
    s << "\"E\": " << E << ", ";
    s << "\"nu\": " << nu << ", ";
    s << "\"sigmaY\": " << sigmaY << ", ";
    s << "\"Hiso\": " << Hiso << ", ";
    s << "\"Hkin\": " << Hkin;
    s << "}";
  }
}

//
//  Parameter interface
//

int
J2BeamThread3d::setParameter(const char** argv, int argc, Parameter& param)
{
  if (strcmp(argv[0], "E") == 0) {
    param.setValue(E);
    return param.addObject(1, this);
  } else if (strcmp(argv[0], "nu") == 0) {
    param.setValue(nu);
    return param.addObject(2, this);
  } else if (strcmp(argv[0], "sigmaY") == 0 || strcmp(argv[0], "fy") == 0 ||
             strcmp(argv[0], "Fy") == 0) {
    param.setValue(sigmaY);
    return param.addObject(5, this);
  } else if (strcmp(argv[0], "Hkin") == 0) {
    param.setValue(Hkin);
    return param.addObject(6, this);
  } else if (strcmp(argv[0], "Hiso") == 0) {
    param.setValue(Hiso);
    return param.addObject(7, this);
  }
  return -1;
}

int
J2BeamThread3d::updateParameter(int parameterID, Information& info)
{
  switch (parameterID) {
  case 1:  E      = info.theDouble; return 0;
  case 2:  nu     = info.theDouble; return 0;
  case 5:  sigmaY = info.theDouble; return 0;
  case 6:  Hkin   = info.theDouble; return 0;
  case 7:  Hiso   = info.theDouble; return 0;
  default: return -1;
  }
}

int
J2BeamThread3d::activateParameter(int paramID)
{
  parameterID = paramID;
  return 0;
}


//
//  On exit the following member variables are set:
//    epsPn1[0..2], alphan1, dg_n1, xsi_n1[0..2], q_n1, sigma[0..2]
//
DomainStatus
J2BeamThread3d::returnMap() noexcept
{
  const double twoG    = E / (1.0 + nu);
  const double G       = 0.5 * twoG;
  const double two3Hkin = two3 * Hkin;

  // Trial (elastic) stress
  Vector3D sig;
  sig[0] = E * (Tepsilon(0) - epsPn[0]);
  sig[1] = G * (Tepsilon(1) - epsPn[1]);
  sig[2] = G * (Tepsilon(2) - epsPn[2]);

  // Trial relative stress  xsi = sig - back-stress
  Vector3D xsi;
  xsi[0] = sig[0] -        Hkin * epsPn[0];
  xsi[1] = sig[1] - one3 * Hkin * epsPn[1];
  xsi[2] = sig[2] - one3 * Hkin * epsPn[2];

  double q = std::sqrt(two3 * xsi[0]*xsi[0]
                     + 2.0  * xsi[1]*xsi[1]
                     + 2.0  * xsi[2]*xsi[2]);

  double F = q - root23 * (sigmaY + Hiso * alphan);

  // Elastic step
  if (F < 0.0) { // -100.0 * DBL_EPSILON) {
    epsPn1[0] = epsPn[0];
    epsPn1[1] = epsPn[1];
    epsPn1[2] = epsPn[2];
    alphan1   = alphan;
    dg_n1     = 0.0;

    xsi_n1    = xsi;
    q_n1      = q;

    sigma(0)  = sig[0];
    sigma(1)  = sig[1];
    sigma(2)  = sig[2];
    return DomainStatus::Success;
  }

  // Plastic step — local Newton on  {xsi, dg}  system
  double dg = 0.0;

  VectorND<4> R{0.0, 0.0, 0.0, F};
  VectorND<4> x{xsi[0], xsi[1], xsi[2], dg};

  MatrixND<4,4> J{};
  VectorND<4>   dx{};

  int iter    = 0;
  int maxIter = 25;
  constexpr double tol = 1.0e-12; // 10 14

  const double tolFy = sigmaY * tol;

  do {
    J(0,0) = 1.0 + dg * two3 * (E + Hkin);
    J(0,1) = 0.0;
    J(0,2) = 0.0;
    J(1,0) = 0.0;
    J(1,1) = 1.0 + dg * (twoG + two3Hkin);
    J(1,2) = 0.0;
    J(2,0) = 0.0;
    J(2,1) = 0.0;
    J(2,2) = 1.0 + dg * (twoG + two3Hkin);

    J(0,3) = two3 * (E + Hkin) * x(0);
    J(1,3) = (twoG + two3Hkin) * x(1);
    J(2,3) = (twoG + two3Hkin) * x(2);

    double c = (q > 0.0) ? (1.0 - two3 * Hiso * dg) / q : 0.0;
    J(3,0) = c * x[0] * two3;
    J(3,1) = c * x[1] * 2.0;
    J(3,2) = c * x[2] * 2.0;
    J(3,3) = -two3 * Hiso * q;

    J.rsolve(R, dx);
    x.addVector(1.0, dx, -1.0);

    dg = x(3);
    q  = std::sqrt(two3 * x(0)*x(0) + 2.0*x(1)*x(1) + 2.0*x(2)*x(2));

    R(0) = x(0) - xsi[0] + dg * two3 * (E + Hkin)  * x(0);
    R(1) = x(1) - xsi[1] + dg * (twoG + two3Hkin)  * x(1);
    R(2) = x(2) - xsi[2] + dg * (twoG + two3Hkin)  * x(2);
    R(3) = q - root23 * (sigmaY + Hiso * (alphan + dg * root23 * q));

  } while ((R.norm() > tolFy) && (iter++ < maxIter));

  if (R.norm() > tolFy) {
    return DomainStatus::MaterialFailedToConverge;
  }


  // Store converged state
  dg_n1     = dg;
  alphan1   = alphan + dg * root23 * q;

  epsPn1[0] = epsPn[0] + dg * two3 * x(0);
  epsPn1[1] = epsPn[1] + dg * 2.0  * x(1);
  epsPn1[2] = epsPn[2] + dg * 2.0  * x(2);

  xsi_n1[0] = x(0);
  xsi_n1[1] = x(1);
  xsi_n1[2] = x(2);
  q_n1      = q;

  // Stress = relative stress + back-stress
  sigma(0) = x(0) +        Hkin * epsPn1[0];
  sigma(1) = x(1) + one3 * Hkin * epsPn1[1];
  sigma(2) = x(2) + one3 * Hkin * epsPn1[2];

  return DomainStatus::Success;
}
