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
#include <ElasticIsotropicBeamFiber.h>           
#include <Channel.h>
#include <string.h>
#include <cassert>
using namespace OpenSees;

Matrix ElasticIsotropicBeamFiber::D(3,3);

ElasticIsotropicBeamFiber::ElasticIsotropicBeamFiber
(int tag, double E, double nu, double rho):
  ElasticIsotropicMaterial (tag, ND_TAG_ElasticIsotropicBeamFiber, E, nu, rho),
  Tepsilon(),
  retStrain(Tepsilon)
{

}

ElasticIsotropicBeamFiber::ElasticIsotropicBeamFiber():
  ElasticIsotropicMaterial (0, ND_TAG_ElasticIsotropicBeamFiber, 0.0, 0.0),
  Tepsilon(),
  retStrain(Tepsilon)
{

}

ElasticIsotropicBeamFiber::~ElasticIsotropicBeamFiber ()
{

}

int
ElasticIsotropicBeamFiber::setTrialStrain(const Vector &strain)
{
  assert(strain.Size() == 3);

  for (int i=0; i<3; i++)
    Tepsilon[i] = strain(i);

  return 0;
}

int
ElasticIsotropicBeamFiber::setTrialStrain(const Vector &strain, const Vector &rate)
{
  assert(strain.Size() == 3);
  for (int i=0; i<3; i++)
    Tepsilon[i] = strain(i);

  return 0;
}

int
ElasticIsotropicBeamFiber::setTrialStrainIncr(const Vector &strain)
{
  return 0;
}

int
ElasticIsotropicBeamFiber::setTrialStrainIncr(const Vector &strain, const Vector &rate)
{
  return 0;
}

const Matrix&
ElasticIsotropicBeamFiber::getTangent()
{
  const double G = 0.5*E/(1.0 + v);

  D(0,0) = E;
  D(1,1) = G;
  D(2,2) = G;
  return D;
}

const Matrix&
ElasticIsotropicBeamFiber::getInitialTangent()
{
  double mu = 0.5*E/(1.0+v);

  D(0,0) = E;
  D(1,1) = mu;
  D(2,2) = mu;
  
  return D;
}

const Vector&
ElasticIsotropicBeamFiber::getStress()
{
  double mu = 0.5*E/(1.0+v);
  
  static Vector sigma(3);
  sigma(0) =  E*Tepsilon(0);
  sigma(1) = mu*Tepsilon(1);
  sigma(2) = mu*Tepsilon(2);
  
  return sigma;
}

const Vector&
ElasticIsotropicBeamFiber::getStrain()
{
  return retStrain;
}

int
ElasticIsotropicBeamFiber::commitState()
{
  return 0;
}

int
ElasticIsotropicBeamFiber::revertToLastCommit()
{
  return 0;
}

int
ElasticIsotropicBeamFiber::revertToStart()
{
  return 0;
}


NDMaterial*
ElasticIsotropicBeamFiber::getCopy()
{
  ElasticIsotropicBeamFiber *theCopy =
    new ElasticIsotropicBeamFiber(this->getTag(), E, v, rho);

  return theCopy;
}

const char*
ElasticIsotropicBeamFiber::getType() const
{
  return "BeamFiber";
}

int
ElasticIsotropicBeamFiber::getOrder() const
{
  return 3;
}

const Vector&
ElasticIsotropicBeamFiber::getStressSensitivity(int gradIndex, bool conditional)
{
  static Vector sigma(3);
  sigma(0) = 0.0;
  sigma(1) = 0.0;
  sigma(2) = 0.0;

  if (parameterID == 1) { // E
    //double mu = 0.5*E/(1.0+v);
    double dmudE = 0.5/(1.0+v);
    sigma(0) = Tepsilon(0);
    sigma(1) = dmudE*Tepsilon(1);
    sigma(2) = dmudE*Tepsilon(2);
  }
  if (parameterID == 2) { // nu
    //double mu = 0.5*E/(1.0+v);
    double dmudnu = -0.5*E/(1.0 + 2*v + v*v);

    sigma(0) = 0.0;
    sigma(1) = dmudnu*Tepsilon(1);
    sigma(2) = dmudnu*Tepsilon(2);
  }

  return sigma;
}
