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
#include <ElasticIsotropicPlaneStrain2D.h>                                                                        
#include <Channel.h>
Vector ElasticIsotropicPlaneStrain2D::sigma(3);
Matrix ElasticIsotropicPlaneStrain2D::D(3,3);

ElasticIsotropicPlaneStrain2D::ElasticIsotropicPlaneStrain2D
(int tag, double E, double nu, double rho) :
 ElasticIsotropicMaterial (tag, ND_TAG_ElasticIsotropicPlaneStrain2d, E, nu, rho),
 epsilon(3), Cepsilon(3)
{
  epsilon.Zero();
  Cepsilon.Zero();
}

ElasticIsotropicPlaneStrain2D::ElasticIsotropicPlaneStrain2D():
 ElasticIsotropicMaterial (0, ND_TAG_ElasticIsotropicPlaneStrain2d, 0.0, 0.0),
 epsilon(3), Cepsilon(3)
{
  epsilon.Zero();
  Cepsilon.Zero();
}

ElasticIsotropicPlaneStrain2D::~ElasticIsotropicPlaneStrain2D ()
{

}

int
ElasticIsotropicPlaneStrain2D::setTrialStrain(const Vector &strain)
{
  epsilon = strain;
  return 0;
}

int
ElasticIsotropicPlaneStrain2D::setTrialStrain(const Vector &strain, const Vector &rate)
{
  epsilon = strain;
  return 0;
}

int
ElasticIsotropicPlaneStrain2D::setTrialStrainIncr (const Vector &strain)
{
  epsilon += strain;
  return 0;
}

int
ElasticIsotropicPlaneStrain2D::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  epsilon += strain;
  return 0;
}

const Matrix&
ElasticIsotropicPlaneStrain2D::getTangent()
{
  return getInitialTangent();
}

const Matrix&
ElasticIsotropicPlaneStrain2D::getInitialTangent()
{
  double mu2 = E/(1.0+v);
  double lam = v*mu2/(1.0-2.0*v);
  double mu = 0.50*mu2;
  
  D(0,0) = D(1,1) = mu2+lam;
  D(0,1) = D(1,0) = lam;
  D(2,2) = mu;

  return D;
}

const Vector&
ElasticIsotropicPlaneStrain2D::getStress()
{
  double mu2 = E/(1.0+v);
  double lam = v*mu2/(1.0-2.0*v);
  double mu = 0.50*mu2;

  double eps0 = epsilon(0);
  double eps1 = epsilon(1);

  mu2 += lam;

  //sigma = D*epsilon;
  sigma(0) = mu2*eps0 + lam*eps1;
  sigma(1) = lam*eps0 + mu2*eps1;
  sigma(2) = mu*epsilon(2);
	
  return sigma;
}

const Vector&
ElasticIsotropicPlaneStrain2D::getStrain()
{
  return epsilon;
}

int
ElasticIsotropicPlaneStrain2D::commitState()
{
  Cepsilon=epsilon;
  return 0;
}

int
ElasticIsotropicPlaneStrain2D::revertToLastCommit()
{
  epsilon=Cepsilon;
  return 0;
}

int
ElasticIsotropicPlaneStrain2D::revertToStart()
{
  epsilon.Zero();
  Cepsilon.Zero();
  return 0;
}

NDMaterial*
ElasticIsotropicPlaneStrain2D::getCopy()
{
  ElasticIsotropicPlaneStrain2D *theCopy =
    new ElasticIsotropicPlaneStrain2D (this->getTag(), E, v, rho);
  
  theCopy->epsilon = epsilon;
  theCopy->Cepsilon = Cepsilon;
  return theCopy;
}

const char*
ElasticIsotropicPlaneStrain2D::getType() const
{
  return "PlaneStrain";
}

int
ElasticIsotropicPlaneStrain2D::getOrder() const
{
  return 3;
}

int 
ElasticIsotropicPlaneStrain2D::sendSelf(int commitTag, Channel &theChannel)
{
  
  static Vector data(7);
  
  data(0) = this->getTag();
  data(1) = E;
  data(2) = v;
  data(3) = rho;
  data(4) = Cepsilon(0);
  data(5) = Cepsilon(1);
  data(6) = Cepsilon(2);
  
  int res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    // opserr << "ElasticIsotropicPlaneStrain2D::sendSelf -- could not send Vector\n";
    return res;
  }
  
  return res;
}

int 
ElasticIsotropicPlaneStrain2D::recvSelf(int commitTag, Channel &theChannel, 
					FEM_ObjectBroker &theBroker)
{
  static Vector data(7);
  
  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    // opserr << "ElasticIsotropicPlaneStrain2D::sendSelf -- could not send Vector\n";
    return res;
  }

  this->setTag((int)data(0));
  E = data(1);
  v = data(2);
  rho = data(3);
  epsilon(0)=data(4);
  epsilon(1)=data(5);
  epsilon(2)=data(6);
  Cepsilon = epsilon;

  return res;
}
