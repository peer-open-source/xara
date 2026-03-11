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
//===----------------------------------------------------------------------===//
//
// Description: Elastic isotropic model where stress components 22, 33, and 23
// are condensed out.
//
#include <ElasticIsotropicBeamThread.h>           
#include <Channel.h>
#include <string.h>
#include <cassert>
using namespace OpenSees;

ElasticIsotropicBeamThread::ElasticIsotropicBeamThread(ElasticIsotropicMaterial& builder_):
  NDMaterial(builder_.getTag(), ND_TAG_ElasticIsotropicBeamFiber),
  builder(builder_),
  rho(builder_.rho),
  stress{},
  retStress(stress),
  Tepsilon{},
  retStrain(Tepsilon),
  D(std::make_shared<Matrix3D>()),
  retTangent(*D),
  parameterID(0)
{
  const double E = builder_.E;
  const double nu = builder_.v;
  const double G = 0.5*E/(1.0 + nu);

  (*D)(0,0) = E;
  (*D)(1,1) = G;
  (*D)(2,2) = G;
}


ElasticIsotropicBeamThread::ElasticIsotropicBeamThread(const ElasticIsotropicBeamThread &other)
  : NDMaterial(other.getTag(), ND_TAG_ElasticIsotropicBeamFiber),
    rho(other.rho),
    builder(other.builder),
    stress{},
    retStress(stress),
    Tepsilon{},
    retStrain(Tepsilon),
    D(other.D),
    retTangent(*D),
    parameterID(0)
{

}


ElasticIsotropicBeamThread::~ElasticIsotropicBeamThread()
{

}


int
ElasticIsotropicBeamThread::setTrialStrain(const Vector &strain)
{
  assert(strain.Size() == 3);

  for (int i=0; i<3; i++)
    Tepsilon[i] = strain(i);

  return 0;
}

int
ElasticIsotropicBeamThread::setTrialStrainIncr(const Vector &strain)
{
  return 0;
}

const Matrix&
ElasticIsotropicBeamThread::getTangent()
{
  return retTangent;
}


const Matrix&
ElasticIsotropicBeamThread::getInitialTangent()
{
  return retTangent;
}


const Vector&
ElasticIsotropicBeamThread::getStress()
{
  const double G = (*D)(1,1);
  const double E = (*D)(0,0);

  stress(0) = E*Tepsilon(0);
  stress(1) = G*Tepsilon(1);
  stress(2) = G*Tepsilon(2);
  return retStress;
}


const Vector&
ElasticIsotropicBeamThread::getStrain()
{
  return retStrain;
}

int
ElasticIsotropicBeamThread::commitState()
{
  return 0;
}

int
ElasticIsotropicBeamThread::revertToLastCommit()
{
  return 0;
}

int
ElasticIsotropicBeamThread::revertToStart()
{
  return 0;
}


NDMaterial*
ElasticIsotropicBeamThread::getCopy()
{
  return new ElasticIsotropicBeamThread(*this);
}



NDMaterial*
ElasticIsotropicBeamThread::getCopy(const char *type)
{
  return builder.getCopy(type);
}


const char*
ElasticIsotropicBeamThread::getType() const
{
  return "BeamFiber";
}

int
ElasticIsotropicBeamThread::getOrder() const
{
  return 3;
}


int
ElasticIsotropicBeamThread::activateParameter(int paramID)
{
  parameterID = paramID;
  return 0;
}


int
ElasticIsotropicBeamThread::setParameter(const char **argv, int argc, Parameter &param)
{
  const double E = (*D)(0,0);
  const double G = (*D)(1,1);
  const double v = 0.5*E/G - 1.0;

  if (strcmp(argv[0],"E") == 0) {
    param.setValue(E);
    return param.addObject(1, this);
  }
  else if (strcmp(argv[0],"nu") == 0 || strcmp(argv[0],"v") == 0) {
    param.setValue(v);
    return param.addObject(2, this);
  }
  else if (strcmp(argv[0],"rho") == 0) {
    param.setValue(rho);
    return param.addObject(3, this);
  }

  return -1;
}

const Vector&
ElasticIsotropicBeamThread::getStressSensitivity(int gradIndex, bool conditional)
{
  static Vector sigma(3);
  sigma(0) = 0.0;
  sigma(1) = 0.0;
  sigma(2) = 0.0;
  const double E = (*D)(0,0);
  const double G = (*D)(1,1);
  const double v = 0.5*E/G - 1.0;

  if (parameterID == 1) { // E
    //double mu = 0.5*E/(1.0+v);
    double dmudE = G/E;
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
