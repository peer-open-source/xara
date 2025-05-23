//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
//
// This material implements a Bouc-Wen hysteretic material with both Pinching 
// and Strength degradation as described by Foliente. This is essentially a 
// combination of the preexisting BWBN and BoucWenMaterial classes.
//
// References
//
// - Prof. Haukaas' notes
//   https://civil-terje.sites.olt.ubc.ca/files/2023/03/Bouc-Wen-Material-Model.pdf
//
// Written: cmp
// Created: May 2025
//
#include "BWBF.h"

#include <Vector.h>
#include <Channel.h>
#include <Logging.h>
#include <math.h>
#include <Matrix.h>
#include <Information.h>
#include <Parameter.h>

#define MAT_TAG_BWBF  1000


BWBF::BWBF(int tag, 
          double E,
          double Fy,
          double alpha,
          double n,
          double beta,
          double delta_a,
          double delta_v,
          double delta_n,
          double pinch_slope,
          double pinch_slip,
          double pinch_start,
          double pinch_rate,
          double pinch_size,
          double pinch_lamda,
          double tolerance,
          int MaxNumIter)
 : UniaxialMaterial(tag,MAT_TAG_BWBF),
  Fy(Fy), Ko(E), alpha(alpha), 
  n(n), beta(beta), gamma(1-beta), Ao(1.0), 
  delta_a(delta_a), delta_v(delta_v), delta_n(delta_n),
  pinch(pinch_start, pinch_slip, pinch_slope, pinch_size, pinch_rate, pinch_lamda),
  tolerance(tolerance),
  maxNumIter(MaxNumIter)
{
  this->revertToStart();
}


BWBF::~BWBF()
{

}


int 
BWBF::setTrialStrain(double strain, double strainRate)
{

  pres.strain = strain;
  const double dStrain = pres.strain - past.strain;
  const double sgn  = signum(dStrain);

  // Newton-Raphson scheme to solve for z(strain)
  double z = past.z; // 0.01;
  int count = 0;
  
  double W, a, h, u, az, hz, ue, e;
  double Ae, nu, eta;
  double g,  Dz = 0;
  do {
    z -= Dz;

    e   =  past.energy + (1.0-alpha)*(Ko/Fy)*z*dStrain;
    Ae  =  Ao - delta_a * e;
    nu  = 1.0 + delta_v * e;
    eta = 1.0 + delta_n * e;

    double Psi = gamma + beta*signum(dStrain*z);
    W = wen(z, Psi, Ae, nu);
    a = W/eta;
    h = pinch.update(sgn, z, e, u);
    u = std::pow(Ae/(nu*(beta + gamma)), 1./n);
    
    double ve = delta_v;
    ue = -(beta + gamma)/n*std::pow(Ae/(nu*(beta + gamma)), 1+1./n)*ve;

    // set hz, az;
    {
      double ez   = (1.0 - alpha)*(Ko/Fy)*dStrain;
      double vz   =   delta_v * ez;
      double nz   =   delta_n * ez;
      double Az   = - delta_a * ez;
      double pow1 = (z==0.0) ? 0.0 : std::pow(std::fabs(z), (n-1));
      double Wz = dwen(z, Psi, Ae, nu, 1, Az, vz);
      hz = pinch.tangent(ez,ue,1);
      az = (Wz*eta - W*nz)/(eta*eta);
    }

    // Newton-Raphson objective function
    g = (z - past.z) - (h*a)*dStrain;

    double gz = 1.0 - (az*h + a*hz)*dStrain;

    // Issue warning if derivative is zero
    if (std::fabs(gz) < 1.0e-10 ) {
        opserr << "WARNING: zero derivative in BWBF Newton-Raphson scheme\n";
        return -1;
    }

    // Take a Newton step
    Dz = g/gz;

    // Update counter
    count++;

    if (count == maxNumIter) {
        opserr << "WARNING: Failed to find root after " << maxNumIter << " iterations"
                << " with norm: " << std::fabs(g) << "\n";
        return -2;
    }

  } while ( ( std::fabs(g) > tolerance ) && count<maxNumIter);


  pres.z = z;
  pres.energy = e;
  //
  // Compute stress
  //
  Tstress = alpha*Ko*pres.strain + (1-alpha)*Fy*z;

  //
  // Compute tangent
  //
  double dzdx;
  if (z == 0.0) 
    dzdx = Ko/Fy;
  else {
    double Psi = gamma + beta*signum(dStrain*pres.z);

    double ax, hx;
    {
      double ex = (1-alpha)*(Ko/Fy)*pres.z; // b1
      double vx =  delta_v * ex;
      double Ax = -delta_a * ex;
      double wx = dwen(z, Psi, Ae, nu,  0, Ax, vx);
      double nx = delta_n * ex;
      ax = (eta*wx - W*nx)/(eta*eta);
      hx = pinch.tangent(ex,ue,0);
    }

    double f  = h*a;
    double fx = hx*a + h*ax;
    double fz = az*h + a*hz;

    dzdx = (f + fx*dStrain)/(1.0 - fz*dStrain);
  }
  pres.tangent = alpha*Ko + (1-alpha)*Fy*dzdx;

  return 0;
}

double 
BWBF::getStress()
{
    return Tstress;
}

double 
BWBF::getInitialTangent()
{
  return alpha*Ko + (1-alpha)*Ko*Ao;
}


double 
BWBF::getTangent()
{
  return pres.tangent;
}

double 
BWBF::getStrain()
{
  return pres.strain;
}

int 
BWBF::commitState()
{
  // Commit trial history variables
  past = pres;
  return 0;
}

int 
BWBF::revertToLastCommit()
{
  return 0;
}

int 
BWBF::revertToStart()
{
  past.z = 0.0;
  past.energy = 0.0;
  past.tangent = Ko;
  pres.strain = 0.0;
  past.strain = 0.0;
  pres.z = 0.01;
  pres.energy = 0.0;
  //
  Tstress = 0.0;
  pres.tangent = alpha*Ko + (1-alpha)*(Ko/Fy)*Ao;

  return 0;
}

UniaxialMaterial *
BWBF::getCopy()
{
  BWBF *theCopy =
  new BWBF(this->getTag(), Ko, Fy, alpha, n,
            beta, 
            delta_a, delta_v, delta_n,
            0,0,0,0,0,0,
            tolerance, maxNumIter);
      
  (*theCopy) = (*this);

  return theCopy;
}


int
BWBF::setParameter(const char **argv, int argc, Parameter &param)
{
  if (strcmp(argv[0],"alpha") == 0)
    return param.addObject(1, this);
  
  if (strcmp(argv[0],"E") == 0)
    return param.addObject(2, this);
  
  if (strcmp(argv[0],"n") == 0)
    return param.addObject(3, this);
  
  if (strcmp(argv[0],"gamma") == 0)
    return param.addObject(4, this);
    
  if (strcmp(argv[0],"beta") == 0)
    return param.addObject(5, this);
  
  if (strcmp(argv[0],"Ao") == 0)
    return param.addObject(6, this);
  
  if (strcmp(argv[0],"deltaA") == 0)
    return param.addObject(7, this);
    
  if (strcmp(argv[0],"deltaNu") == 0)
    return param.addObject(8, this);
  
  if (strcmp(argv[0],"deltaEta") == 0)
    return param.addObject(9, this);

  if (strcmp(argv[0],"Fy") == 0)
    return param.addObject(10, this);

  return -1;
}


int
BWBF::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
  case 1:
      this->alpha = info.theDouble;
      return 0;
  case 2:
      this->Ko = info.theDouble;
      return 0;
  case 3:
      this->n = info.theDouble;
      return 0;
  case 4:
      this->gamma = info.theDouble;
      return 0;
  case 5:
      this->beta = info.theDouble;
      return 0;
  case 6:
      this->Ao = info.theDouble;
      return 0;
  case 7:
      this->delta_a = info.theDouble;
      return 0;
  case 8:
      this->delta_v = info.theDouble;
      return 0;
  case 9:
      this->delta_n = info.theDouble;
      return 0;
  case 10:
      this->Fy = info.theDouble;
      return 0;
  default:
      return -1;
  }
}

int
BWBF::activateParameter(int p)
{
  parameterID = p;

  return 0;
}

int 
BWBF::sendSelf(int cTag, Channel &)
{
  return -1;
}

int 
BWBF::recvSelf(int cTag, Channel &, FEM_ObjectBroker &)
{
  return -1;
}

void 
BWBF::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": " << "\"BWBF\", ";
    s << "\"Fy\": " << Fy << ", ";
    s << "\"alpha\": " << alpha << ", ";
    s << "\"E\": " << Ko << ", ";
    s << "\"n\": " << n << ", ";
    s << "\"gamma\": " << gamma << ", ";
    s << "\"beta\": " << beta << ", ";
    s << "\"Ao\": " << Ao << ", ";
    
    s << "\"delta_a\": " << delta_a << ", ";
    s << "\"delta_v\": " << delta_v << ", ";
    s << "\"delta_n\": " << delta_n << ", ";

    s << "\"pinch\": {";
    s << "\"slope\": " << pinch.p << ", ";
    s << "\"slip\": "  << pinch.zetas << ", ";
    s << "\"start\": " << pinch.q << ", ";
    s << "\"size\": "  << pinch.psi << ", ";
    s << "\"rate\": "  << pinch.delta_psi << ", ";
    s << "\"lamda\": " << pinch.lamda;
    s << "}";
    s << "}";
    return;

  }
}


