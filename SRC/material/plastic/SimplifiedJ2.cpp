//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// 3D J2 plasticity model with linear isotropic and kinematic hardening
//
// Written: Quan Gu and Zhijian Qiu
// Created: 2013.7
// Reference:
//   JP Conte, MK. Jagannath, \Seismic relibility analysis of concrete
//     gravity dams, A Report on Research, Rice University, 1995.
//
// -------------------
//
#include <math.h>
#include <stdlib.h>
#include <ID.h>
#include <Matrix.h>
#include <MatrixND.h>
#include <T2Vector.h>
#include <Vector.h>
#include <Projector.hh>
#include <SimplifiedJ2.h>
#include <Parameter.h>
#include <Information.h>
#include <MaterialResponse.h>

using namespace OpenSees;

Vector SimplifiedJ2::tmpVector(6);

// --- element: eps(1,1),eps(2,2),eps(3,3),2*eps(1,2),2*eps(2,3),2*eps(1,3) ----
// --- material strain: eps(1,1),eps(2,2),eps(3,3),eps(1,2),eps(2,3),eps(1,3) , same sign ----


SimplifiedJ2::SimplifiedJ2(int pTag,
                           int pNd,
                           double pG,
                           double pK,
                           double pSigmaY0,
                           double pH_kin,
                           double pH_iso,
                           double density)
 : NDMaterial(pTag, ND_TAG_SimplifiedJ2),
   density(density),
   stress(6),
   strain(6),
   Cstress(6),
   Cstrain(6),
   theTangent(tangent),
   CplastStrainDev(6),
   CbackStress(6),
   plastStrainDev(6),
   backStress(6)
{
  this->ndm     = pNd;
  this->G       = pG;
  this->K       = pK;
  this->sigmaY0 = pSigmaY0;
  this->H_kin   = pH_kin;
  this->H_iso   = pH_iso;

  stress.Zero();
  strain.Zero();
  sigmaY = pSigmaY0;

  Cstress.Zero();
  Cstrain.Zero();
  CsigmaY = pSigmaY0;

  lambda = 0.;
}

SimplifiedJ2::SimplifiedJ2()
 : NDMaterial(0, ND_TAG_SimplifiedJ2),
   stress(6),
   strain(6),
   Cstress(6),
   Cstrain(6),
   theTangent(tangent),
   CplastStrainDev(6),
   CbackStress(6),
   plastStrainDev(6),
   backStress(6)
{
}

SimplifiedJ2::SimplifiedJ2(const SimplifiedJ2& a)
 : NDMaterial(a.getTag(), ND_TAG_SimplifiedJ2),
   stress(6),
   strain(6),
   Cstress(6),
   Cstrain(6),
   theTangent(tangent),
   CplastStrainDev(6),
   CbackStress(6),
   plastStrainDev(6),
   backStress(6)
{
  this->ndm     = a.ndm;
  this->G       = a.G;
  this->K       = a.K;
  this->sigmaY0 = a.sigmaY0;
  this->H_kin   = a.H_kin;
  this->H_iso   = a.H_iso;

  stress.Zero();
  strain.Zero();
  sigmaY = a.sigmaY0;
  lambda = 0.;
  Cstress.Zero();
  Cstrain.Zero();
  CsigmaY = a.sigmaY0;
}


SimplifiedJ2::~SimplifiedJ2() { return; }


int
SimplifiedJ2::plastIntegrator()
{
  double trace = strain(0) + strain(1) + strain(2);


  VectorND<6> I2{1.0, 1.0, 0, 0, 0, 0};

  VectorND<6> strainDev;
  strainDev = strain;
  strainDev.addVector(1.0, I2, -trace / 3.0);

  //Vector TbackStress(6);
  //TbackStress = CbackStress;

  VectorND<6> TstressDev;
  TstressDev.addVector(0.0, strainDev, 2. * G);
  TstressDev.addVector(1.0, CplastStrainDev, -2. * G);

  VectorND<6> Teta = TstressDev;
  Teta.addVector(1.0, CbackStress, -1.0);


  // --- check elastic or plastic--
  double yieldFunction =
      pow((Teta && Teta), 0.5) - pow(2. / 3, 0.5) * CsigmaY; // to replace Yn=(2/3)^.5*sigmaYn

  if (yieldFunction > 0) { // plastic corrector

    lambda = yieldFunction / (2. * G + 2. / 3. * (H_iso + H_kin));

    sigmaY = CsigmaY +
             pow(2. / 3., 0.5) * H_iso * lambda; //  Note:to replace Y_n+1 = Yn=(2/3)^.5*sigmaYn

    VectorND<6> n;
    n.addVector(0, Teta, 1. / pow((Teta && Teta), 0.5));

    //Vector eta(6);
    //eta.addVector(0, n, pow( (Teta &&  Teta),0.5)-(2.*G+2./3.*H_kin)*lambda);

    backStress.addVector(0.0, CbackStress, 1.0);
    backStress.addVector(1.0, n, 2. / 3. * H_kin * lambda);

    plastStrainDev.addVector(0.0, CplastStrainDev, 1.0);
    plastStrainDev.addVector(1.0, n, lambda);

    // cumPlastStrainDev = CcumPlastStrainDev + pow(2./3.,0.5)*lambda;

    // sigmaY = CsigmaY + H_iso*pow(2./3., 0.5) * lambda;

    stress.addVector(0.0, TstressDev, 1.0);
    stress.addVector(1.0, n, -2. * G * lambda);
    //stress.addVector(1.0, I2, 1./3*K*trace);

    stress.addVector(1.0, I2, K * trace);

    //double qiu;
    //    qiu=lambda+CplastStrainDev(0);


    // -------- consistent tangent modulus ------------

    double A = 2. * G / (2. * G + 2. / 3. * H_kin + 2. / 3. * H_iso);

    double C = 2. * G * lambda / pow(Teta && Teta, 0.5);

    //  double D = 2./3.*H_kin*lambda/pow(Teta&&Teta, 0.5);


    //----------------//////////////////////////////////////////////////--------------------------------

    //    static constexpr MatrixND<6,6> I_dev;
    //    I_dev.Zero();

    //    for (int i=0; i<6; i++)
    //      I_dev(i,i) =1.0;

    //    for (int i=0; i<3; i++)
    //        for (int j=0; j<3; j++)
    //            I_dev(i,j) -= 1.0/3.0;

    //    static constexpr VectorND<6> I2{1.0, 1.0, 0, 0, 0, 0};  // unit vector order 2

    //    tmpMatrix.Zero();

    //    for(int i=0; i<3; i++)       // I2@I2
    //        for(int j=0; j<3; j++)
    //            tmpMatrix(i,j)=1.0;

    tangent.zero();
    tangent.addMatrix(IoI, K);

    tangent.addMatrix(IIdevMix, 2. * G * (1 - C));

    MatrixND<6, 6> non{};

    for (int i = 0; i < 6; i++) {
      for (int j = 0; j < 3; j++)
        non(i, j) = n[i] * n[j];
      for (int j = 3; j < 6; j++)
        non(i, j) =
            n[i] * n[j] *
            2.0; // To be consistent with the transformation between 4th order tensor and matrix
    }

    tangent.addMatrix(non, 2. * G * (C - A));
  }

  else {
    //
    // elastic
    //


    sigmaY = CsigmaY;
    backStress.addVector(0.0, CbackStress, 1.0);
    plastStrainDev.addVector(0.0, CplastStrainDev, 1.0);
    //cumPlastStrainDev = CcumPlastStrainDev;
    // sigmaY = CsigmaY;
    Vector n(6);
    n.addVector(0, Teta, 1. / pow((Teta && Teta), 0.5));

    //Vector eta(6);
    //eta.addVector(0, n, pow( (Teta &&  Teta),0.5));

    stress.addVector(0.0, TstressDev, 1.0);
    //stress.addVector(1.0, I2, 1./3*K*trace);
    stress.addVector(1.0, I2, K * trace);


    theTangent.Zero();

    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        theTangent(i, j) = K - 2.0 / 3.0 * G;

    for (int i = 0; i < 6; i++)
      theTangent(i, i) += 2.0 * G;
  }


  for (int i = 0; i < 6; i++)
    for (int j = 3; j < 6; j++)
      theTangent(i, j) /= 2.0;


  return 0;
}


int
SimplifiedJ2::setTrialStrain(const Vector& pStrain)
{

  if (ndm == 3 && pStrain.Size() == 6)
    strain = pStrain;
  else if (ndm == 2 && pStrain.Size() == 3) {
    strain[0] = pStrain[0];
    strain[1] = pStrain[1];
    strain[2] = 0.0;
    strain[3] = pStrain[2];
    strain[4] = 0.0;
    strain[5] = 0.0;
  } else {
    opserr << "Fatal:SimplifiedJ2:: Material dimension is: " << ndm << endln;
    opserr << "But strain vector size is: " << pStrain.Size() << endln;
    return -1;
  }

  // ----- change to real strain instead of eng. strain

  for (int i = 3; i < 6; i++)
    strain[i] /= 2.0;

  this->plastIntegrator();


  return 0;
}

int
SimplifiedJ2::setTrialStrain(const Vector& v, const Vector& r)
{
  return this->setTrialStrain(v);
}

int
SimplifiedJ2::setTrialStrainIncr(const Vector& v)
{
  // ----- change to real strain instead of eng. strain
  // ---- since all strain in material is the true strain, not eng.strain.

  for (int i = 0; i < 3; i++) {
    tmpVector(i)     = v(i);
    tmpVector(i + 3) = v(i + 3) / 2.0;
  }

  if (ndm == 3 && v.Size() == 6)
    strain = Cstrain + v;

  else if (ndm == 2 && v.Size() == 3) {
    strain[0] = Cstrain[0] + v[0];
    strain[1] = Cstrain[1] + v[1];
    strain[2] = 0.0;
    strain[3] = Cstrain[2] + v[2];
    strain[4] = 0.0;
    strain[5] = 0.0;
  } else {
    opserr << "Fatal:SimplifiedJ2:: Material dimension is: " << ndm << endln;
    opserr << "But strain vector size is: " << v.Size() << endln;
    return -1;
  }

  return this->plastIntegrator();
}


int
SimplifiedJ2::setTrialStrainIncr(const Vector& v, const Vector& r)
{
  return this->setTrialStrainIncr(v);
}


// Calculates current tangent stiffness.
const Matrix&
SimplifiedJ2::getTangent()
{
  if (ndm == 3)
    return theTangent;
  else {
    static Matrix workM(3, 3);
    workM(0, 0) = theTangent(0, 0);
    workM(0, 1) = theTangent(0, 1);
    workM(0, 2) = theTangent(0, 3);
    workM(1, 0) = theTangent(1, 0);
    workM(1, 1) = theTangent(1, 1);
    workM(1, 2) = theTangent(1, 3);
    workM(2, 0) = theTangent(3, 0);
    workM(2, 1) = theTangent(3, 1);
    workM(2, 2) = theTangent(3, 3);
    return workM;
  }
}

const Matrix&
SimplifiedJ2::getInitialTangent()
{
  // TODO: Initial tangent
  return this->getTangent();
}


const Vector&
SimplifiedJ2::getStress()
{
  return stress;
}

const Vector&
SimplifiedJ2::getStrain()
{
  return strain;
}

const Vector&
SimplifiedJ2::getCommittedStress()
{
  return Cstress;
}

const Vector&
SimplifiedJ2::getCommittedStrain()
{
  return Cstrain;
}


int
SimplifiedJ2::commitState()
{
  Cstress         = stress;
  Cstrain         = strain;
  CplastStrainDev = plastStrainDev;
  CbackStress     = backStress;
  CsigmaY         = sigmaY;
  //CcumPlastStrainDev = cumPlastStrainDev;

  return 0;
}

int
SimplifiedJ2::revertToLastCommit()
{
  // -- to be implemented.
  return 0;
}

int
SimplifiedJ2::revertToStart()
{
  // -- to be implemented.
  return 0;
}


NDMaterial*
SimplifiedJ2::getCopy()
{
  return new SimplifiedJ2(*this);
}


NDMaterial*
SimplifiedJ2::getCopy(const char* code)
{
  if (strcmp(code, "PlaneStress") == 0 || strcmp(code, "PlaneStrain") == 0 ||
      strcmp(code, "ThreeDimensional") == 0) {
    SimplifiedJ2* theJ2 = new SimplifiedJ2(*this);
    return theJ2;
  }
  return NDMaterial::getCopy(code);
}


int
SimplifiedJ2::sendSelf(int commitTag, Channel& theChannel)
{
  static Vector data(7 + 6 + 6 + 6 + 6 + 1);
  data(0) = this->getTag();
  data(1) = ndm;
  data(2) = G;
  data(3) = K;
  data(4) = sigmaY0;
  data(5) = H_kin;
  data(6) = H_iso;

  data(7) = CsigmaY;

  for (int i = 0; i < 6; i++) {
    data(8 + i)      = Cstress(i);
    data(8 + 6 + i)  = Cstrain(i);
    data(8 + 12 + i) = CplastStrainDev(i);
    data(8 + 18 + i) = CbackStress(i);
  }

  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "SimplifiedJ2::sendSelf - failed to send vector to channel\n";
    return -1;
  }

  return 0;
}

int
SimplifiedJ2::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
  static Vector data(7 + 6 + 6 + 6 + 6 + 1);
  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "SimplifiedJ2::recvSelf - failed to recv vector from channel\n";
    return -1;
  }

  this->setTag((int)data(0));
  ndm     = (int)data(1);
  G       = data(2);
  K       = data(3);
  sigmaY0 = data(4);
  H_kin   = data(5);
  H_iso   = data(6);

  CsigmaY = data(7);

  for (int i = 0; i < 6; i++) {
    Cstress(i)         = data(8 + i);
    Cstrain(i)         = data(8 + 6 + i);
    CplastStrainDev(i) = data(8 + 12 + i);
    CbackStress(i)     = data(8 + 18 + i);
  }

  return 0;
}


Response*
SimplifiedJ2::setResponse(const char** argv, int argc, OPS_Stream& s)
{

  if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "stresses") == 0)
    return new MaterialResponse(this, 1, stress);

  else if (strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "strains") == 0)
    return new MaterialResponse(this, 2, strain);

  else if (strcmp(argv[0], "tangent") == 0 || strcmp(argv[0], "Tangent") == 0)
    return new MaterialResponse(this, 3, theTangent);


  else if (strcmp(argv[0], "plasticStrainDev") == 0 || strcmp(argv[0], "plasticStrainDevs") == 0)
    return new MaterialResponse(this, 4, plastStrainDev);

  else
    return 0;
}


int
SimplifiedJ2::getResponse(int responseID, Information& matInfo)
{
  switch (responseID) {
  case -1: return -1;
  case 1:
    if (matInfo.theVector != 0)
      *(matInfo.theVector) = stress;
    return 0;

  case 2:
    if (matInfo.theVector != 0)
      *(matInfo.theVector) = strain;
    return 0;

  case 3:
    if (matInfo.theMatrix != 0)
      *(matInfo.theMatrix) = theTangent;
    return 0;

  case 4:
    if (matInfo.theVector != 0)
      *(matInfo.theVector) = plastStrainDev;
    return 0;
  }

  return 0;
}

void
SimplifiedJ2::Print(OPS_Stream& s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"" << this->getClassType() << "\", ";
    s << "\"G\": " << G << ", ";
    s << "\"K\": " << K << ", ";
    s << "\"Fy\": " << sigmaY0 << ", ";
    s << "\"Hiso\": " << H_iso << ", ";
    s << "\"Hkin\": " << H_kin << ", ";
    s << "\"density\": " << density;
    s << "}";
  }
  return;
}


int
SimplifiedJ2::setParameter(const char** argv, int argc, Parameter& param)
{
  // -- to be implemented.
  return 0;
}

int
SimplifiedJ2::updateParameter(int responseID, Information& eleInformation)
{
  // -- to be implemented.
  /*switch (passedParameterID) {
    case -1:
        return -1;

    case 1:
        this->Nd= info.theDouble; // 
        break;

    case 2:
        this->G  = info.theDouble;
        break;

    case 3:
        this->K = info.theDouble; // 
        break;

    case 4:
        this->SigmaY0  = info.theDouble;
        break;


    case 5:
        this->H_kin= info.theDouble; // 
        break;

    case 6:
        this->H_iso= info.theDouble; // 
        break;
    }
    */
  return 0;
}
