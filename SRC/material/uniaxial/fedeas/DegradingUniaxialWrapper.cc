#include <stdlib.h>
#include <string.h>
#include <unordered_map>
#include <functional>
#include <tcl.h>

#include <DegradingUniaxialWrapper.hh>
#include <ID.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

#include <OPS_Globals.h>

#include <elementAPI.h>

// TODO
#define MAT_TAG_FedeasDeg 90909090

DegradingUniaxialWrapper::DegradingUniaxialWrapper(int tag,
                                                   UniaxialMaterial &material,
                                                   double min, double max)
    : UniaxialMaterial(tag, MAT_TAG_FedeasDeg), theMaterial(0), minStrain(min),
      maxStrain(max), m_stress(0.0)
{
  theMaterial = material.getCopy();
  m_tangent = theMaterial->getInitialTangent();

  if (theMaterial == 0) {
    opserr << "DegradingUniaxialWrapper::DegradingUniaxialWrapper -- failed to "
              "get copy of material\n";
    exit(-1);
  }
}

DegradingUniaxialWrapper::DegradingUniaxialWrapper()
    : UniaxialMaterial(0, MAT_TAG_FedeasDeg), theMaterial(0), minStrain(0.0),
      maxStrain(0.0)
{
}

DegradingUniaxialWrapper::~DegradingUniaxialWrapper()
{
  if (theMaterial)
    delete theMaterial;
}

int
DegradingUniaxialWrapper::setCoupling(double coupling_param){return 0;}

int 
DegradingUniaxialWrapper::setDamageWrapper(Tcl_Interp *interp, std::string tag) {

    typedef std::function<int(void*,void*)> degrade_t;
    // std::string tag = Tcl_GetString(objv[1]);
    
    auto map = (std::unordered_map<std::string, degrade_t(*)(Tcl_Interp*, std::string)>*)
      Tcl_GetAssocData(interp, "elle::libdmg::DamageGenerators", NULL);
    if (!map){
      printf("Failed to get constructor map");
      return -1;
    }

    this->degrade = (degrade_t) ((*map)["Uniaxial"](interp, tag));

    return 1;
}
  
int
DegradingUniaxialWrapper::setTrialStrain(double strain, double temp, double strainRate)
{

    double pastStrain = theMaterial->getStrain();
    theMaterial->setTrialStrain(strain, temp, strainRate);
    double trialStrain = theMaterial->getStrain();

    double strain_incr = trialStrain - pastStrain;

    if (degrade) {//  && abs(strain_incr) > m_rate_tol){

      UniaxialState pres, past = {
        .e   = pastStrain,
        .ep  = trialStrain,
        .De  = strain_incr,
        .se  = theMaterial->getStress(),
        .kt  = theMaterial->getTangent(),
        .ke  = theMaterial->getInitialTangent()
      };

      this->degrade((void*)&past, (void*)&pres);

      this->m_stress   = pres.se;
      this->m_tangent  = pres.kt;

    } else {
      this->m_stress   = theMaterial->getStress();
      this->m_tangent  = theMaterial->getTangent();
    }
    return 0;
}

int
DegradingUniaxialWrapper::setTrialStrain(double trialStrain, double strainRate)
{

  return this->setTrialStrain(trialStrain, 0.0, strainRate);
  /*
    double pastStrain = theMaterial->getStrain();
    theMaterial->setTrialStrain(trialStrain, strainRate);
    // double trialStrain = theMaterial->getStrain();
    if (degrade){
      // double currStrain = BaseMaterial::getStrain();
      // opserr << "#Degrading\n";

      UniaxialState pres, past = {
        .e   = pastStrain,
        .ep  = trialStrain,
        .De  = trialStrain - pastStrain,
        .se  = theMaterial->getStress(),
        .kt  = theMaterial->getTangent(),
        .ke  = theMaterial->getInitialTangent()
      };

      degrade((void*)&past, (void*)&pres);

      this->m_stress  = pres.se;
      this->m_tangent = pres.kt;
    } else {
      // opserr << "#Not Degrading\n";
    }
    return 0;
*/
}

double
DegradingUniaxialWrapper::getTangent(void)
{
  // return theMaterial->getTangent();
  
  if (degrade)
    return this->m_tangent;
  else
    return theMaterial->getTangent();
  
}

double
DegradingUniaxialWrapper::getStress(void)
{
  if (degrade)
    return m_stress;
  else
    return theMaterial->getStress();
}

double
DegradingUniaxialWrapper::getInitialTangent(void){return theMaterial->getInitialTangent();}

double
DegradingUniaxialWrapper::getDampTangent(void){return theMaterial->getDampTangent();}

double
DegradingUniaxialWrapper::getStrain(void)
{return theMaterial->getStrain();}

double
DegradingUniaxialWrapper::getStrainRate(void)
{return theMaterial->getStrainRate();}

int
DegradingUniaxialWrapper::commitState(void)
{return theMaterial->commitState();}

int
DegradingUniaxialWrapper::revertToLastCommit(void)
{return theMaterial->revertToLastCommit();}

int
DegradingUniaxialWrapper::revertToStart(void)
{return theMaterial->revertToStart();}

UniaxialMaterial *
DegradingUniaxialWrapper::getCopy(void)
{
  DegradingUniaxialWrapper *theCopy = new DegradingUniaxialWrapper(
      this->getTag(), *theMaterial, minStrain, maxStrain);
  theCopy->degrade = this->degrade;

  return theCopy;
}

int
DegradingUniaxialWrapper::sendSelf(int cTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();

  static ID dataID(3);
  dataID(0) = this->getTag();
  dataID(1) = theMaterial->getClassTag();
  int matDbTag = theMaterial->getDbTag();
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    theMaterial->setDbTag(matDbTag);
  }
  dataID(2) = matDbTag;
  if (theChannel.sendID(dbTag, cTag, dataID) < 0) {
    opserr << "DegradingUniaxialWrapper::sendSelf() - failed to send the ID\n";
    return -1;
  }

  static Vector dataVec(3);
  dataVec(0) = minStrain;
  dataVec(1) = maxStrain;
  dataVec(2) = 0.0;

  if (theChannel.sendVector(dbTag, cTag, dataVec) < 0) {
    opserr
        << "DegradingUniaxialWrapper::sendSelf() - failed to send the Vector\n";
    return -2;
  }

  if (theMaterial->sendSelf(cTag, theChannel) < 0) {
    opserr << "DegradingUniaxialWrapper::sendSelf() - failed to send the "
              "Material\n";
    return -3;
  }

  return 0;
}

int
DegradingUniaxialWrapper::recvSelf(int cTag, Channel &theChannel,
                                   FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static ID dataID(3);
  if (theChannel.recvID(dbTag, cTag, dataID) < 0) {
    opserr << "DegradingUniaxialWrapper::recvSelf() - failed to get the ID\n";
    return -1;
  }
  this->setTag(int(dataID(0)));

  // as no way to change material, don't have to check classTag of the material
  if (theMaterial == 0) {
    int matClassTag = int(dataID(1));
    theMaterial = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theMaterial == 0) {
      opserr << "DegradingUniaxialWrapper::recvSelf() - failed to create "
                "Material with classTag "
             << dataID(0) << endln;
      return -2;
    }
  }
  theMaterial->setDbTag(dataID(2));

  static Vector dataVec(3);
  if (theChannel.recvVector(dbTag, cTag, dataVec) < 0) {
    opserr
        << "DegradingUniaxialWrapper::recvSelf() - failed to get the Vector\n";
    return -3;
  }

  minStrain = dataVec(0);
  maxStrain = dataVec(1);



  if (theMaterial->recvSelf(cTag, theChannel, theBroker) < 0) {
    opserr << "DegradingUniaxialWrapper::recvSelf() - failed to get the "
              "Material\n";
    return -4;
  }
  return 0;
}

void
DegradingUniaxialWrapper::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
    s << "DegradingUniaxialWrapper, tag: " << this->getTag() << endln;
    s << "  material: " << theMaterial->getTag() << endln;
    s << "  min strain: " << minStrain << endln;
    s << "  max strain: " << maxStrain << endln;
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"DegradingUniaxialWrapper\", ";
    s << "\"material\": \"" << theMaterial->getTag() << "\", ";
    s << "\"epsMin\": " << minStrain << ", ";
    s << "\"epsMax\": " << maxStrain << "}";
  }
}

int
DegradingUniaxialWrapper::setParameter(const char **argv, int argc,
                                       Parameter &param)
{
  //
  // I suppose epsMin and epsMax could be parameters, but not for now -- MHS
  //
  return theMaterial->setParameter(argv, argc, param);
}

int
DegradingUniaxialWrapper::updateParameter(int parameterID, Information &info)
{
  return 0;
}

double
DegradingUniaxialWrapper::getStressSensitivity(int gradIndex, bool conditional)
{
    return theMaterial->getStressSensitivity(gradIndex, conditional);
}

double
DegradingUniaxialWrapper::getStrainSensitivity(int gradIndex)
{
  return theMaterial->getStrainSensitivity(gradIndex);
}

double
DegradingUniaxialWrapper::getInitialTangentSensitivity(int gradIndex)
{
  return theMaterial->getInitialTangentSensitivity(gradIndex);
}

double
DegradingUniaxialWrapper::getDampTangentSensitivity(int gradIndex)
{
  return theMaterial->getDampTangentSensitivity(gradIndex);
}

double
DegradingUniaxialWrapper::getRhoSensitivity(int gradIndex)
{
  return theMaterial->getRhoSensitivity(gradIndex);
}

int
DegradingUniaxialWrapper::commitSensitivity(double strainGradient,
                                            int gradIndex, int numGrads)
{
  return theMaterial->commitSensitivity(strainGradient, gradIndex, numGrads);
}
