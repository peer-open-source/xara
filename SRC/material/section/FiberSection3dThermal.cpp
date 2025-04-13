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
// Written: fmk
// Created: 04/04
//
// Description: This file contains the class implementation of FiberSection2d.
// Modified for SIF modelling by Jian Jiang,Liming Jiang [http://openseesforfire.github.io]


#include <stdlib.h>

#include <Channel.h>
#include <Vector.h>
#include <Matrix.h>
#include <classTags.h>
#include <Logging.h>
#include <OPS_Stream.h>
#include <FiberSection3dThermal.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <MaterialResponse.h>
#include <UniaxialMaterial.h>
#include <ElasticMaterial.h>
#include <math.h>

namespace OpenSees {
ID FiberSection3dThermal::code(4);

#if 0
#include <elementAPI.h>
void* OPS_FiberSection3dThermal()
{
    int numData = OPS_GetNumRemainingInputArgs();
    if(numData < 1) {
            opserr<<"insufficient arguments for FiberSection3d\n";
            return 0;
    }
    
    numData = 1;
    int tag;
    if (OPS_GetIntInput(&numData, &tag) < 0) return 0;

    if (OPS_GetNumRemainingInputArgs() < 2) {
      opserr << "WARNING torsion not specified for FiberSection\n";
      opserr << "Use either -GJ $GJ or -torsion $matTag\n";
      opserr << "\nFiberSection3d section: " << tag << "\n";
      return 0;
    }
    
    UniaxialMaterial *torsion = 0;
    bool deleteTorsion = false;
    bool computeCentroid = true;
    while (OPS_GetNumRemainingInputArgs() > 0) {
      const char* opt = OPS_GetString();
      if (strcmp(opt,"-noCentroid") == 0) {
        computeCentroid = false;
      }
      if (strcmp(opt, "-GJ") == 0 && OPS_GetNumRemainingInputArgs() > 0) {
        numData = 1;
        double GJ;
        if (OPS_GetDoubleInput(&numData, &GJ) < 0) {
          opserr << "WARNING: failed to read GJ\n";
          return 0;
        }
        torsion = new ElasticMaterial(0,GJ);
        deleteTorsion = true;
      }
      if (strcmp(opt, "-torsion") == 0 && OPS_GetNumRemainingInputArgs() > 0) {
        numData = 1;
        int torsionTag;
        if (OPS_GetIntInput(&numData, &torsionTag) < 0) {
          opserr << "WARNING: failed to read torsion\n";
          return 0;
        }
        torsion = OPS_getUniaxialMaterial(torsionTag);
      }
    }

    if (torsion == 0) {
      opserr << "WARNING torsion not specified for FiberSection\n";
      opserr << "\nFiberSection3d section: " << tag << "\n";
      return 0;
    }
    
    int num = 30;
    FrameSection *section = new FiberSection3dThermal(tag, num, *torsion, computeCentroid);
    if (deleteTorsion)
      delete torsion;
    return section;
}
#endif

#if 0
// constructors:
FiberSection3dThermal::FiberSection3dThermal(int tag, int num, Fiber **fibers,
                                             UniaxialMaterial &torsion,  bool compCentroid):
  FrameSection(tag, SEC_TAG_FiberSection3dThermal),
  sizeFibers(num), theMaterials(0),
  QzBar(0.0), QyBar(0.0), ABar(0.0), yBar(0.0), zBar(0.0), computeCentroid(compCentroid),
  e(4), eCommit(4), s(0), ks(0), theTorsion(0), sT(3), Fiber_T(0), Fiber_TMax(0),
  parameterID(0), SHVs(0), AverageThermalElong(4)
{
  if (fibers.size() > 0) {
    
    for (int i = 0; i < fibers.size(); i++) {
      Fiber *theFiber = fibers[i];
      double yLoc, zLoc, Area;
      theFiber->getFiberLocation(yLoc, zLoc);
      Area = theFiber->getArea();
      QzBar += yLoc*Area;
      QyBar += zLoc*Area;
      ABar  += Area;

      fibers[i].y = -yLoc;
      fibers[i].z = zLoc;
      fibers[i].area = Area;
      UniaxialMaterial *theMat = theFiber->getMaterial();
      fibers[i].material = theMat->getCopy();

      if (fibers[i].material == 0) {
        opserr << "FiberSection3dThermal::FiberSection3dThermal -- failed to get copy of a Material\n";
        exit(-1);
      }

      Fiber_T[i] = 0.0;
      Fiber_TMax[i] = 0.0;
    }

    if (computeCentroid) {
      yBar = QzBar/ABar;
      zBar = QyBar/ABar;
    }
  }

  theTorsion = torsion.getCopy();
  if (theTorsion == 0)
    opserr << "FiberSection3d::FiberSection3d -- failed to get copy of torsion material\n";

  s = new Vector(sData, 3);
  ks = new Matrix(kData, 3, 3);

  sData[0] = 0.0;
  sData[1] = 0.0;
  sData[2] = 0.0;

  for (int i=0; i<16; i++)
    kData[i] = 0.0;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_MY;
  code(3) = SECTION_RESPONSE_T;

 // AddingSensitivity:BEGIN ////////////////////////////////////
  parameterID = 0;
  SHVs=0;
  // AddingSensitivity:END //////////////////////////////////////
}
#endif

FiberSection3dThermal::FiberSection3dThermal(int tag, int num, UniaxialMaterial &torsion, bool compCentroid):
  FrameSection(tag, SEC_TAG_FiberSection3dThermal),
  sizeFibers(num), 
  QzBar(0.0), QyBar(0.0), ABar(0.0), yBar(0.0), zBar(0.0), 
  computeCentroid(compCentroid),
  e(4), eCommit(4), s(0), ks(0), theTorsion(nullptr),
  sT(3),
  parameterID(0), SHVs(0)
{
    fibers.reserve(num);

    theTorsion = torsion.getCopy();
    
    s = new Vector(sData, 4);
    ks = new Matrix(kData, 4, 4);

    sData[0] = 0.0;
    sData[1] = 0.0;
    sData[2] = 0.0;

    for (int i=0; i<16; i++)
      kData[i] = 0.0;

    code(0) = SECTION_RESPONSE_P;
    code(1) = SECTION_RESPONSE_MZ;
    code(2) = SECTION_RESPONSE_MY;
    code(3) = SECTION_RESPONSE_T;

   // AddingSensitivity:BEGIN ////////////////////////////////////
    parameterID = 0;
    SHVs=0;
    // AddingSensitivity:END //////////////////////////////////////
}

// constructor for blank object that recvSelf needs to be invoked upon
FiberSection3dThermal::FiberSection3dThermal():
  FrameSection(0, SEC_TAG_FiberSection3dThermal),
  sizeFibers(0),
  QzBar(0.0), QyBar(0.0), ABar(0.0), yBar(0.0), zBar(0.0), computeCentroid(true),
  e(4), eCommit(4), s(0), ks(0), theTorsion(nullptr),
  sT(3),
  parameterID(0), SHVs(0)
{
  s = new Vector(sData, 4);
  ks = new Matrix(kData, 4, 4);

  sData[0] = 0.0;
  sData[1] = 0.0;
  sData[2] = 0.0;

  for (int i=0; i<16; i++)
    kData[i] = 0.0;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_MY;
  code(3) = SECTION_RESPONSE_T;

 // AddingSensitivity:BEGIN ////////////////////////////////////
  parameterID = 0;
  SHVs=0;
  // AddingSensitivity:END //////////////////////////////////////
}

int
FiberSection3dThermal::addFiber(UniaxialMaterial &theMat, 
                                double area, 
                                double yLoc, 
                                double zLoc)
{

  fibers.push_back({
      .area=area,
      .y=yLoc,
      .z=zLoc,
      .temp=0.0,
      .temp_max=0.0,
      .material=theMat.getCopy(),
  });


  // Recompute centroid
  if (computeCentroid) {
    ABar  += area;
    QzBar += yLoc*area;
    QyBar += zLoc*area;
    
    yBar = QzBar/ABar;
    zBar = QyBar/ABar;
  }
  
  return 0;
}



// destructor:
FiberSection3dThermal::~FiberSection3dThermal()
{
  for (int i = 0; i < fibers.size(); i++)
    if (fibers[i].material != nullptr)
      delete fibers[i].material;

  if (s != 0)
    delete s;

  if (ks != 0)
    delete ks;

  //if (TemperatureTangent != 0)
    //delete [] TemperatureTangent;

  if (theTorsion != 0)
    delete theTorsion;  
}

int
FiberSection3dThermal::setTrialSectionDeformation (const Vector &deforms)
{
  int res = 0;
  e = deforms;

  for (int i = 0; i < 4; i++)
      sData[i] = 0.0;
  for (int i = 0; i < 16; i++)
      kData[i] = 0.0;

  int loc = 0;

  double d0 = deforms(0);
  double d1 = deforms(1);
  double d2 = deforms(2);
  double d3 = deforms(3);

  for (FiberData& fiber : fibers) {
      UniaxialMaterial *theMat = fiber.material;
      double y = fiber.y - yBar;
      double z = fiber.z - zBar;
      double A = fiber.area;

      double FiberTemperature = fiber.temp; // Added by Liming to obtain fiber T;
      double FiberTempMax= fiber.temp_max;  // Maximum Temp;

      //---Calculating the Fiber Temperature---end

      double strain = d0 + y*d1 + z*d2;  //axial strain d0, rotational degree d1,d2;
      double tangent =0.0;
      double stress = 0.0;
      double ThermalElongation = 0.0;
      static Vector tData(4);
      static Information iData(tData);
      tData(0) = FiberTemperature;
      tData(1) = tangent;
      tData(2) = ThermalElongation;
      tData(3) = FiberTempMax;
      iData.setVector(tData);
      theMat->getVariable("ElongTangent", iData);
      tData = iData.getData();
      tangent = tData(1);
      ThermalElongation = tData(2);

      // determine material strain and set it
      strain = d0 + y*d1 + z*d2 - ThermalElongation;
      res += theMat->setTrial(strain, FiberTemperature, stress, tangent, ThermalElongation);

      double EA = tangent * A;
      double vas1 = y*EA;
      double vas2 = z*EA;
      double vas1as2 = vas1*z;

      kData[0] += EA;
      kData[1] += vas1;
      kData[2] += vas2;

      kData[5] += vas1 * y;
      kData[6] += vas1as2;

      kData[10] += vas2 * z;

      double fs0 = stress * A;

      sData[0] += fs0;
      sData[1] += fs0 * y;
      sData[2] += fs0 * z;
  }

  kData[4] = kData[1];
  kData[8] = kData[2];
  kData[9] = kData[6];

  if (theTorsion != 0) {
    double stress, tangent;
    res += theTorsion->setTrial(d3, stress, tangent);
    sData[3] = stress;
    kData[15] = tangent;
  }

  return res;
}

const Matrix&
FiberSection3dThermal::getInitialTangent()
{
  static double kInitialData[16];
  static Matrix kInitial(kInitialData, 4, 4);
  
  kInitial.Zero();

  int loc = 0;

  for (int i = 0; i < fibers.size(); i++) {
    UniaxialMaterial *theMat = fibers[i].material;
    double y = fibers[i].y - yBar;
    double z = fibers[i].z - zBar;
    double A = fibers[i].area;

    double tangent = theMat->getInitialTangent();

    double value = tangent * A;
    double vas1 = y*value;
    double vas2 = z*value;
    double vas1as2 = vas1*z;

    kInitialData[0] += value;
    kInitialData[1] += vas1;
    kInitialData[2] += vas2;

    kInitialData[5] += vas1 * y;
    kInitialData[6] += vas1as2;

    kInitialData[10] += vas2 * z;
  }

  kInitialData[4] = kInitialData[1];
  kInitialData[8] = kInitialData[2];
  kInitialData[9] = kInitialData[6];

  if (theTorsion != 0)
    kInitialData[15] = theTorsion->getInitialTangent();

  return kInitial;
}

const Vector&
FiberSection3dThermal::getSectionDeformation()
{
  return e;
}

const Matrix&
FiberSection3dThermal::getSectionTangent()
{
  return *ks;
}

const Vector&
FiberSection3dThermal::getStressResultant()
{
  return *s;
}


//JJadd--12.2010---to get section force due to thermal load----start-----
const Vector&
FiberSection3dThermal::getTemperatureStress(const Vector& dataMixed)
{
  sT.Zero();
  AverageThermalElong.zero();
  //JJadd, 12/2010, updata yBar = Ai*Ei*yi/(Ai*E*)  start
  double ThermalTangent[1000];
  double ThermalElong[1000];
  for (int i = 0; i < fibers.size(); i++) {
        ThermalTangent[i] = 0.0;
        ThermalElong[i] = 0.0;
  }

  for (int i=0; i< fibers.size(); i++) {

    FiberData& fiber = fibers[i];
    UniaxialMaterial *theMat = fiber.material;

    double FiberTemperature= this->determineFiberTemperature( dataMixed, -fiber.y, fiber.z);

    // determine material strain and set it
    double tangent =0.0;
    double ThermalElongation =0.0;
    static Vector tData(4);
    static Information iData(tData);
    tData(0) = FiberTemperature;
    tData(1) = tangent;
    tData(2) = ThermalElongation;
    tData(3) = 0.0;
    iData.setVector(tData);
    theMat->getVariable("ElongTangent", iData);
    tData = iData.getData();
    FiberTemperature = tData(0);
    tangent = tData(1);
    ThermalElongation = tData(2);
    double FiberTempMax = tData(3);

    //  double strain = -ThermalElongation;
    //  theMat->setTrialTemperature(strain, FiberTemperature, stress, tangent, ThermalElongation);
    fiber.temp     = FiberTemperature;
    fiber.temp_max = FiberTempMax;
    ThermalTangent[i]  = tangent;
    ThermalElong[i]    = ThermalElongation;

  }

 // calculate section resisting force due to thermal load

  double FiberForce;
  double SectionArea = 0;
  double ThermalForce = 0;
  double ThermalMomentY = 0; double ThermalMomentZ = 0;
  double SectionMomofAreaY = 0; double SectionMomofAreaZ = 0;

  for (int i = 0; i < fibers.size(); i++) {
      FiberForce = ThermalTangent[i]*fibers[i].area*ThermalElong[i];
      sT(0) += FiberForce;
      sT(1) += FiberForce*(fibers[i].y - yBar);
      sT(2) += FiberForce*(fibers[i].z - zBar);
      // added GR
      SectionArea += fibers[i].area;
      SectionMomofAreaY += (fibers[i].area * (fibers[i].y - yBar) * (fibers[i].y - yBar));
      SectionMomofAreaZ += (fibers[i].area * (fibers[i].z - zBar) * (fibers[i].z - zBar));
      ThermalForce += ThermalElong[i] * fibers[i].area;
      ThermalMomentY += ThermalElong[i] * fibers[i].area * (fibers[i].y - yBar);
      ThermalMomentZ += ThermalElong[i] * fibers[i].area * (fibers[i].z - zBar);
  }
  //double ThermalMoment;
  //ThermalMoment = abs(sTData[1]);
 // sTData[1] = ThermalMoment;
  AverageThermalElong(0) = ThermalForce / SectionArea;
  AverageThermalElong(1) = ThermalMomentY / SectionMomofAreaY;
  AverageThermalElong(2) = ThermalMomentZ / SectionMomofAreaZ;
  AverageThermalElong(3) = 0.0; // no contribution in torsion

  return sT;
}
//JJadd--12.2010---to get section force due to thermal load----end-----

//UoE group///Calculating Thermal stresses at each /////////////////////////////////////////////////////end
const Vector&
FiberSection3dThermal::getThermalElong()
{
    static Vector wrapper(4);
    wrapper.setData(AverageThermalElong);
    return wrapper;
}

//Retuning ThermalElongation

FrameSection*
FiberSection3dThermal::getFrameCopy()
{
  FiberSection3dThermal *theCopy = new FiberSection3dThermal();
  theCopy->setTag(this->getTag());

  if (fibers.size() > 0) {

    theCopy->fibers = fibers;

    for (int i = 0; i < fibers.size(); i++) {
      theCopy->fibers[i].y = fibers[i].y;
      theCopy->fibers[i].z = fibers[i].z;
      theCopy->fibers[i].area = fibers[i].area;
      theCopy->fibers[i].material = fibers[i].material->getCopy();
      theCopy->fibers[i].temp = fibers[i].temp;
      theCopy->fibers[i].temp_max = fibers[i].temp_max;
    }
  }

  theCopy->eCommit = eCommit;
  theCopy->e = e;
  theCopy->QzBar = QzBar;
  theCopy->QyBar = QyBar;
  theCopy->ABar = ABar;  
  theCopy->yBar = yBar;
  theCopy->zBar = zBar;
  theCopy->computeCentroid = computeCentroid;
  
  for (int i=0; i<16; i++)
    theCopy->kData[i] = kData[i];

  theCopy->sData[0] = sData[0];
  theCopy->sData[1] = sData[1];
  theCopy->sData[2] = sData[2];
  theCopy->sData[3] = sData[3];
  theCopy->sT = sT;
  
  if (theTorsion != 0)
    theCopy->theTorsion = theTorsion->getCopy();
  else
    theCopy->theTorsion = 0;

  return theCopy;
}

const ID&
FiberSection3dThermal::getType ()
{
  return code;
}

int
FiberSection3dThermal::getOrder () const
{
  return 4;
}

int
FiberSection3dThermal::commitState()
{
  int err = 0;

  for (int i = 0; i < fibers.size(); i++)
    err += fibers[i].material->commitState();

  if (theTorsion != 0)
    err += theTorsion->commitState();
  
  eCommit = e;

  return err;
}

int
FiberSection3dThermal::revertToLastCommit()
{
  int err = 0;

  // Last committed section deformations
  e = eCommit;


  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0; kData[3] = 0.0;
  kData[4] = 0.0; kData[5] = 0.0; kData[6] = 0.0; kData[7] = 0.0;
  kData[8] = 0.0;
  kData[15] = 0.0;
  sData[0] = 0.0; sData[1] = 0.0;  sData[2] = 0.0; sData[3] = 0.0;

  int loc = 0;

  for (int i = 0; i < fibers.size(); i++) {
    UniaxialMaterial *theMat = fibers[i].material;
    double y = fibers[i].y - yBar;
    double z = fibers[i].z - zBar;
    double A = fibers[i].area;

    // invoke revertToLast on the material
    err += theMat->revertToLastCommit();

    double tangent = theMat->getTangent();
    double stress = theMat->getStress();

    double value = tangent * A;
    double vas1 = y*value;
    double vas2 = z*value;
    double vas1as2 = vas1*z;

    kData[0] += value;
    kData[1] += vas1;
    kData[2] += vas2;

    kData[5] += vas1 * y;
    kData[6] += vas1as2;

    kData[10] += vas2 * z;

    double fs0 = stress * A;
    sData[0] += fs0;
    sData[1] += fs0 * y;
    sData[2] += fs0 * z;
  }

  kData[4] = kData[1];
  kData[8] = kData[2];
  kData[9] = kData[6];

  if (theTorsion != 0) {
    err += theTorsion->revertToLastCommit();
    kData[15] = theTorsion->getTangent();
  } else
    kData[15] = 0.0;

  return err;
}

int
FiberSection3dThermal::revertToStart()
{
  // revert the fibers to start
  int err = 0;


  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0; kData[3] = 0.0;
  kData[4] = 0.0; kData[5] = 0.0; kData[6] = 0.0; kData[7] = 0.0;
  kData[8] = 0.0; kData[15] = 0.0;
  sData[0] = 0.0; sData[1] = 0.0;  sData[2] = 0.0; sData[3] = 0.0;

  int loc = 0;

  for (int i = 0; i < fibers.size(); i++) {
    UniaxialMaterial *theMat = fibers[i].material;
    double y = fibers[i].y - yBar;
    double z = fibers[i].z - zBar;
    double A = fibers[i].area;

    // invoke revertToStart on the material
    err += theMat->revertToStart();

    double tangent = theMat->getTangent();
    double stress = theMat->getStress();

    double value = tangent * A;
    double vas1 = y*value;
    double vas2 = z*value;
    double vas1as2 = vas1*z;

    kData[0] += value;
    kData[1] += vas1;
    kData[2] += vas2;

    kData[5] += vas1 * y;
    kData[6] += vas1as2;

    kData[10] += vas2 * z;

    double fs0 = stress * A;
    sData[0] += fs0;
    sData[1] += fs0 * y;
    sData[2] += fs0 * z;
  }

  kData[4] = kData[1];
  kData[8] = kData[2];
  kData[9] = kData[6];

  if (theTorsion != 0) {
    err += theTorsion->revertToStart();
    kData[15] = theTorsion->getTangent();
    sData[3] = theTorsion->getStress();
  } else {
    kData[15] = 0.0;
    sData[3] = 0.0;
  }

  return err;
}

int
FiberSection3dThermal::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;

  // create an id to send objects tag and fibers.size(),
  //     size 3 so no conflict with matData below if just 1 fiber
  static ID data(9);
  data(0) = this->getTag();
  data(1) = fibers.size();
  //data(2) = computeCentroid ? 1 : 0; // Now the ID data is really 3
  data(2) = (theTorsion != 0) ? 1 : 0;
  if (theTorsion != 0) {
    data(3) = theTorsion->getClassTag();
    int torsionDbTag = theTorsion->getDbTag();
    if (torsionDbTag == 0) {
      torsionDbTag = theChannel.getDbTag();
      if (torsionDbTag != 0)
        theTorsion->setDbTag(torsionDbTag);
    }
    data(4) = torsionDbTag;
  }
  data(5) = computeCentroid ? 1 : 0; // Now the ID data is really 5


  int dbTag = this->getDbTag();  
  res += theChannel.sendID(dbTag, commitTag, data);
  if (res < 0) {
    opserr << "FiberSection3d::sendSelf - failed to send ID data\n";
    return res;
  }    

  if (theTorsion != 0)
    theTorsion->sendSelf(commitTag, theChannel);

  
  if (fibers.size() != 0) {
    
    // create an id containingg classTag and dbTag for each material & send it
    ID materialData(2*fibers.size());
    for (int i=0; i<fibers.size(); i++) {
      UniaxialMaterial *theMat = fibers[i].material;
      materialData(2*i) = theMat->getClassTag();
      int matDbTag = theMat->getDbTag();
      if (matDbTag == 0) {
        matDbTag = theChannel.getDbTag();
        if (matDbTag != 0)
          theMat->setDbTag(matDbTag);
      }
      materialData(2*i+1) = matDbTag;
    }    
    
    res += theChannel.sendID(dbTag, commitTag, materialData);
    if (res < 0) {
     opserr << "FiberSection3d::sendSelf - failed to send material data\n";
     return res;
    }    

    // send the fiber data, i.e. area and loc, T, and Tmax
    Vector fiberData(5*fibers.size());
    for (int i = 0; i < fibers.size(); i++) {
      fiberData(0 + i*5) = fibers[i].y;
      fiberData(1 + i*5) = fibers[i].z;
      fiberData(2 + i*5) = fibers[i].area;
      fiberData(3 + i*5) = fibers[i].temp;
      fiberData(4 + i*5) = fibers[i].temp_max;
    }    
    res += theChannel.sendVector(dbTag, commitTag, fiberData);
    if (res < 0) {
     opserr << "FiberSection2d::sendSelf - failed to send material data\n";
     return res;
    }

    // now invoke send(0 on all the materials
    for (int j=0; j<fibers.size(); j++) {
      fibers[j].material->sendSelf(commitTag, theChannel);
      if (res < 0) {
        opserr << "FiberSection3d::sendSelf - failed to send material with tag "
               << fibers[j].material->getTag() << "\n";
        return res;
      }
    }
  }

  return res;
}

int
FiberSection3dThermal::recvSelf(int commitTag, Channel &theChannel,
                         FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static ID data(9);
  
  int dbTag = this->getDbTag();
  res += theChannel.recvID(dbTag, commitTag, data);
  if (res < 0) {
   opserr << "FiberSection3d::recvSelf - failed to recv ID data\n";
   return res;
  } 
  this->setTag(data(0));

  if (data(2) == 1 && theTorsion == 0) {        
    int torsionClassTag = data(3);
    int torsionDbTag = data(4);
    theTorsion = theBroker.getNewUniaxialMaterial(torsionClassTag);
    if (theTorsion == 0) {
      opserr << "FiberSection3d::recvSelf - failed to get torsion material \n";
      return -1;
    }
    theTorsion->setDbTag(torsionDbTag);
  }

  if (theTorsion->recvSelf(commitTag, theChannel, theBroker) < 0) {
           opserr << "FiberSection3d::recvSelf - torsion failed to recvSelf \n";
       return -2;
  }

  
  // recv data about materials objects, classTag and dbTag
  if (data(1) != 0) {
    ID materialData(2*data(1));
    res += theChannel.recvID(dbTag, commitTag, materialData);
    if (res < 0) {
     opserr << "FiberSection3d::recvSelf - failed to recv material data\n";
     return res;
    }    

    // if current arrays not of correct size, release old and resize
    if (fibers.size() != data(1)) {
      // delete old stuff if outa date
      fibers.clear();
    }

    Vector fiberData(5*fibers.size());
    res += theChannel.recvVector(dbTag, commitTag, fiberData);
    if (res < 0) {
       opserr << "FiberSection3d::recvSelf - failed to recv fiber data\n";
       return res;
    }
    for (int i = 0; i < fibers.size(); i++) {
      fibers[i].y        = fiberData(0+i*5);
      fibers[i].z        = fiberData(1+i*5);
      fibers[i].area     = fiberData(2+i*5);
      fibers[i].temp     = fiberData(3+i*5);
      fibers[i].temp_max = fiberData(4+i*5);
    }       

    for (int i=0; i<fibers.size(); i++) {
      int classTag = materialData(2*i);
      int dbTag = materialData(2*i+1);

      // if material pointed to is blank or not of corrcet type, 
      // release old and create a new one
      if (fibers[i].material == 0)
        fibers[i].material = theBroker.getNewUniaxialMaterial(classTag);

      else if (fibers[i].material->getClassTag() != classTag) {
        delete fibers[i].material;
        fibers[i].material = theBroker.getNewUniaxialMaterial(classTag);      
      }

      if (fibers[i].material == 0) {
        opserr << "FiberSection3d::recvSelf -- failed to allocate double array for material data\n";
        exit(-1);
      }

      fibers[i].material->setDbTag(dbTag);
      res += fibers[i].material->recvSelf(commitTag, theChannel, theBroker);
    }

    QzBar = 0.0;
    QyBar = 0.0;
    ABar = 0.0;

    computeCentroid = data(5) ? true : false;

    // Recompute centroid
    for (int i = 0; computeCentroid && i < fibers.size(); i++) {
      double yLoc = fibers[i].y;
      double zLoc = fibers[i].z;
      double Area = fibers[i].area;
      ABar  += Area;
      QzBar += yLoc*Area;
      QyBar += zLoc*Area;
    }

    if (computeCentroid) {
      yBar = QzBar/ABar;
      zBar = QyBar/ABar;
    } else {
      yBar = 0.0;
      zBar = 0.0;      
    }
  }   

  return res;
}

void
FiberSection3dThermal::Print(OPS_Stream &s, int flag)
{
  if (flag == 2) {
    for (int i = 0; i < fibers.size(); i++) {
      s << -fibers[i].y << " "  << fibers[i].z << " "  << fibers[i].area << " " ;
      s << fibers[i].material->getStress() << " "  << fibers[i].material->getStrain() << "\n";
    }
  } 


  else if (flag == 3) {
    for (int i = 0; i < fibers.size(); i++) {
      s << fibers[i].material->getTag() << " " << fibers[i].y << " "  << fibers[i].z << " "  << fibers[i].area << " " ;
      s << fibers[i].material->getStress() << " "  << fibers[i].material->getStrain() << "\n";
    } 
  }
    
  else if (flag == 4) {
    for (int i = 0; i < fibers.size(); i++) {
      s << "add fiber # " << i+1 << " using material # " << fibers[i].material->getTag() << " to section # 1\n";
      s << "fiber_cross_section = " << fibers[i].area << "*m^2\n";
      s << "fiber_location = (" << fibers[i].y << "*m, " << fibers[i].z << "*m);\n\n";
    }
  }

  else if (flag == OPS_PRINT_PRINTMODEL_JSON) { 
      s << TaggedObject::JsonPropertyIndent << "{";
      s << "\"name\": \"" << this->getTag() << "\", ";
      s << "\"type\": \"" << this->getClassType() << "\", ";

      if (theTorsion != 0)
        s << "\"torsion\": " << theTorsion->getInitialTangent() << ", ";

      s << "\"fibers\": [\n";
      for (int i = 0; i < fibers.size(); i++) {
            s << TaggedObject::JsonPropertyIndent 
              << "\t{\"coord\": [" << fibers[i].y << ", " 
                                   << fibers[i].z << "], ";
            s << "\"area\": " << fibers[i].area << ", ";
            s << "\"material\": " << fibers[i].material->getTag();
            if (i < fibers.size() - 1)
                  s << "},\n";
            else
                  s << "}\n";
      }
      s << TaggedObject::JsonPropertyIndent << "]}";
      return;
  }

  else {
    s << "\nFiberSection3dThermal, tag: " << this->getTag() << "\n";
    s << "\tSection code: " << code;
    s << "\tNumber of Fibers: " << static_cast<int>(fibers.size()) << "\n";
    s << "\tCentroid: (" << yBar << ", " << zBar << ')' << "\n";
    if (theTorsion != 0)
        theTorsion->Print(s, flag); 

    if (flag == 1) {
      for (int i = 0; i < fibers.size(); i++) {
        s << "\nLocation (y, z) = (" << -fibers[i].y << ", " << fibers[i].z << ")";
        s << "\nArea = " << fibers[i].area << "\n";
      fibers[i].material->Print(s, flag);
      }
    }
  }
}

Response*
FiberSection3dThermal::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;
  
  if (argc > 2 && strcmp(argv[0],"fiber") == 0) {

    int key = fibers.size();
    int passarg = 2;
    
    if (argc <= 3)        {  // fiber number was input directly
      
      key = atoi(argv[1]);
      
    } else if (argc > 4) {         // find fiber closest to coord. with mat tag
      int matTag = atoi(argv[3]);
      double yCoord = atof(argv[1]);
      double zCoord = atof(argv[2]);
      double closestDist = 0.0;
      double ySearch, zSearch, dy, dz;
      double distance;
      int j;
      
      // Find first fiber with specified material tag
      for (j = 0; j < fibers.size(); j++) {
        if (matTag == fibers[j].material->getTag()) {
          ySearch = -fibers[j].y;
          zSearch =  fibers[j].z;
          dy = ySearch-yCoord;
          dz = zSearch-zCoord;
          closestDist = sqrt(dy*dy + dz*dz);
          key = j;
          break;
        }
      }
      
      // Search the remaining fibers
      for ( ; j < fibers.size(); j++) {
        if (matTag == fibers[j].material->getTag()) {
          ySearch = -fibers[j].y;
          zSearch =  fibers[j].z;
          dy = ySearch-yCoord;
          dz = zSearch-zCoord;
          distance = sqrt(dy*dy + dz*dz);
          if (distance < closestDist) {
            closestDist = distance;
            key = j;
          }
        }
      }
      passarg = 4;
    }
    
    else {                  // fiber near-to coordinate specified
      double yCoord = atof(argv[1]);
      double zCoord = atof(argv[2]);
      double closestDist;
      double ySearch, zSearch, dy, dz;
      double distance;
      ySearch = -fibers[0].y;
      zSearch =  fibers[0].z;
      dy = ySearch-yCoord;
      dz = zSearch-zCoord;
      closestDist = sqrt(dy*dy + dz*dz);
      key = 0;
      for (int j = 1; j < fibers.size(); j++) {
        ySearch = -fibers[j].y;
        zSearch =  fibers[j].z;
        dy = ySearch-yCoord;
        dz = zSearch-zCoord;
        distance = sqrt(dy*dy + dz*dz);
        if (distance < closestDist) {
          closestDist = distance;
          key = j;
        }
      }
      passarg = 3;
    }
    
    if (key < fibers.size() && key >= 0) {
      output.tag("FiberOutput");
      output.attr("yLoc", fibers[key].y);
      output.attr("zLoc", fibers[key].z);
      output.attr("area", fibers[key].area);
      
      theResponse = fibers[key].material->setResponse(&argv[passarg], argc-passarg, output);
      
      output.endTag();
    }
  
  } else if (strcmp(argv[0],"fiberData") == 0) {
    int numData = fibers.size()*5;
    for (int j = 0; j < fibers.size(); j++) {
      output.tag("FiberOutput");
      output.attr("yLoc", fibers[j].y);
      output.attr("zLoc", fibers[j].z);
      output.attr("area", fibers[j].area);    
      output.tag("ResponseType","yCoord");
      output.tag("ResponseType","zCoord");
      output.tag("ResponseType","area");
      output.tag("ResponseType","stress");
      output.tag("ResponseType","strain");
      output.endTag();
    }
    Vector theResponseData(numData);
    theResponse = new MaterialResponse(this, 5, theResponseData);
  }

  if (theResponse == 0)
    return FrameSection::setResponse(argv, argc, output);

  return theResponse;
}


int
FiberSection3dThermal::getResponse(int responseID, Information &sectInfo)
{
  if (responseID == 5) {
    int numData = 5*fibers.size();
    Vector data(numData);
    int count = 0;
    for (int j = 0; j < fibers.size(); j++) {
      double yLoc, zLoc, A, stress, strain;
      yLoc = -fibers[j].y;
      zLoc = fibers[j].z;
      A = fibers[j].area;
      stress = fibers[j].material->getStress();
      strain = fibers[j].material->getStrain();
      data(count) = yLoc; data(count+1) = zLoc; data(count+2) = A;
      data(count+3) = stress; data(count+4) = strain;
      count += 5;
    }
    return sectInfo.setVector(data);
  } else
    return FrameSection::getResponse(responseID, sectInfo);
}

int
FiberSection3dThermal::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 3)
    return -1;


  int result = -1;

  // A material parameter
  if (strstr(argv[0],"material") != 0) {

    // Get the tag of the material
    int paramMatTag = atoi(argv[1]);

    // Loop over fibers to find the right material(s)
    int ok = 0;
    for (int i = 0; i < fibers.size(); i++)
      if (paramMatTag == fibers[i].material->getTag()) {
        ok = fibers[i].material->setParameter(&argv[2], argc-2, param);
        if (ok != -1)
          result = ok;
      }
    
    if (paramMatTag == theTorsion->getTag()) {
        ok = theTorsion->setParameter(&argv[2], argc-2, param);
        if (ok != -1)
          result = ok;
    }
    return result;
  }    

  // Check if it belongs to the section integration
  else if (strstr(argv[0],"integration") != 0) {
    return -1;
  }

  int ok = 0;

  // loop over every material
  for (int i = 0; i < fibers.size(); i++) {
    ok = fibers[i].material->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }

  // Don't really need to do this in "default" mode
  //ok = theTorsion->setParameter(argv, argc, param);
  //if (ok != -1)
  //  result = ok;

  return result;
}

const Vector &
FiberSection3dThermal::getSectionDeformationSensitivity(int gradIndex)
{
        static Vector dummy(3);
        dummy.Zero();
        if (SHVs !=0) {
                dummy(0) = (*SHVs)(0,gradIndex);
                dummy(1) = (*SHVs)(1,gradIndex);
                dummy(2) = (*SHVs)(2,gradIndex);
        }
        return dummy;
}


const Vector &
FiberSection3dThermal::getStressResultantSensitivity(int gradIndex, bool conditional)
{
  static Vector ds(4);
  
  ds.Zero();

  double  stressGradient;
  int loc = 0;


  for (int i = 0; i < fibers.size(); i++) {
    double y = fibers[i].y - yBar;
    double z = fibers[i].z - zBar;
    double A = fibers[i].area;
    stressGradient = fibers[i].material->getStressSensitivity(gradIndex,conditional);
    stressGradient *=  A;
    ds(0) += stressGradient;
    ds(1) += stressGradient * y;
    ds(2) += stressGradient * z;
  }
  
  ds(3) = theTorsion->getStressSensitivity(gradIndex, conditional);

  return ds;
}

const Matrix &
FiberSection3dThermal::getSectionTangentSensitivity(int gradIndex)
{
  static Matrix something(4,4);

  something.Zero();

  something(3,3) = theTorsion->getTangentSensitivity(gradIndex);
  
  return something;
}

int
FiberSection3dThermal::commitSensitivity(const Vector& defSens, int gradIndex, int numGrads)
{

  // here add SHVs to store the strain sensitivity.

  if (SHVs == 0) {
    SHVs = new Matrix(4,numGrads);
  }

  (*SHVs)(0,gradIndex) = defSens(0);
  (*SHVs)(1,gradIndex) = defSens(1);
  (*SHVs)(2,gradIndex) = defSens(2);
  (*SHVs)(3,gradIndex) = defSens(3);
  int loc = 0;

  double d0 = defSens(0);
  double d1 = defSens(1);
  double d2 = defSens(2);
  double d3 = defSens(3);
  for (int i = 0; i < fibers.size(); i++) {
      double y = fibers[i].y - yBar;
      double z = fibers[i].z - zBar;

      double strainSens = d0 + y*d1 + z*d2;

      fibers[i].material->commitSensitivity(strainSens,gradIndex,numGrads);
  }

  theTorsion->commitSensitivity(d3, gradIndex, numGrads);

  return 0;
}

// AddingSensitivity:END ///////////////////////////////////


double
FiberSection3dThermal::determineFiberTemperature(const Vector& DataMixed, double fiberLocy, double fiberLocz)
{
    double FiberTemperature = 0;
    if (DataMixed.Size()==18) {
      if ( fabs(DataMixed(1)) <= 1e-10 && fabs(DataMixed(17)) <= 1e-10 ) //no tempe load
      {
              return 0 ;
      }

      double dataTempe[18]; //PK changed 18 to 27 to pass max temps
      for (int i = 0; i < 18; i++) {
              dataTempe[i] = DataMixed(i);
      }

      if (  fiberLocy <= dataTempe[1])
              opserr <<"FiberSection2dThermal::setTrialSectionDeformationTemperature -- fiber loc is out of the section";
      else if (fiberLocy <= dataTempe[3])
              FiberTemperature = dataTempe[0] - (dataTempe[1] - fiberLocy) * (dataTempe[0] - dataTempe[2])/(dataTempe[1] - dataTempe[3]);

      else if (   fiberLocy <= dataTempe[5] )
      {
              FiberTemperature = dataTempe[2] - (dataTempe[3] - fiberLocy) * (dataTempe[2] - dataTempe[4])/(dataTempe[3] - dataTempe[5]);
      }
      else if ( fiberLocy <= dataTempe[7] )
      {
              FiberTemperature = dataTempe[4] - (dataTempe[5] - fiberLocy) * (dataTempe[4] - dataTempe[6])/(dataTempe[5] - dataTempe[7]);
      }
      else if ( fiberLocy <= dataTempe[9] )
      {
              FiberTemperature = dataTempe[6] - (dataTempe[7] - fiberLocy) * (dataTempe[6] - dataTempe[8])/(dataTempe[7] - dataTempe[9]);
      }
      else if (fiberLocy <= dataTempe[11] )
      {
              FiberTemperature = dataTempe[8] - (dataTempe[9] - fiberLocy) * (dataTempe[8] - dataTempe[10])/(dataTempe[9] - dataTempe[11]);
      }
      else if (fiberLocy <= dataTempe[13] )
      {
              FiberTemperature = dataTempe[10] - (dataTempe[11] - fiberLocy) * (dataTempe[10] - dataTempe[12])/(dataTempe[11] - dataTempe[13]);
      }
      else if (fiberLocy <= dataTempe[15] )
      {
              FiberTemperature = dataTempe[12] - (dataTempe[13] - fiberLocy) * (dataTempe[12] - dataTempe[14])/(dataTempe[13] - dataTempe[15]);
      }
      else if ( fiberLocy <= dataTempe[17] )
          FiberTemperature = dataTempe[14] - (dataTempe[15] - fiberLocy) * (dataTempe[14] - dataTempe[16])/(dataTempe[15] - dataTempe[17]);
      else
          opserr <<"FiberSection3dThermal::setTrialSectionDeformation -- fiber loc " <<fiberLocy<<" is out of the section"<<"\n";
  }
  else if(DataMixed.Size()==25) {
  //---------------if temperature Data has 25 elements--------------------
      double dataTempe[25]; //
      for (int i = 0; i < 25; i++) { //
              dataTempe[i] = DataMixed(i);
      }

      if ( fabs(dataTempe[0]) <= 1e-10 && fabs(dataTempe[2]) <= 1e-10 && fabs(dataTempe[10]) <= 1e-10 && fabs(dataTempe[11]) <= 1e-10) //no tempe load
      {
              return 0;
      }

// calculate the fiber tempe, T=T1-(Y-Y1)*(T1-T2)/(Y1-Y2)
// first for bottom flange if existing
      if (  fiberLocy <= dataTempe[1]) {
            if (fiberLocz <= dataTempe[12]){
            opserr<<"WARNING: FiberSection3dThermal failed to find the fiber with locy: "<<fiberLocy <<" , locZ: "<<fiberLocz <<"\n";
            }
            else if (fiberLocz<= dataTempe[15]){
            FiberTemperature = dataTempe[10] - (dataTempe[10] - dataTempe[13])*(dataTempe[12] - fiberLocz) /(dataTempe[12] - dataTempe[15]);
            }
            else if (fiberLocz<= dataTempe[18]){
            FiberTemperature = dataTempe[13] - (dataTempe[13] - dataTempe[16])*(dataTempe[15] - fiberLocz) /(dataTempe[15] - dataTempe[18]);
            }
            else if (fiberLocz<= dataTempe[21]){
            FiberTemperature = dataTempe[16] - (dataTempe[16] - dataTempe[19])*(dataTempe[18] - fiberLocz) /(dataTempe[18] - dataTempe[21]);
            }
            else if (fiberLocz<= dataTempe[24]){
            FiberTemperature = dataTempe[19] - (dataTempe[19] - dataTempe[22])*(dataTempe[21] - fiberLocz) /(dataTempe[21] - dataTempe[24]);
            }
            else {
            opserr<<"WARNING: FiberSection3dThermal failed to find the fiber with locy: "<<fiberLocy <<" , locZ: "<<fiberLocz <<"\n";
            }
      }
      else if (fiberLocy <= dataTempe[3])
      {
          FiberTemperature = dataTempe[0] - (dataTempe[1] - fiberLocy) * (dataTempe[0] - dataTempe[2])/(dataTempe[1] - dataTempe[3]);
      }
      else if (   fiberLocy <= dataTempe[5] )
      {
          FiberTemperature = dataTempe[2] - (dataTempe[3] - fiberLocy) * (dataTempe[2] - dataTempe[4])/(dataTempe[3] - dataTempe[5]);
      }
      else if ( fiberLocy <= dataTempe[7] )
      {
          FiberTemperature = dataTempe[4] - (dataTempe[5] - fiberLocy) * (dataTempe[4] - dataTempe[6])/(dataTempe[5] - dataTempe[7]);
      }
      else if ( fiberLocy <= dataTempe[9] )
          FiberTemperature = dataTempe[6] - (dataTempe[7] - fiberLocy) * (dataTempe[6] - dataTempe[8])/(dataTempe[7] - dataTempe[9]);
      else {
          if (fiberLocz <= dataTempe[12]){
            opserr<<"WARNING: FiberSection3dThermal failed to find the fiber with locy: "<<fiberLocy <<" , locZ: "<<fiberLocz <<"\n";
          }
          else if (fiberLocz<= dataTempe[15]){
            FiberTemperature = dataTempe[11] - (dataTempe[11] - dataTempe[14])*(dataTempe[12] - fiberLocz) /(dataTempe[12] - dataTempe[15]);
          }
          else if (fiberLocz<= dataTempe[18]){
            FiberTemperature = dataTempe[14] - (dataTempe[14] - dataTempe[17])*(dataTempe[15] - fiberLocz) /(dataTempe[15] - dataTempe[18]);
          }
          else if (fiberLocz<= dataTempe[21]){
            FiberTemperature = dataTempe[17] - (dataTempe[17] - dataTempe[20])*(dataTempe[18] - fiberLocz) /(dataTempe[18] - dataTempe[21]);
          }
          else if (fiberLocz<= dataTempe[24]){
            FiberTemperature = dataTempe[20] - (dataTempe[20] - dataTempe[23])*(dataTempe[21] - fiberLocz) /(dataTempe[21] - dataTempe[24]);
          }
          else {
             opserr<<"WARNING: FiberSection3dThermal failed to find the fiber with locy: "<<fiberLocy <<" , locZ: "<<fiberLocz <<"\n";
          }
      }
  }
  else if (DataMixed.Size() == 35) {
      //---------------GR mod - if temperature Data has 35 elements--------------------
      double py[5], pz[5];
      for (int i = 0; i < 5; i++) {
          py[i] = DataMixed(i);
          pz[i] = DataMixed(5 + i);
      }
      double dataTempe[5][5];
      for (int i = 0; i < 5; i++) {
          for (int j = 0; j < 5; j++) {
              dataTempe[i][j] = DataMixed(10 + 5 * i + j);
          }
      }
      // check grid corners
      if (fabs(dataTempe[0][0]) <= 1e-10 && fabs(dataTempe[4][4]) <= 1e-10 && fabs(dataTempe[4][0]) <= 1e-10 && fabs(dataTempe[0][4]) <= 1e-10) // no tempe load
          return 0;

      // calculate the fiber temperature, weighted with inverse of the distance from grid nodes
      // check if coords are inside the grid, otherwise first or last temperature is returned
      if (fiberLocy < py[0] || fiberLocz < pz[0])      return dataTempe[0][0];
      if (fiberLocy > py[4] || fiberLocz > pz[4])      return dataTempe[4][4];
      // first, find nearest grid points
      for (int i = 1; i < 5; i++) {
          for (int j = 1; j < 5; j++) {
              if ((fiberLocy >= py[i - 1] && fiberLocy <= py[i]) && (fiberLocz >= pz[j - 1] && fiberLocz <= pz[j])) {
                  double sy[4], sz[4], sT[4], d[4], dsum = 0, Tdsum = 0;
                  // selecting the 4 points of the grid, ordered anti-clockwise
                  sy[0] = py[i - 1]; sy[1] = py[i - 1]; sy[2] = py[i]; sy[3] = py[i];
                  sz[0] = pz[j - 1]; sz[1] = pz[j]; sz[2] = pz[j]; sz[3] = pz[j - 1];
                  // select grid temperatures
                  sT[0] = dataTempe[i - 1][j - 1]; sT[1] = dataTempe[i - 1][j]; sT[2] = dataTempe[i][j]; sT[3] = dataTempe[i][j - 1];
                  // calculate distance
                  for (int k = 0; k < 4; k++) {
                      d[k] = sqrt(pow(fiberLocy - sy[k], 2) + pow(fiberLocz - sz[k], 2));
                      if (d[k] == 0) d[k] = 1.e-6;
                      dsum += 1.0 / d[k]; Tdsum += sT[k] / d[k];
                  }
                  if (dsum == 0) return 0;
                  return Tdsum / dsum;
              }
          }
      }

  }
  return FiberTemperature;
}

} // namespace
