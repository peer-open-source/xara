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
// Description: This file contains the class implementation of FrameFiberSection3d.
//
// Written: cmp
// Created: Summer 2024
//
#include <memory>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <Channel.h>
#include <Vector.h>
#include <Matrix.h>
#include <classTags.h>
#include <FrameFiberSection3d.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>
#include <SensitiveResponse.h>
typedef SensitiveResponse<FrameSection> SectionResponse;
#include <UniaxialMaterial.h>
#include <ElasticMaterial.h>
#include "ElasticLinearFrameSection3d.h"

#include "FiberResponse.h"

ID FrameFiberSection3d::code(FrameFiberSection3d::nsr);

FrameFiberSection3d::FrameFiberSection3d(int tag,
                                         int num,
                                         const Frame::Shape& shape_data,
                                         double mass,
                                         bool use_mass)
  : FrameSection(tag, SEC_TAG_FrameFiberSection3d, mass, use_mass)
  , es{}, sr{}, tangent{}
  , e(es), s(sr), K_wrap(tangent.matrix)
  , shape(std::make_shared<Frame::Shape>(shape_data))
  , fibers(std::make_shared<std::vector<FiberData>>())
{
  code(inx) = FrameStress::N;
  code(iny) = FrameStress::Vy;
  code(inz) = FrameStress::Vz;
  code(imx) = FrameStress::T;
  code(imy) = FrameStress::My;
  code(imz) = FrameStress::Mz;
  code(iwx) = FrameStress::Bimoment;
  code(iwy) = FrameStress::By;
  code(iwz) = FrameStress::Bz;
  code(ivx) = FrameStress::Bishear;
  code(ivy) = FrameStress::Qy;
  code(ivz) = FrameStress::Qz;

  fibers->reserve(num);
  materials.reserve(num);

  this->initializeShear(tangent, *shape);
}


FrameFiberSection3d::FrameFiberSection3d(const FrameFiberSection3d &other)
  : FrameSection(other.getTag(), other.getClassTag())
  , es{}, sr{}, tangent{}
  , e(es), s(sr), K_wrap(tangent.matrix)
  , shape(other.shape)
  , fibers(other.fibers)
{
  materials.reserve(other.materials.size());
  for (auto& mat : other.materials)
    materials.push_back(mat->getCopy());
  this->initializeShear(tangent, *shape);
}

FrameFiberSection3d::~FrameFiberSection3d()
{
  // Clean up the materials
  for (auto& material : materials)
    delete material;
}


int
FrameFiberSection3d::addFiber(UniaxialMaterial &theMat, 
                              double Area, double yLoc, double zLoc)
{
  FiberData fiber{yLoc, zLoc, Area};
  fibers->push_back(fiber);
  materials.push_back(theMat.getCopy());
  return materials.size() - 1;
}



int
FrameFiberSection3d::setTrialSectionDeformation(const Vector &deforms)
{
  e = deforms;
  return this->updateAxial(State::Pres, tangent, sr);
}

int
FrameFiberSection3d::updateAxial(const State state_flag, 
                                 Tangent& tangent, VectorND<nsr>& sr) const
{
  sr.zero();
  tangent.zeroAxial();
  MatrixND<nsr,nsr> &ks = tangent.matrix;

  const double e0 = e(inx), // u'
               kz = e(imz),
               ky = e(imy);

  int res = 0;
  const int nf = fibers->size();
  for (int i = 0; i < nf; i++) {
    const FiberData& fiber = (*fibers)[i];

    const double y  = fiber.y;
    const double z  = fiber.z;
    const double A  = fiber.area;

    // Determine material strain and set it
    double strain = e0 - y*kz + z*ky;
    double tangent, stress;
    res += materials[i]->setTrial(strain, stress, tangent);

    double EA = tangent * A;

    ks( inx, inx) +=     EA;
    ks( inx, imz) +=  -y*EA;
    ks( inx, imy) +=   z*EA;

    ks( imz, imz) +=  y*y*EA; // 5
    ks( imy, imy) +=  z*z*EA; // 10
    ks( imz, imy) += -y*z*EA;

    double fs0 = stress * A;
    sr[inx] +=    fs0;  // N
    sr[imz] += -y*fs0;  // Mz
    sr[imy] +=  z*fs0;  // My
  }

  // Fill in the symmetric terms
  ks(imz, inx) = ks(inx, imz);
  ks(imy, inx) = ks(inx, imy);
  ks(imy, imz) = ks(imz, imy);

  return res;
}

int 
FrameFiberSection3d::initializeShear(Tangent& tangent, Frame::Shape& shape)
{
  tangent.zeroShear();
  ElasticLinearFrameSection3d elastic(0, shape, 0.0, false);
  tangent.matrix = elastic.getFullTangent(State::Init);
  return 0;
}


const Matrix&
FrameFiberSection3d::getInitialTangent()
{
  static Tangent T{};
  static Matrix kInitial(T.matrix);
  this->initializeShear(T, *shape);
  this->updateAxial(State::Init, T, sr);
  return kInitial;
}

const Vector&
FrameFiberSection3d::getSectionDeformation()
{
  return e;
}

const Matrix&
FrameFiberSection3d::getSectionTangent()
{
  return K_wrap;
}

const Vector&
FrameFiberSection3d::getStressResultant()
{
  return s;
}

FrameSection*
FrameFiberSection3d::getFrameCopy()
{
  return new FrameFiberSection3d(*this);
}


const ID&
FrameFiberSection3d::getType()
{
  return code;
}

int
FrameFiberSection3d::getOrder() const
{
  return nsr;
}

int
FrameFiberSection3d::commitState()
{
  int err = 0;
  const int numFibers = fibers->size();
  for (int i = 0; i < numFibers; i++)
    err += materials[i]->commitState();

  return err;
}


int
FrameFiberSection3d::revertToLastCommit()
{
  int err = 0;
  for (auto& material : materials)
    err += material->revertToLastCommit();

  return err;
}

int
FrameFiberSection3d::revertToStart()
{
  // revert the fibers to start    
  int err = 0;
  for (auto& material : materials)
    err += material->revertToStart();

  return err;
}

double
FrameFiberSection3d::getEnergy() const
{
  double energy = 0;
  const int nf = fibers->size();
  for (int i = 0; i < nf; i++) {
      double A = (*fibers)[i].area;
      energy += A * materials[i]->getEnergy();
  }
  return energy;
}

int
FrameFiberSection3d::getIntegral(Field field, State state, double& value) const
{
  value = 0.0;
  const int numFibers = fibers->size();
  switch (field) {
    case Field::Unit:
      for (int i=0; i<numFibers; i++) {
        const double A  = (*fibers)[i].area;
        value += A;
      }
      return 0;

    case Field::Density:
      // First check if density has been specified for the section
      if (this->FrameSection::getIntegral(field, state, value) == 0) 
        return 0;

      for (int i=0; i<numFibers; i++) {
        double density;
        const double A  = (*fibers)[i].area;
        if (materials[i]->getRho() != 0)
          value += A*density;
        else
          return -1;
      }
      return 0;

    case Field::UnitY: // TODO: Centroid
      for (int i=0; i<numFibers; i++) {
        const double A  = (*fibers)[i].area;
        const double y  = (*fibers)[i].y;
        value += A*y;
      }
      return 0;

    case Field::UnitZ: // TODO: Centroid
      for (int i=0; i<numFibers; i++) {
        const double A  = (*fibers)[i].area;
        const double z  = (*fibers)[i].z;
        value += A*z;
      }
      return 0;



    case Field::UnitYY:
      for (int i=0; i<numFibers; i++) {
        const double A  = (*fibers)[i].area;
        const double y  = (*fibers)[i].y;
        value += A*y*y;
      }
      return 0;

    case Field::UnitZZ:
      for (int i=0; i<numFibers; i++) {
        const double A  = (*fibers)[i].area;
        const double z  = (*fibers)[i].z;
        value += A*z*z;
      }

    default:
      return -1;
  }
  return -1;
}

int
FrameFiberSection3d::sendSelf(int commitTag, Channel &)
{
  return -1;
}

int
FrameFiberSection3d::recvSelf(int commitTag, Channel &, FEM_ObjectBroker &theBroker)
{
  return -1;
}


Response*
FrameFiberSection3d::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = nullptr;
  const int numFibers = fibers->size();
  
  if (argc > 2 && strcmp(argv[0],"fiber") == 0) {

    int key = numFibers;
    int passarg = 2;
    
    if (argc <= 3)      {  // fiber number was input directly
      key = atoi(argv[1]);

    } else if (argc > 4) {         // find fiber closest to coord. with mat tag
      int matTag = atoi(argv[3]);
      double yCoord = atof(argv[1]);
      double zCoord = atof(argv[2]);
      double closestDist = 0.0;
      double ySearch, zSearch, dy, dz;
      double distance;

      // Find first fiber with specified material tag
      int j;
      for (j = 0; j < numFibers; j++) {
        if (matTag == materials[j]->getTag()) {
          ySearch = (*fibers)[j].y;
          zSearch = (*fibers)[j].z;
          dy = ySearch-yCoord;
          dz = zSearch-zCoord;
          closestDist = dy*dy + dz*dz;
          key = j;
          break;
        }
      }

      // Search the remaining fibers
      for ( ; j < numFibers; j++) {
        if (matTag == materials[j]->getTag()) {
          ySearch = (*fibers)[j].y;
          zSearch = (*fibers)[j].z;

          dy = ySearch - yCoord;
          dz = zSearch - zCoord;
          distance = dy*dy + dz*dz;
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
      ySearch = (*fibers)[0].y;
      zSearch = (*fibers)[0].z;
      // ySearch = yLocs[0];
      // zSearch = zLocs[0];
      dy = ySearch-yCoord;
      dz = zSearch-zCoord;
      closestDist = sqrt(dy*dy + dz*dz);
      key = 0;
      for (int j = 1; j < numFibers; j++) {
        ySearch = (*fibers)[j].y;
        zSearch = (*fibers)[j].z;                      
        dy = ySearch - yCoord;
        dz = zSearch - zCoord;
        distance = sqrt(dy*dy + dz*dz);
        if (distance < closestDist) {
          closestDist = distance;
          key = j;
        }
      }
      passarg = 3;
    }
    
    if (key < numFibers && key >= 0) {
      output.tag("FiberOutput");
      output.attr("yLoc",(*fibers)[key].y);
      output.attr("zLoc",(*fibers)[key].z);
      output.attr("area",(*fibers)[key].area);
      
      theResponse = materials[key]->setResponse(&argv[passarg], argc-passarg, output);
      
      output.endTag();
    }
  
  } else if (strcmp(argv[0],"fiberData") == 0) {
    int numData = numFibers*5;
    for (int j = 0; j < numFibers; j++) {
      output.tag("FiberOutput");
      output.attr("yLoc", (*fibers)[j].y);
      output.attr("zLoc", (*fibers)[j].z);
      output.attr("area", (*fibers)[j].area);    
      output.tag("ResponseType","yCoord");
      output.tag("ResponseType","zCoord");
      output.tag("ResponseType","area");
      output.tag("ResponseType","stress");
      output.tag("ResponseType","strain");
      output.endTag();
    }
    Vector theResponseData(numData);
    theResponse = new SectionResponse(*this, FiberResponse::FiberData, theResponseData);

  } else if (strcmp(argv[0],"fiberData2") == 0) {
    int numData = numFibers*6;
    for (int j = 0; j < numFibers; j++) {
      output.tag("FiberOutput");
      output.attr("yLoc", (*fibers)[j].y);
      output.attr("zLoc", (*fibers)[j].z);
      output.attr("area", (*fibers)[j].area);    
      output.attr("material", materials[j]->getTag());
      output.tag("ResponseType","yCoord");
      output.tag("ResponseType","zCoord");
      output.tag("ResponseType","area");
      output.tag("ResponseType","material");
      output.tag("ResponseType","stress");
      output.tag("ResponseType","strain");
      output.endTag();
    }
    Vector theResponseData(numData);
    theResponse = new SectionResponse(*this, FiberResponse::FiberData02, theResponseData);

  } else if ((strcmp(argv[0],"numFailedFiber") == 0) || 
             (strcmp(argv[0],"numFiberFailed") == 0)) {
    int count = 0;
    theResponse = new SectionResponse(*this, 6, count);

  } else if ((strcmp(argv[0],"sectionFailed") == 0) ||
             (strcmp(argv[0],"hasSectionFailed") == 0) ||
             (strcmp(argv[0],"hasFailed") == 0)) {

    int count = 0;
    theResponse = new SectionResponse(*this, 7, count);
  }
  //by SAJalali
  else if ((strcmp(argv[0], "energy") == 0) || (strcmp(argv[0], "Energy") == 0)) {
    output.tag("SectionOutput");
    output.attr("secType", this->getClassType());
    output.attr("secTag", this->getTag());
    output.tag("ResponseType", "energy");
    theResponse = new SectionResponse(*this, 10, getEnergy());
    output.endTag();
  }

  if (theResponse == nullptr)
    return FrameSection::setResponse(argv, argc, output);

  return theResponse;
}


int 
FrameFiberSection3d::getResponse(int responseID, Information &sectInfo)
{
  const int numFibers = fibers->size();
  // Just call the base class method ... don't need to define
  // this function, but keeping it here just for clarity
  if (responseID == FiberResponse::FiberData) {
    int numData = 5*numFibers;
    Vector data(numData);
    int count = 0;
    for (int j = 0; j < numFibers; j++) {
      double yLoc = (*fibers)[j].y;
      double zLoc = (*fibers)[j].z;
      double A = (*fibers)[j].area;
      double stress = materials[j]->getStress();
      double strain = materials[j]->getStrain();
      data(count) = yLoc; data(count+1) = zLoc; data(count+2) = A;
      data(count+3) = stress; data(count+4) = strain;
      count += 5;
    }
    return sectInfo.setVector(data);
  } 
  else if (responseID == FiberResponse::FiberData02) {
    int numData = 6*numFibers;
    Vector data(numData);
    int count = 0;
    for (int j = 0; j < numFibers; j++) {
      data(count)   = (*fibers)[j].y; // y
      data(count+1) = (*fibers)[j].z; // z
      data(count+2) = (*fibers)[j].area; // A
      data(count+3) = (double)materials[j]->getTag();
      data(count+4) = materials[j]->getStress();
      data(count+5) = materials[j]->getStrain();	    
      count += 6;
    }
    return sectInfo.setVector(data);		  

  } else  if (responseID == 6) {
    int count = 0;
    for (int j = 0; j < numFibers; j++) {    
      if (materials[j]->hasFailed() == true)
      count++;
    }
    return sectInfo.setInt(count);
  } else  if (responseID == 7) {
    int count = 0;
    for (int j = 0; j < numFibers; j++) {    
      if (materials[j]->hasFailed() == true) {
      count+=1;
      }
    }
    if (count == numFibers)
      count = 1;
    else
      count = 0;

    return sectInfo.setInt(count);
  } 
  else  if (responseID == 10) {
    return sectInfo.setDouble(getEnergy());
  }

  return FrameSection::getResponse(responseID, sectInfo);
}

int
FrameFiberSection3d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  int result = -1;
  const int numFibers = fibers->size();

  // A material parameter
  if (strstr(argv[0],"material") != 0) {

    // Get the tag of the material
    int paramMatTag = atoi(argv[1]);

    // Loop over fibers to find the right material(s)
    int ok = 0;
    for (int i = 0; i < numFibers; i++)
      if (paramMatTag == materials[i]->getTag()) {
      ok = materials[i]->setParameter(&argv[2], argc-2, param);
      if (ok != -1)
        result = ok;
      }

    return result;
  }

  else if (strstr(argv[0],"fiber") != 0) {
    // ... fiber $fiberID $field
    if (argc < 3)
      return -1;
    int fiberID = atoi(argv[1]);
    if (fiberID < 0 || fiberID >= numFibers)
      return -1;
    int field;
    if (strcmp(argv[2],"y") == 0)
      field = Param::FiberY;
    else if (strcmp(argv[2],"z") == 0)
      field = Param::FiberZ;
    else if (strcmp(argv[2],"area") == 0)
      field = Param::FiberArea;
    else
      return -1;

    return param.addObject(Param::FiberFieldBase + fiberID*100+ field, this);
  }

  // Check if it belongs to the section integration
  else if (strstr(argv[0],"integration") != 0)
    return -1;

  int ok = 0;
  
  // loop over every material
  for (int i = 0; i < numFibers; i++) {
    ok = materials[i]->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }

  // Don't really need to do this in "default" mode
  //ok = shear->setParameter(argv, argc, param);
  //if (ok != -1)
  //  result = ok;

  return result;
}


int 
FrameFiberSection3d::updateParameter(int parameterID, Information &info)
{
  if (parameterID >= Param::FiberFieldBase) {
    int fiberID = (parameterID - Param::FiberFieldBase) / 100;
    int field = (parameterID - Param::FiberFieldBase) % 100;
    if (fiberID < 0 || fiberID >= fibers->size())
      return -1;

    double value = info.theDouble;

    switch (field) {
      case Param::FiberY:
        (*fibers)[fiberID].y = value;
        break;
      case Param::FiberZ:
        (*fibers)[fiberID].z = value;
        break;
      case Param::FiberArea:
        (*fibers)[fiberID].area = value;
        break;
      default:
        return -1;
    }
    return 0;
  }

  // loop over every material
  for (int i = 0; i < materials.size(); i++)
    materials[i]->updateParameter(parameterID, info);

  return 0;
}


const Vector &
FrameFiberSection3d::getSectionDeformationSensitivity(int gradIndex)
{
  static Vector dummy(nsr);
  
  dummy.Zero();
  
  return dummy;
}

   
const Vector &
FrameFiberSection3d::getStressResultantSensitivity(int gradIndex, bool conditional)
{
  static Vector ds(nsr);
  
  ds.Zero();
  const int numFibers = fibers->size();
  
  double stress = 0;
  double dsigdh = 0;
  double sig_dAdh = 0;
  double tangent = 0;

  for (int i = 0; i < numFibers; i++) {
#if 1
    double y  = (*fibers)[i].y;
    double z  = (*fibers)[i].z;
    double A  = (*fibers)[i].area;
    double dA=0, dy=0, dz=0;
    this->FiberGrad(i, gradIndex, dA, dy, dz);
    
    dsigdh = materials[i]->getStressSensitivity(gradIndex, conditional);

    ds(0) += dsigdh*A;
    ds(1) += -y*dsigdh*A;
    ds(2) +=  z*dsigdh*A;

    if (dA != 0.0 || dy != 0.0 ||  dz != 0.0)
      stress = materials[i]->getStress();

    if (dy != 0.0 || dz != 0.0)
      tangent = materials[i]->getTangent();

    if (dA != 0.0) {
      sig_dAdh = stress*dA;
      
      ds(0) += sig_dAdh;
      ds(1) += -y*sig_dAdh;
      ds(2) +=  z*sig_dAdh;
    }

    if (dy != 0.0)
      ds(1) += -dy * (stress*A);

    if (dz != 0.0)
      ds(2) +=  dz * (stress*A);

    static Matrix as(1,3);
    as(0,0) = 1;
    as(0,1) = -y;
    as(0,2) =  z;

    static Matrix dasdh(1,3);
    dasdh(0,1) = -dy;
    dasdh(0,2) =  dz;

    static Matrix tmpMatrix(3,3);
    tmpMatrix.addMatrixTransposeProduct(0.0, as, dasdh, tangent);

    //ds.addMatrixVector(1.0, tmpMatrix, e, A);
    ds(inx) += (tmpMatrix(0,0)*e(0) + tmpMatrix(0,1)*e(1) + tmpMatrix(0,2)*e(2))*A;
    ds(imz) += (tmpMatrix(1,0)*e(0) + tmpMatrix(1,1)*e(1) + tmpMatrix(1,2)*e(2))*A;
    ds(imy) += (tmpMatrix(2,0)*e(0) + tmpMatrix(2,1)*e(1) + tmpMatrix(2,2)*e(2))*A;
#endif
  }

  return ds;
}

const Matrix &
FrameFiberSection3d::getSectionTangentSensitivity(int gradIndex)
{
  static Matrix something(nsr,nsr);
  
  something.Zero();
  
  return something;
}


int
FrameFiberSection3d::commitSensitivity(const Vector& defSens, int gradIndex, int numGrads)
{

  double d0 = defSens(0);
  double d1 = defSens(1);
  double d2 = defSens(2);
  double d3 = defSens(3);

  //dedh = defSens;

  static double dydh[10000]{};
  static double dzdh[10000]{};

  const int numFibers = fibers->size();

  double depsdh = 0;

  for (int i = 0; i < numFibers; i++) {
    FiberData& fiber = (*fibers)[i];
    const double y = fiber.y;
    const double z = fiber.z;

    // determine material strain and set it
    depsdh = d0 - y*d1 + z*d2 - dydh[i]*e(1) + dzdh[i]*e(2);

    materials[i]->commitSensitivity(depsdh,gradIndex,numGrads);
  }

  return 0;
}


void
FrameFiberSection3d::Print(OPS_Stream &s, int flag)
{
  const int numFibers = fibers->size();
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"" << this->getClassType() << "\", ";

    double mass;
    if (this->FrameSection::getIntegral(Field::Density, State::Init, mass) == 0)
      s << "\"mass\": " << mass;

    s << "\"fibers\": [\n";
    for (int i = 0; i < numFibers; i++) {
      const FiberData& fiber = (*fibers)[i];
      s << OPS_PRINT_JSON_MATE_INDENT << "\t{";
      s << "\"coord\": [" << fiber.y << ", " << fiber.z << "], ";
      s << "\"area\": " << fiber.area << ", ";
      s << "\"material\": " << materials[i]->getTag();
      if (i < numFibers - 1)
          s << "},\n";
      else
          s << "}\n";
    }
    s << OPS_PRINT_JSON_MATE_INDENT << "]}";
    return;
  }

  if (flag == OPS_PRINT_PRINTMODEL_SECTION || flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
    s << "\nFrameFiberSection3d, tag: " << this->getTag() << "\n";
    s << "\tSection code: " << code;
    s << "\tNumber of Fibers: " << numFibers << "\n"; 

    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
      for (int i = 0; i < numFibers; i++) {
        FiberData& fiber = (*fibers)[i];
        s << "\nLocation (y, z) = (" << fiber.y << ", " << fiber.z << ")";
        s << "\nArea = " << fiber.area << "\n";
        materials[i]->Print(s, flag);
      }
    }
  }
}