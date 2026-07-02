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
// Description: This file contains the class implementation of FiberSection2d.
//
// Written: fmk
// Created: 04/04
//
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <Channel.h>
#include <Vector.h>
#include <Matrix.h>
#include <classTags.h>
#include <FiberSection2d.h>
#include <ID.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <SensitiveResponse.h>
typedef SensitiveResponse<FrameSection> SectionResponse;
#include <UniaxialMaterial.h>

#include "FiberResponse.h"

ID FiberSection2d::code(2);


// allocate memory for fibers
FiberSection2d::FiberSection2d(int tag, int num, bool compCentroid): 
  FrameSection(tag, SEC_TAG_FiberSection2d),
  theMaterials(0), 
  fibers(std::make_shared<std::vector<FiberData>>(0)),
  QzBar(0.0), ABar(0.0), yBar(0.0), computeCentroid(compCentroid),
  e(2), s(sData,2), ks(kData,2,2), dedh(2)
{
  fibers->reserve(num);
  theMaterials.reserve(num);

  sData[0] = 0.0;
  sData[1] = 0.0;

  kData[0] = 0.0;
  kData[1] = 0.0;
  kData[2] = 0.0;
  kData[3] = 0.0;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
}

// constructor for blank object that recvSelf needs to be invoked upon
FiberSection2d::FiberSection2d():
  FrameSection(0, SEC_TAG_FiberSection2d),
  theMaterials(0),
  QzBar(0.0), ABar(0.0), yBar(0.0), computeCentroid(true),
  e(2), s(sData,2), ks(kData,2,2), dedh(2)
{
  // s = new Vector(sData, 2);
  // ks = new Matrix(kData, 2, 2);

  sData[0] = 0.0;
  sData[1] = 0.0;

  kData[0] = 0.0;
  kData[1] = 0.0;
  kData[2] = 0.0;
  kData[3] = 0.0;

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
}


FiberSection2d::FiberSection2d(const FiberSection2d &other):
  FrameSection(other)
  , fibers(other.fibers)
  , QzBar(other.QzBar)
  , ABar(other.ABar)
  , yBar(other.yBar)
  , computeCentroid(other.computeCentroid),
  e(other.e), s(sData,2), ks(kData,2,2), dedh(other.dedh)
{
  theMaterials.reserve(other.theMaterials.size());
  for (int i=0; i<other.theMaterials.size(); i++)
    theMaterials.push_back(other.theMaterials[i]->getCopy());

  s = other.s;
  ks = other.ks;
}

int
FiberSection2d::addFiber(UniaxialMaterial &theMat, const double Area, const double yLoc)
{
  
  fibers->emplace_back(FiberData{Area, yLoc});
  theMaterials.push_back(theMat.getCopy());
  // Recompute centroid
  if (computeCentroid && ABar != 0.0) {
    ABar += Area;
    QzBar += yLoc*Area;
    yBar = QzBar/ABar;
  }
  
  return theMaterials.size() - 1;
}


// destructor:
FiberSection2d::~FiberSection2d()
{
  for (UniaxialMaterial* m : theMaterials)
    delete m;
}


int
FiberSection2d::getIntegral(Field field, State state, double& value) const
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
        double density = theMaterials[i]->getRho();
        const double A  = (*fibers)[i].area;
        if (density != 0)
          value += A*density;
        else
          return -1;
      }
      return 0;


    case Field::UnitYY:
    case Field::UnitCentroidYY:
      for (int i=0; i<numFibers; i++) {
        const double A  = (*fibers)[i].area;
        const double y  = (*fibers)[i].y
                        - yBar*(Field::UnitCentroidYY == field);
        value += A*y*y;
      }
      return 0;


    default:
      return -1;
  }
  return -1;
}

int
FiberSection2d::setTrialSectionDeformation (const Vector &deforms)
{

  e = deforms;

  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0; kData[3] = 0.0;
  sData[0] = 0.0; sData[1] = 0.0;

  const double d0 = deforms(0),
               d1 = deforms(1);

  
  int res = 0;
  const int numFibers = fibers->size();
  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    const double y = (*fibers)[i].y - yBar;
    const double A = (*fibers)[i].area;

    // determine material strain and set it
    double strain = d0 - y*d1;
    double tangent, stress;
    res += theMat->setTrial(strain, stress, tangent);

    double ks0 = tangent * A;
    double ks1 = ks0 * -y;
    kData[0]  += ks0;
    kData[1]  += ks1;
    kData[3]  += ks1 * -y;

    double fs0 = stress * A;
    sData[0] += fs0;
    sData[1] += fs0 * -y;
  }

  kData[2] = kData[1];

  return res;
}

const Vector&
FiberSection2d::getSectionDeformation()
{
  return e;
}

const Matrix&
FiberSection2d::getInitialTangent()
{
  static double kInitial[4];
  static Matrix kInitialMatrix(kInitial, 2, 2);
  kInitial[0] = 0.0; kInitial[1] = 0.0; kInitial[2] = 0.0; kInitial[3] = 0.0;


  const int numFibers = fibers->size();
  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    const double y = (*fibers)[i].y - yBar;
    const double A = (*fibers)[i].area;

    double tangent = theMat->getInitialTangent();

    double ks0 = tangent * A;
    double ks1 = ks0 * -y;
    kInitial[0] += ks0;
    kInitial[1] += ks1;
    kInitial[3] += ks1 * -y;
  }

  kInitial[2] = kInitial[1];

  return kInitialMatrix;
}

const Matrix&
FiberSection2d::getSectionTangent()
{
  return ks;
}

const Vector&
FiberSection2d::getStressResultant()
{
  return s;
}

FrameSection*
FiberSection2d::getFrameCopy()
{
  return new FiberSection2d(*this);
}

const ID&
FiberSection2d::getType()
{
  return code;
}

int
FiberSection2d::getOrder () const
{
  return 2;
}

int
FiberSection2d::commitState()
{
  int err = 0;

  const int numFibers = fibers->size();
  for (int i = 0; i < numFibers; i++)
    err += theMaterials[i]->commitState();

  return err;
}

int
FiberSection2d::revertToLastCommit()
{
  int err = 0;

  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0; kData[3] = 0.0;
  sData[0] = 0.0; sData[1] = 0.0;

  const int numFibers = fibers->size();
  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    const double y = (*fibers)[i].y - yBar;
    const double A = (*fibers)[i].area;

    // invoke revertToLast on the material
    err += theMat->revertToLastCommit();

    // get material stress & tangent for this strain and determine ks and fs
    double tangent = theMat->getTangent();
    double stress = theMat->getStress();
    double ks0 = tangent * A;
    double ks1 = ks0 * -y;
    kData[0] += ks0;
    kData[1] += ks1;
    kData[3] += ks1 * -y;

    double fs0 = stress * A;
    sData[0] += fs0;
    sData[1] += fs0 * -y;
  }

  kData[2] = kData[1];

  return err;
}

int
FiberSection2d::revertToStart()
{
  // revert the fibers to start    
  int err = 0;

  kData[0] = 0.0; kData[1] = 0.0; kData[2] = 0.0; kData[3] = 0.0;
  sData[0] = 0.0; sData[1] = 0.0;

  const int numFibers = fibers->size();
  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    double y = (*fibers)[i].y - yBar;
    double A = (*fibers)[i].area;

    // invoke revertToLast on the material
    err += theMat->revertToStart();

    // get material stress & tangent for this strain and determine ks and fs
    double tangent = theMat->getTangent();
    double stress = theMat->getStress();
    double ks0 = tangent * A;
    double ks1 = ks0 * -y;
    kData[0] += ks0;
    kData[1] += ks1;
    kData[3] += ks1 * -y;

    double fs0 = stress * A;
    sData[0] = fs0;
    sData[1] = fs0 * -y;
  }

  kData[2] = kData[1];

  return err;
}

int
FiberSection2d::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int
FiberSection2d::recvSelf(int commitTag, Channel &theChannel,
                         FEM_ObjectBroker &theBroker)
{
  return -1;
}

void
FiberSection2d::Print(OPS_Stream &s, int flag)
{
  const int numFibers = fibers->size();
  if (flag == OPS_PRINT_PRINTMODEL_SECTION || flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
    s << "\nFiberSection2d, tag: " << this->getTag() << "\n";
    s << "\tSection code: " << code;
    s << "\tNumber of Fibers: " << numFibers << "\n";
    s << "\tCentroid: " << yBar << "\n";
    
    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
      for (int i = 0; i < numFibers; i++) {
        s << "\nLocation (y) = (" << (*fibers)[i].y << ")";
        s << "\nArea = " << (*fibers)[i].area << "\n";
        theMaterials[i]->Print(s, flag);
      }
    }
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << OPS_PRINT_JSON_MATE_INDENT << "{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"FiberSection2d\", ";
    s << "\"computeCentroid\": " << (computeCentroid ? "true" : "false") << ", ";
    s << "\"fibers\": [\n";
    for (int i = 0; i < numFibers; i++) {
      s << OPS_PRINT_JSON_MATE_INDENT << "  ";
      s << "{\"coord\": [" << (*fibers)[i].y << ", 0.0], ";
      s << "\"area\": " << (*fibers)[i].area << ", ";
      s << "\"material\": " << theMaterials[i]->getTag();
      if (i < numFibers-1)
        s << "},\n";
      else
        s << "}\n";	
    }
    s << OPS_PRINT_JSON_MATE_INDENT << "]}";
  }
}


Response*
FiberSection2d::setResponse(const char **argv, int argc,
                            OPS_Stream &output)
{
  Response *theResponse = 0;
  
  const int numFibers = fibers->size();

  if (argc > 2 && strcmp(argv[0],"fiber") == 0) {
  
    int key = numFibers;
    int passarg = 2;
    
    if (argc <= 3) {		  // fiber number was input directly
      
      key = atoi(argv[1]);
      
    } else if (argc > 4) {  // find fiber closest to coord. with mat tag
      
      int matTag = atoi(argv[3]);
      double yCoord = atof(argv[1]);
      
      double closestDist = 0;
      double ySearch, dy;
      double distance;
      int j;
      // Find first fiber with specified material tag
      for (j = 0; j < numFibers; j++) {
        if (matTag == theMaterials[j]->getTag()) {
          //ySearch = matData[2*j];
          ySearch = (*fibers)[j].y;
          dy = ySearch-yCoord;
          closestDist = dy*dy;
          key = j;
          break;
        }
      }
      // Search the remaining fibers
      for ( ; j < numFibers; j++) {
        if (matTag == theMaterials[j]->getTag()) {
          //ySearch = matData[2*j];
          ySearch = (*fibers)[j].y;	  
          dy = ySearch-yCoord;
          distance = dy*dy;
          if (distance < closestDist) {
            closestDist = distance;
            key = j;
          }
        }
      }
      passarg = 4;
    }
    
    else { // fiber near-to coordinate specified
      
      double yCoord = atof(argv[1]);
      double closestDist;
      double ySearch, dy;
      double distance;
      
      //ySearch = matData[0];
      ySearch = (*fibers)[0].y;
      dy = ySearch-yCoord;
      closestDist = dy*dy;
      key = 0;
      for (int j = 1; j < numFibers; j++) {
        //ySearch = matData[2*j];
        ySearch = (*fibers)[j].y;
        dy = ySearch-yCoord;
        
        distance = dy*dy;
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
      output.attr("zLoc",0.0);
      output.attr("area",(*fibers)[key].area);
      
      theResponse = theMaterials[key]->setResponse(&argv[passarg], argc-passarg, output);
      
      output.endTag();
    }

  } else if (strcmp(argv[0],"fiberData") == 0) {
    
    int numData = numFibers*5;
    for (int j = 0; j < numFibers; j++) {
      output.tag("FiberOutput");
      output.attr("yLoc",(*fibers)[j].y);
      output.attr("zLoc", 0.0);
      output.attr("area",(*fibers)[j].area);    
      output.tag("ResponseType","yCoord");
      output.tag("ResponseType","zCoord");
      output.tag("ResponseType","area");
      output.tag("ResponseType","stress");
      output.tag("ResponseType","strain");
      output.endTag();
    }
    Vector theResponseData(numData);
    theResponse = new SectionResponse(*this, FiberResponse::FiberData, theResponseData);
  }

  else if (strcmp(argv[0],"fiberData2") == 0) {
    
    int numData = numFibers*6;
    for (int j = 0; j < numFibers; j++) {
      output.tag("FiberOutput");
      output.attr("yLoc",(*fibers)[j].y);
      output.attr("zLoc", 0.0);
      output.attr("area",(*fibers)[j].area);    
      output.tag("ResponseType","yCoord");
      output.tag("ResponseType","zCoord");
      output.tag("ResponseType","area");
      output.tag("ResponseType","stress");
      output.tag("ResponseType","strain");
      output.endTag();
    }
    Vector theResponseData(numData);
    theResponse = new SectionResponse(*this, FiberResponse::FiberData02, theResponseData);
  }

  else if ((strcmp(argv[0],"numFailedFiber") == 0) || (strcmp(argv[0],"numFiberFailed") == 0)) {
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
	  theResponse = new SectionResponse(*this, 8, getEnergy());
  }

  if (theResponse == 0)
    return FrameSection::setResponse(argv, argc, output);

  return theResponse;
}


int 
FiberSection2d::getResponse(int responseID, Information &sectInfo)
{
  const int numFibers = fibers->size();
  if (responseID == FiberResponse::FiberData) {
    int numData = 5*numFibers;
    Vector data(numData);
    int count = 0;
    for (int j = 0; j < numFibers; j++) {
      double yLoc = (*fibers)[j].y;
      double zLoc = 0.0;
      double A    = (*fibers)[j].area;
      double stress = theMaterials[j]->getStress();
      double strain = theMaterials[j]->getStrain();
      data(count)   = yLoc;
      data(count+1) = zLoc;
      data(count+2) = A;
      data(count+3) = stress;
      data(count+4) = strain;
      count += 5;
    }
    return sectInfo.setVector(data);	
  }

  if (responseID == FiberResponse::FiberData02) {
    int numData = 6*numFibers;
    Vector data(numData);
    int count = 0;
    for (int j = 0; j < numFibers; j++) {
      data(count)   = (*fibers)[j].y;   // y
      data(count+1) = 0.0;            // z
      data(count+2) = (*fibers)[j].area; // A
      data(count+3) = (double)theMaterials[j]->getTag();
      data(count+4) = theMaterials[j]->getStress();
      data(count+5) = theMaterials[j]->getStrain();
      count += 6;      
    }
    return sectInfo.setVector(data);

  } else  if (responseID == 6) {
    
    int count = 0;
    for (int j = 0; j < numFibers; j++) {    
      if (theMaterials[j]->hasFailed() == true)
        count++;
    }
    return sectInfo.setInt(count);

  } else  if (responseID == 7) {
    int count = 0;
    for (int j = 0; j < numFibers; j++) {    
      if (theMaterials[j]->hasFailed() == true) {
        count+=1;
      }
    }
    if (count == numFibers)
      count = 1;
    else
      count = 0;

    return sectInfo.setInt(count);
  } 
  //by SAJalali
  else if (responseID == 8) {
	  return sectInfo.setDouble(getEnergy());
  }

  return FrameSection::getResponse(responseID, sectInfo);
}



int
FiberSection2d::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  int result = -1;

  // Check if the parameter belongs to the material
  if (strstr(argv[0],"material") != 0) {
    
    if (argc < 3)
      return 0;

    // Get the tag of the material
    int materialTag = atoi(argv[1]);
    
    // Loop over fibers to find the right material
    const int numFibers = fibers->size();
    for (int i = 0; i < numFibers; i++)
      if (materialTag == theMaterials[i]->getTag()) {
        int ok = theMaterials[i]->setParameter(&argv[2], argc-2, param);
        if (ok != -1)
          result = ok;
      }
    return result;
  }

  // Check if the parameter belongs to a fiber
  // unlike setResponse, only allowing 'fiber y z matTag ...' because
  // the setResponse logic breaks down with the trailing arguments
  if (strstr(argv[0],"fiber") != 0) {
    
    const int numFibers = fibers->size();
    int key = numFibers;
    int passarg = 2;
    
    if (argc < 5)
      return 0;

    int matTag = atoi(argv[3]);
    double yCoord = atof(argv[1]);
      
    double closestDist = 0;
    double ySearch, dy;
    double distance;
    int j;
    // Find first fiber with specified material tag
    for (j = 0; j < numFibers; j++) {
      if (matTag == theMaterials[j]->getTag()) {
        ySearch = (*fibers)[j].y;
        dy = ySearch-yCoord;
        closestDist = dy*dy;
        key = j;
        break;
      }
    }
    // Search the remaining fibers
    for ( ; j < numFibers; j++) {
      if (matTag == theMaterials[j]->getTag()) {
        ySearch = (*fibers)[j].y;
        dy = ySearch-yCoord;
        distance = dy*dy;
        if (distance < closestDist) {
          closestDist = distance;
          key = j;
        }
      }
      passarg = 4;
    }
    
    // Finally, call setParameter
    if (key >= 0 && key < numFibers)
      return theMaterials[key]->setParameter(&argv[passarg], argc-passarg, param);
  }

  // Check if it belongs to the section integration
  if (strstr(argv[0],"integration") != 0) {
      return -1;
  }

  int ok = 0;

  const int numFibers = fibers->size();
  for (int i = 0; i < numFibers; i++) {
    ok = theMaterials[i]->setParameter(argv, argc, param);
    if (ok != -1)
      result = ok;
  }

  return result;
}

const Vector &
FiberSection2d::getSectionDeformationSensitivity(int gradIndex)
{
  static Vector dummy(2);

  return dummy;
}

const Vector &
FiberSection2d::getStressResultantSensitivity(int gradIndex, bool conditional)
{
  static Vector ds(2);
  
  ds.Zero();

  double tangent = 0.0;
  double sig_dAdh = 0.0;

  static double locsDeriv[10000]{};
  static double areaDeriv[10000]{};

  const int numFibers = fibers->size();
  
  for (int i = 0; i < numFibers; i++) {
    const double y = (*fibers)[i].y - yBar;
    const double A = (*fibers)[i].area;
    
    double stressGradient = A*theMaterials[i]->getStressSensitivity(gradIndex,true);

    ds(0) += stressGradient;
    ds(1) += stressGradient * -y;

    double stress = 0.0;
    if (areaDeriv[i] != 0.0 || locsDeriv[i] != 0.0)
      stress = theMaterials[i]->getStress();

    if (areaDeriv[i] != 0.0) {
      sig_dAdh = stress*areaDeriv[i];
      
      ds(0) += sig_dAdh;
      ds(1) += sig_dAdh * -y;
    }

    if (locsDeriv[i] != 0.0) {
      ds(1) += (stress*A) * -locsDeriv[i];
      
      tangent = theMaterials[i]->getTangent();
      tangent = tangent * A * e(1);
      
      ds(0) +=  -locsDeriv[i]*tangent;
      ds(1) += y*locsDeriv[i]*tangent;
    }
  }
  
  return ds;
}

const Matrix &
FiberSection2d::getInitialTangentSensitivity(int gradIndex)
{
  static Matrix dksdh(2,2);
  
  dksdh.Zero();

  double dydh, dAdh;
  double tangent = 0.0;
  double dtangentdh = 0.0;

  static double locsDeriv[10000]{};
  static double areaDeriv[10000]{};


  const int numFibers = fibers->size();

  for (int i = 0; i < numFibers; i++) {
    const double y = (*fibers)[i].y - yBar;
    const double A = (*fibers)[i].area;
    dydh = locsDeriv[i];
    dAdh = areaDeriv[i];
    
    tangent = theMaterials[i]->getInitialTangent();
    dtangentdh = theMaterials[i]->getInitialTangentSensitivity(gradIndex);

    dksdh(0,0) += dtangentdh*A + tangent*dAdh;

    dksdh(0,1) += -y*(dtangentdh*A+tangent*dAdh) - dydh*(tangent*A);

    dksdh(1,1) += 2*(y*dydh*tangent*A) + y*y*(dtangentdh*A+tangent*dAdh);
  }

  dksdh(1,0) = dksdh(0,1);

  return dksdh;
}

int
FiberSection2d::commitSensitivity(const Vector& defSens,
				  int gradIndex, int numGrads)
{
  double d0 = defSens(0);
  double d1 = defSens(1);

  dedh = defSens;

  static double locsDeriv[10000]{};
  static double areaDeriv[10000]{};

  const int numFibers = fibers->size();
  for (int i = 0; i < numFibers; i++) {
    locsDeriv[i] = 0.0;
    areaDeriv[i] = 0.0;
  }

  double kappa = e(1);

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    const double y = (*fibers)[i].y - yBar;

    // determine material strain and set it
    double strainSens = d0 - y*d1 - locsDeriv[i]*kappa;
    theMat->commitSensitivity(strainSens,gradIndex,numGrads);
  }

  return 0;
}


// by SAJalali
double FiberSection2d::getEnergy() const
{
  double energy = 0;
  const int numFibers = fibers->size();
  for (int i = 0; i < numFibers; i++) {
    const double A = (*fibers)[i].area;
    energy += A * theMaterials[i]->getEnergy();
  }
  return energy;
}
