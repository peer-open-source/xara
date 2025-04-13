//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
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
#include <SensitiveResponse.h>
typedef SensitiveResponse<FrameSection> SectionResponse;
#include <UniaxialMaterial.h>
#include <ElasticMaterial.h>

#include "FiberResponse.h"

// #include <threads/thread_pool.hpp>
// #define N_FIBER_THREADS 6

ID FrameFiberSection3d::code(4);

FrameFiberSection3d::FrameFiberSection3d(int tag, int num, UniaxialMaterial &torsion, bool compCentroid, 
                                         double mass, bool use_mass)
  : FrameSection(tag, SEC_TAG_FrameFiberSection3d, mass, use_mass),
    numFibers(0), sizeFibers(num), 
    theMaterials(nullptr), matData(new double [num*3]{}),
    QzBar(0.0), QyBar(0.0), Abar(0.0), 
    yBar(0.0), zBar(0.0), computeCentroid(compCentroid),
    theTorsion(0),
#ifdef N_FIBER_THREADS
    pool((void*)new OpenSees::thread_pool{N_FIBER_THREADS}),
#endif
    e(es), s(sr), K_wrap(ks)
{
    if (sizeFibers != 0) {
      theMaterials = new UniaxialMaterial *[sizeFibers]{};
      // matData.reset();
    }

    theTorsion = torsion.getCopy();
    if (theTorsion == nullptr)
      opserr << "FrameFiberSection3d::FrameFiberSection3d -- failed to get copy of torsion material\n";

    es.zero();
    sr.zero();
    ks.zero();

    code(0) = SECTION_RESPONSE_P;  // 0  0
    code(1) = SECTION_RESPONSE_MZ; // 1  5
    code(2) = SECTION_RESPONSE_MY; // 2  4
    code(3) = SECTION_RESPONSE_T;  // 3  3
}


// constructor for blank object that recvSelf needs to be invoked upon
FrameFiberSection3d::FrameFiberSection3d():
  FrameSection(0, SEC_TAG_FrameFiberSection3d, 0, false),
  numFibers(0), sizeFibers(0), 
  theMaterials(0), 
  matData(0),
  QzBar(0.0), QyBar(0.0), Abar(0.0), 
  yBar(0.0), zBar(0.0), computeCentroid(true),
#ifdef N_FIBER_THREADS
  pool((void*)new OpenSees::thread_pool{N_FIBER_THREADS}),
#endif
  e(es), s(sr), K_wrap(ks),
  theTorsion(nullptr)
{
  es.zero();
  sr.zero();
  ks.zero();

  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_MZ;
  code(2) = SECTION_RESPONSE_MY;
  code(3) = SECTION_RESPONSE_T;
}


FrameFiberSection3d::~FrameFiberSection3d()
{
  if (theMaterials != 0) {
    for (int i = 0; i < numFibers; i++)
      if (theMaterials[i] != nullptr)
      delete theMaterials[i];

    delete [] theMaterials;
  }

  if (theTorsion != nullptr)
    delete theTorsion;
}


double
FrameFiberSection3d::getEnergy() const
{
    double energy = 0;
    for (int i = 0; i < numFibers; i++) {
        double A = matData[3 * i + 2];
        energy += A * theMaterials[i]->getEnergy();
    }
    return energy;
}

int
FrameFiberSection3d::getIntegral(Field field, State state, double& value) const
{
  value = 0.0;

  switch (field) {
    case Field::Unit:
      for (int i=0; i<numFibers; i++) {
        const double A  = matData[3*i+2];
        value += A;
      }
      return 0;

    case Field::Density:
      // First check if density has been specified for the section
      if (this->FrameSection::getIntegral(field, state, value) == 0) 
        return 0;

      for (int i=0; i<numFibers; i++) {
        double density;
        const double A  = matData[3*i+2];
//      if (theMaterials[i]->getIntegral(field, state, density) == 0)
        if (theMaterials[i]->getRho() != 0)
          value += A*density;
        else
          return -1;
      }
      return 0;

    case Field::UnitY: // TODO: Centroid
      for (int i=0; i<numFibers; i++) {
        const double A  = matData[3*i+2];
        const double y  = matData[3*i] - yBar;
        value += A*y;
      }
      return 0;


    case Field::UnitYY:
    case Field::UnitCentroidYY:
      for (int i=0; i<numFibers; i++) {
        const double A  = matData[3*i+2];
        const double y  = matData[3*i]
                        - yBar*(Field::UnitCentroidYY == field);
        value += A*y*y;
      }
      return 0;

    case Field::UnitZZ:
    case Field::UnitCentroidZZ:
      for (int i=0; i<numFibers; i++) {
        const double A  = matData[3*i+2];
        const double z  = matData[3*i+1] 
                        - zBar*(Field::UnitCentroidZZ == field);
        value += A*z*z;
      }

    default:
      return -1;
  }
  return -1;
}


int
FrameFiberSection3d::addFiber(UniaxialMaterial &theMat, 
                              const double Area, const double yLoc, const double zLoc)
{
  // need to create a larger array
  if (numFibers == sizeFibers) {
      int newSize = 2*sizeFibers;
      if (newSize == 0) 
        newSize = 30;
      UniaxialMaterial **newArray = new UniaxialMaterial *[newSize]; 
      std::shared_ptr<double[]> newMatData(new double [3 * newSize]);

      // copy the old pointers
      for (int i = 0; i < numFibers; i++) {
        newArray[i]       = theMaterials[i];
        newMatData[3*i]   = matData[3*i];
        newMatData[3*i+1] = matData[3*i+1];
        newMatData[3*i+2] = matData[3*i+2];
      }

      // initialize new memomry
      for (int i = numFibers; i < newSize; i++) {
        newArray[i]       = nullptr;
        newMatData[3*i]   = 0.0;
        newMatData[3*i+1] = 0.0;
        newMatData[3*i+2] = 0.0;
      }
      sizeFibers = newSize;

      // set new memory
      if (theMaterials != nullptr)
        delete [] theMaterials;

      theMaterials = newArray;
      matData = newMatData;
  }
          
  // set the new pointers
  matData[numFibers*3]   = yLoc;
  matData[numFibers*3+1] = zLoc;
  matData[numFibers*3+2] = Area;
  theMaterials[numFibers] = theMat.getCopy();

  numFibers++;

  // Recompute centroid
  if (computeCentroid) {
    Abar  += Area;
    QzBar += yLoc*Area;
    QyBar += zLoc*Area;
    
    yBar = QzBar/Abar;
    zBar = QyBar/Abar;
  }
  
  return 0;
}


#ifdef N_FIBER_THREADS
#include <mutex>
int
FrameFiberSection3d::setTrialSectionDeformation(const Vector &deforms)
{
  e = deforms;
 
  sr.zero();
  ks.zero();

  const double e0 = deforms(0), // u'
               e1 = deforms(1),
               e2 = deforms(2),
               e3 = deforms(3);

  int res = 0;
  std::mutex resp_mutex;

  ((OpenSees::thread_pool*)pool)->submit_loop<unsigned int>(0, numFibers,
  [&,e0,e1,e2](int i){
    int res = 0;

    const double y  = matData[3*i]   - yBar;
    const double z  = matData[3*i+1] - zBar;
    const double A  = matData[3*i+2];

    // determine material strain and set it
    const double strain = e0 - y*e1 + z*e2;
    double tangent, stress;
    res += theMaterials[i]->setTrial(strain, stress, tangent);

    const double EA = tangent * A;

    const std::lock_guard<std::mutex> lock(resp_mutex);
    ks(0, 0) +=     EA;
    ks(0, 1) +=  -y*EA;
    ks(0, 2) +=   z*EA;

    ks(1, 1) +=  y*y*EA;
    ks(2, 2) +=  z*z*EA; 
    ks(1, 2) += -y*z*EA;

    double fs0 = stress * A;
    sr[ 0] +=    fs0;  // N
    sr[ 1] += -y*fs0;  // Mz
    sr[ 2] +=  z*fs0;  // My

    return res;
  }).wait();

  ks(1, 0) = ks(0, 1);
  ks(2, 0) = ks(0, 2);
  ks(2, 1) = ks(1, 2);
 
  if (theTorsion != nullptr) {
    double stress, tangent;
    res += theTorsion->setTrial(e3, stress, tangent);
    sr[ 3] = stress;
    ks(3, 3) = tangent;
  }

  return res;
}

#else

int
FrameFiberSection3d::setTrialSectionDeformation(const Vector &deforms)
{
  e = deforms;
 
  sr.zero();
  ks.zero();

  const double e0 = deforms(0), // u'
               k1 = deforms(1),
               k2 = deforms(2),
               e3 = deforms(3);

  int res = 0;
  for (int i = 0; i < numFibers; i++) {

    const double y  = matData[3*i]   - yBar;
    const double z  = matData[3*i+1] - zBar;
    const double A  = matData[3*i+2];

    // Determine material strain and set it
    double strain = e0 - y*k1 + z*k2;
    double tangent, stress;
    res += theMaterials[i]->setTrial(strain, stress, tangent);

    double EA     = tangent * A;

    ks( 0, 0) +=     EA;
    ks( 0, 1) +=  -y*EA;
    ks( 0, 2) +=   z*EA;
//  ks( 0, 2) +=   z*EA;

    ks( 1, 1) +=  y*y*EA; // 5
    ks( 2, 2) +=  z*z*EA; // 10
    ks( 1, 2) += -y*z*EA;

    double fs0 = stress * A;
    sr[0] +=    fs0;  // N
    sr[1] += -y*fs0;  // Mz
    sr[2] +=  z*fs0;  // My
  }

  ks(1, 0) = ks(0, 1);
  ks(2, 0) = ks(0, 2);
  ks(2, 1) = ks(1, 2);
 
  if (theTorsion != nullptr) {
    double stress, tangent;
    res += theTorsion->setTrial(e3, stress, tangent);
    sr[ 3] = stress;
    ks(3, 3) = tangent;
  }

  return res;
}
#endif



const Matrix&
FrameFiberSection3d::getInitialTangent()
{
  static double kInitialData[16];
  static Matrix kInitial(kInitialData, 4, 4);
  
  kInitial.Zero();

  for (int i = 0; i < numFibers; i++) {
    const double y = matData[3*i]   - yBar;
    const double z = matData[3*i+1] - zBar;
    const double A = matData[3*i+2];

    double tangent = theMaterials[i]->getInitialTangent();

    const double EA = tangent * A;
    const double vas1    = -y*EA;
    const double vas2    =  z*EA;
    const double vas1as2 =  vas1*z;

    kInitialData[0] +=      EA;
    kInitialData[1] +=   -y*EA;
    kInitialData[2] +=    z*EA;

    kInitialData[5] +=  y*y*EA;
    kInitialData[6] += -y*z*EA;

    kInitialData[10] += vas2 * z; 
  }

  kInitialData[4] = kInitialData[1];
  kInitialData[8] = kInitialData[2];
  kInitialData[9] = kInitialData[6];

  if (theTorsion != nullptr)
    kInitialData[15] = theTorsion->getInitialTangent();

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
  FrameFiberSection3d *theCopy = new FrameFiberSection3d();
  theCopy->setTag(this->getTag());
  theCopy->numFibers  = numFibers;
  theCopy->sizeFibers = numFibers;
  theCopy->pool       = pool;

  if (numFibers != 0) {
    theCopy->theMaterials = new UniaxialMaterial *[numFibers];

    theCopy->matData = matData; // new double [numFibers*3];

    for (int i = 0; i < numFibers; i++) {
      theCopy->theMaterials[i] = theMaterials[i]->getCopy();

      if (theCopy->theMaterials[i] == nullptr) {
        delete theCopy;
        return nullptr;
      }
    }    
  }

  theCopy->e = e;
  theCopy->sr = sr;
  theCopy->ks = ks;
  theCopy->QzBar = QzBar;
  theCopy->QyBar = QyBar;
  theCopy->Abar  = Abar;
  theCopy->yBar  = yBar;
  theCopy->zBar  = zBar;
  theCopy->computeCentroid = computeCentroid;

  if (theTorsion != nullptr)
    theCopy->theTorsion = theTorsion->getCopy();
  else
    theCopy->theTorsion = nullptr;

  return theCopy;
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
  for (int i = 0; i < numFibers; i++)
    err += theMaterials[i]->commitState();

  if (theTorsion != nullptr)
    err += theTorsion->commitState();

  return err;
}

int
FrameFiberSection3d::revertToLastCommit()
{
  int err = 0;

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    // invoke revertToLast on the material
    err += theMat->revertToLastCommit();
  }


  if (theTorsion != nullptr)
    err += theTorsion->revertToLastCommit();

  return err;
}

int
FrameFiberSection3d::revertToStart()
{
  // revert the fibers to start    
  int err = 0;

  for (int i = 0; i < numFibers; i++) {
    UniaxialMaterial *theMat = theMaterials[i];
    // invoke revertToStart on the material
    err += theMat->revertToStart();
  }

  if (theTorsion != nullptr)
    err += theTorsion->revertToStart();

  return err;
}

int
FrameFiberSection3d::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;

  // create an id to send objects tag and numFibers, 
  // size 5 so no conflict with matData below if just 2 fibers
  static ID data(5);
  data(0) = this->getTag();
  data(1) = numFibers;
  data(2) = (theTorsion != 0) ? 1 : 0;
  int dbTag = this->getDbTag();
  if (theTorsion != nullptr) {
    theTorsion->setDbTag(dbTag);
    data(3) = theTorsion->getClassTag();
  }
  data(4) = computeCentroid ? 1 : 0; // Now the ID data is really 5

  res += theChannel.sendID(dbTag, commitTag, data);
  if (res < 0) {
    opserr << "FrameFiberSection3d::sendSelf - failed to send ID data\n";
    return res;
  }    

  if (theTorsion != nullptr)
    theTorsion->sendSelf(commitTag, theChannel);

  if (numFibers != 0) { 
    // create an id containingg classTag and dbTag for each material & send it
    ID materialData(2*numFibers);
    for (int i=0; i<numFibers; i++) {
      UniaxialMaterial *theMat = theMaterials[i];
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
     opserr << "FrameFiberSection3d::sendSelf - failed to send material data\n";
     return res;
    }    

    // send the fiber data, i.e. area and loc
    Vector fiberData(matData, 3*numFibers);
    res += theChannel.sendVector(dbTag, commitTag, fiberData);
    if (res < 0) {
     opserr << "FrameFiberSection3d::sendSelf - failed to send fiber data\n";
     return res;
    }    

    // now invoke send on all the materials
    for (int j=0; j<numFibers; j++)
      theMaterials[j]->sendSelf(commitTag, theChannel);
  }

  return res;
}

int
FrameFiberSection3d::recvSelf(int commitTag, Channel &theChannel,
                   FEM_ObjectBroker &theBroker)
{
  int res = 0;

  static ID data(5);

  int dbTag = this->getDbTag();
  res += theChannel.recvID(dbTag, commitTag, data);

  if (res < 0) {
   opserr << "FrameFiberSection3d::recvSelf - failed to recv ID data\n";
   return res;
  } 

  this->setTag(data(0));

  if (data(2) == 1 && theTorsion == 0) {      
    int cTag = data(3);
    theTorsion = theBroker.getNewUniaxialMaterial(cTag);
    if (theTorsion == 0) {
      opserr << "FrameFiberSection3d::recvSelf - failed to get torsion material \n";
      return -1;
    }
    theTorsion->setDbTag(dbTag);
  }

  if (theTorsion->recvSelf(commitTag, theChannel, theBroker) < 0) {
         opserr << "FrameFiberSection3d::recvSelf - torsion failed to recvSelf \n";
       return -2;
  }
  
  // recv data about materials objects, classTag and dbTag
  if (data(1) != 0) {
    ID materialData(2*data(1));
    res += theChannel.recvID(dbTag, commitTag, materialData);
    if (res < 0) {
     opserr << "FrameFiberSection3d::recvSelf - failed to recv material data\n";
     return res;
    }    

    // if current arrays not of correct size, release old and resize
    if (theMaterials == 0 || numFibers != data(1)) {
      // delete old stuff if outa date
      if (theMaterials != 0) {
        for (int i=0; i<numFibers; i++)
          delete theMaterials[i];
        delete [] theMaterials;
  //    if (matData != 0)
  //      delete [] matData;
  //    matData = 0;
        theMaterials = 0;
      }

      // create memory to hold material pointers and fiber data
      numFibers = data(1);
      sizeFibers = data(1);
      if (numFibers != 0) {
        theMaterials = new UniaxialMaterial *[numFibers];
        
        for (int j=0; j<numFibers; j++)
          theMaterials[j] = 0;
        
//      matData = new double [numFibers*3];
        matData.reset(new double [numFibers*3]);
      }
    }

    Vector fiberData(matData, 3*numFibers);
    res += theChannel.recvVector(dbTag, commitTag, fiberData);
    if (res < 0) {
     opserr << "FrameFiberSection3d::recvSelf - failed to recv fiber data\n";
     return res;
    }    
    
    for (int i=0; i<numFibers; i++) {
      int classTag = materialData(2*i);
      int dbTag = materialData(2*i+1);

      // if material pointed to is blank or not of corrcet type, 
      // release old and create a new one
      if (theMaterials[i] == 0)
        theMaterials[i] = theBroker.getNewUniaxialMaterial(classTag);
        else if (theMaterials[i]->getClassTag() != classTag) {
        delete theMaterials[i];
        theMaterials[i] = theBroker.getNewUniaxialMaterial(classTag);      
      }

      theMaterials[i]->setDbTag(dbTag);
      res += theMaterials[i]->recvSelf(commitTag, theChannel, theBroker);
    }

    QzBar = 0.0;
    QyBar = 0.0;
    Abar  = 0.0;
    double yLoc, zLoc, Area;

    computeCentroid = data(4) ? true : false;
    
    // Recompute centroid
    for (int i = 0; computeCentroid && i < numFibers; i++) {
      yLoc = matData[3*i];
      zLoc = matData[3*i+1];
      Area = matData[3*i+2];
      Abar  += Area;
      QzBar += yLoc*Area;
      QyBar += zLoc*Area;
    }

    if (computeCentroid) {
      yBar = QzBar/Abar;
      zBar = QyBar/Abar;
    } else {
      yBar = 0.0;
      zBar = 0.0;      
    }
  }    

  return res;
}

void
FrameFiberSection3d::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << OPS_PRINT_JSON_MATE_INDENT << "{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"" << this->getClassType() << "\", ";
        if (theTorsion != 0)
          s << "\"torsion\": " << theTorsion->getInitialTangent() << ", ";

        double mass;
        if (this->FrameSection::getIntegral(Field::Density, State::Init, mass) == 0)
          s << "\"mass\": " << mass;

        s << "\"fibers\": [\n";
        for (int i = 0; i < numFibers; i++) {
              s << OPS_PRINT_JSON_MATE_INDENT << "\t{\"coord\": [" << matData[3*i] << ", " << matData[3*i+1] << "], ";
              s << "\"area\": " << matData[3*i+2] << ", ";
              s << "\"material\": " << theMaterials[i]->getTag();
              if (i < numFibers - 1)
                  s << "},\n";
              else
                  s << "}\n";
        }
        s << OPS_PRINT_JSON_MATE_INDENT << "]}";
        return;
  }

  if (flag == OPS_PRINT_PRINTMODEL_SECTION || flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
    s << "\nFrameFiberSection3d, tag: " << this->getTag() << endln;
    s << "\tSection code: " << code;
    s << "\tNumber of Fibers: " << numFibers << endln;
    s << "\tCentroid: (" << yBar << ", " << zBar << ')' << endln;
    if (theTorsion != 0)
        theTorsion->Print(s, flag);    

    if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
      for (int i = 0; i < numFibers; i++) {
      s << "\nLocation (y, z) = (" << matData[3*i] << ", " << matData[3*i+1] << ")";
      s << "\nArea = " << matData[3*i+2] << endln;
      theMaterials[i]->Print(s, flag);
      
      }
    }
  }
  if (flag == 3) {
    for (int i = 0; i < numFibers; i++) {
      s << theMaterials[i]->getTag() << " " << matData[3*i] << " "  << matData[3*i+1] << " "  << matData[3*i+2] << " " ;
      s << theMaterials[i]->getStress() << " "  << theMaterials[i]->getStrain() << endln;
    } 
  }
    
  if (flag == 4) {
    for (int i = 0; i < numFibers; i++) {
      s << "add fiber # " << i+1 << " using material # " << theMaterials[i]->getTag() << " to section # 1\n";
      s << "fiber_cross_section = " << matData[3*i+2] << "*m^2\n";
      s << "fiber_location = (" << matData[3*i] << "*m, " << matData[3*i+1] << "*m);\n\n";
    }
  }
}

Response*
FrameFiberSection3d::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = nullptr;
  
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
        if (matTag == theMaterials[j]->getTag()) {
          ySearch = matData[3*j];
          zSearch = matData[3*j+1];
          dy = ySearch-yCoord;
          dz = zSearch-zCoord;
          closestDist = dy*dy + dz*dz;
          key = j;
          break;
        }
      }

      // Search the remaining fibers
      for ( ; j < numFibers; j++) {
        if (matTag == theMaterials[j]->getTag()) {
          ySearch = matData[3*j];
          zSearch = matData[3*j+1];

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
      ySearch = matData[0];
      zSearch = matData[1];
      // ySearch = yLocs[0];
      // zSearch = zLocs[0];
      dy = ySearch-yCoord;
      dz = zSearch-zCoord;
      closestDist = sqrt(dy*dy + dz*dz);
      key = 0;
      for (int j = 1; j < numFibers; j++) {
        ySearch = matData[3*j];
        zSearch = matData[3*j+1];                      
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
      output.attr("yLoc",matData[3*key]);
      output.attr("zLoc",matData[3*key+1]);
      output.attr("area",matData[3*key+2]);
      
      theResponse = theMaterials[key]->setResponse(&argv[passarg], argc-passarg, output);
      
      output.endTag();
    }
  
  } else if (strcmp(argv[0],"fiberData") == 0) {
    int numData = numFibers*5;
    for (int j = 0; j < numFibers; j++) {
      output.tag("FiberOutput");
      output.attr("yLoc", matData[3*j]);
      output.attr("zLoc", matData[3*j+1]);
      output.attr("area", matData[3*j+2]);    
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
      output.attr("yLoc", matData[3*j]);
      output.attr("zLoc", matData[3*j+1]);
      output.attr("area", matData[3*j+2]);    
      output.attr("material", theMaterials[j]->getTag());
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
  // Just call the base class method ... don't need to define
  // this function, but keeping it here just for clarity
  if (responseID == FiberResponse::FiberData) {
    int numData = 5*numFibers;
    Vector data(numData);
    int count = 0;
    for (int j = 0; j < numFibers; j++) {
      double yLoc = matData[3*j];
      double zLoc = matData[3*j+1];
      double A = matData[3*j+2];
      double stress = theMaterials[j]->getStress();
      double strain = theMaterials[j]->getStrain();
      data(count) = yLoc; data(count+1) = zLoc; data(count+2) = A;
      data(count+3) = stress; data(count+4) = strain;
      count += 5;
    }
    return sectInfo.setVector(data);

  } else if (responseID == FiberResponse::FiberData02) {
    int numData = 6*numFibers;
    Vector data(numData);
    int count = 0;
    for (int j = 0; j < numFibers; j++) {
      data(count)   = matData[3*j  ]; // y
      data(count+1) = matData[3*j+1]; // z
      data(count+2) = matData[3*j+2]; // A
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

  // A material parameter
  if (strstr(argv[0],"material") != 0) {

    // Get the tag of the material
    int paramMatTag = atoi(argv[1]);

    // Loop over fibers to find the right material(s)
    int ok = 0;
    for (int i = 0; i < numFibers; i++)
      if (paramMatTag == theMaterials[i]->getTag()) {
      ok = theMaterials[i]->setParameter(&argv[2], argc-2, param);
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
  else if (strstr(argv[0],"integration") != 0)
      return -1;

  int ok = 0;
  
  // loop over every material
  for (int i = 0; i < numFibers; i++) {
    ok = theMaterials[i]->setParameter(argv, argc, param);
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
FrameFiberSection3d::getSectionDeformationSensitivity(int gradIndex)
{
  static Vector dummy(4);
  
  dummy.Zero();
  
  return dummy;
}

   
const Vector &
FrameFiberSection3d::getStressResultantSensitivity(int gradIndex, bool conditional)
{
  static Vector ds(4);
  
  ds.Zero();
  
  double stress = 0;
  double dsigdh = 0;
  double sig_dAdh = 0;
  double tangent = 0;
  static double dydh[10000];
  static double dzdh[10000];
  static double areaDeriv[10000];
#if 0
  if (sectionIntegr != 0) {
    sectionIntegr->getLocationsDeriv(numFibers, dydh, dzdh);  
    sectionIntegr->getWeightsDeriv(numFibers, areaDeriv);
  } else
#endif
    for (int i = 0; i < numFibers; i++) {
      dydh[i] = 0.0;
      dzdh[i] = 0.0;
      areaDeriv[i] = 0.0;
    }
  
  for (int i = 0; i < numFibers; i++) {

    double y  = matData[3*i]   - yBar;
    double z  = matData[3*i+1] - zBar;
    double A  = matData[3*i+2];
    
    dsigdh = theMaterials[i]->getStressSensitivity(gradIndex, conditional);

    ds(0) += dsigdh*A;
    ds(1) += -y*dsigdh*A;
    ds(2) +=  z*dsigdh*A;

    if (areaDeriv[i] != 0.0 || dydh[i] != 0.0 ||  dzdh[i] != 0.0)
      stress = theMaterials[i]->getStress();

    if (dydh[i] != 0.0 || dzdh[i] != 0.0)
      tangent = theMaterials[i]->getTangent();

    if (areaDeriv[i] != 0.0) {
      sig_dAdh = stress*areaDeriv[i];
      
      ds(0) += sig_dAdh;
      ds(1) += -y*sig_dAdh;
      ds(2) +=  z*sig_dAdh;
    }

    if (dydh[i] != 0.0)
      ds(1) += -dydh[i] * (stress*A);

    if (dzdh[i] != 0.0)
      ds(2) +=  dzdh[i] * (stress*A);

    static Matrix as(1,3);
    as(0,0) = 1;
    as(0,1) = -y;
    as(0,2) = z;
    
    static Matrix dasdh(1,3);
    dasdh(0,1) = -dydh[i];
    dasdh(0,2) = dzdh[i];
    
    static Matrix tmpMatrix(3,3);
    tmpMatrix.addMatrixTransposeProduct(0.0, as, dasdh, tangent);
    
    //ds.addMatrixVector(1.0, tmpMatrix, e, A);
    ds(0) += (tmpMatrix(0,0)*e(0) + tmpMatrix(0,1)*e(1) + tmpMatrix(0,2)*e(2))*A;
    ds(1) += (tmpMatrix(1,0)*e(0) + tmpMatrix(1,1)*e(1) + tmpMatrix(1,2)*e(2))*A;
    ds(2) += (tmpMatrix(2,0)*e(0) + tmpMatrix(2,1)*e(1) + tmpMatrix(2,2)*e(2))*A;
  }

  ds(3) = theTorsion->getStressSensitivity(gradIndex, conditional);

  return ds;
}

const Matrix &
FrameFiberSection3d::getSectionTangentSensitivity(int gradIndex)
{
  static Matrix something(nsr,nsr);
  
  something.Zero();

  something(3,3) = theTorsion->getTangentSensitivity(gradIndex);
  
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

  static double yLocs[10000];
  static double zLocs[10000];

  { // TODO
    for (int i = 0; i < numFibers; i++) {
      yLocs[i] = matData[3*i];
      zLocs[i] = matData[3*i+1];
    }
  }

  static double dydh[10000];
  static double dzdh[10000];

  { // TODO
    for (int i = 0; i < numFibers; i++) {
      dydh[i] = 0.0;
      dzdh[i] = 0.0;
    }
  }


  double depsdh = 0;

  for (int i = 0; i < numFibers; i++) {
    double y = yLocs[i] - yBar;
    double z = zLocs[i] - zBar;

    // determine material strain and set it
    depsdh = d0 - y*d1 + z*d2 - dydh[i]*e(1) + dzdh[i]*e(2);

    theMaterials[i]->commitSensitivity(depsdh,gradIndex,numGrads);
  }

  theTorsion->commitSensitivity(d3, gradIndex, numGrads);

  return 0;
}

