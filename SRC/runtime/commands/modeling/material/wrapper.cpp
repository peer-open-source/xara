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
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
//
#include <tcl.h>
#include <set>
#include <vector>
#include <Logging.h>
#include <Parsing.h>
#include <ArgumentTracker.h>

#include <ModelRegistry.h>
#include <NDMaterial.h>
#include <UniaxialMaterial.h>

#include <InitStressMaterial.h>
#include <InitStrainMaterial.h>

#include <InitStrainNDMaterial.h>
#include <InitStressNDMaterial.h>
#include <ContinuumUniaxial.h>

#include <PlaneStressMaterial.h>
#include <PlaneStrainMaterial.h>
#include <PlaneStressRebarMaterial.h>

#include <ParallelMaterial.h>
#include <FatigueMaterial.h>

#include <PlateRebarMaterial.h>
#include <PlateFiberMaterial.h>

#include <PlateFiberMaterial.h>

#include <BeamFiberMaterial.h>
#include <BeamFiberMaterial2d.h>
#include <BeamFiberMaterial2dPS.h>

#include <Parallel3DMaterial.h>
#include <OrthotropicMaterial.h>

#ifdef _MSC_VER 
#  include <string.h>
#  define strcasecmp _stricmp
#else
#  include <strings.h>
#endif

int
TclCommand_addWrappingMaterial(ClientData clientData, Tcl_Interp* interp,
                               int argc, TCL_Char** const argv)
{
  //
  //
  //
  ModelRegistry *builder = static_cast<ModelRegistry*>(clientData);

  if (argc < 4) {
    opserr << OpenSees::PromptValueError << " insufficient arguments\n";
    return TCL_ERROR;
  }

  int tago, tagi;
  if (Tcl_GetInt(interp, argv[2], &tago) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "failed to read tag\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3], &tagi) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "failed to read tag\n";
    return TCL_ERROR;
  }

  if (strcmp(argv[1], "ContinuumWrapper") == 0 || strcmp(argv[1], "Continuum") == 0) {
    NDMaterial* inside = builder->getTypedObject<NDMaterial>(tagi);
    if (inside == nullptr) 
      return TCL_ERROR;
    NDMaterial*test = inside->getCopy("ThreeDimensional");
    if (test == nullptr) {
      delete test;
      opserr << OpenSees::PromptValueError << "ContinuumUniaxial only works with 3D materials\n";
      return TCL_ERROR;
    }
    delete test;
    return builder->addTypedObject<UniaxialMaterial>(tago, new ContinuumUniaxial(tago, *inside));
  }

  else if (strcmp(argv[0], "uniaxialMaterial") == 0 && (
            strstr(argv[1], "InitStress") != 0 || 
            strstr(argv[1], "InitialStress") != 0 || 
            strstr(argv[1], "InitialStrain") != 0 || 
            strstr(argv[1], "InitStrain") != 0)) {
          
      double initial;
      if (Tcl_GetDouble(interp, argv[4], &initial) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "failed to read initial value\n";
          return TCL_ERROR;
      }
      UniaxialMaterial* inside = builder->getTypedObject<UniaxialMaterial>(tagi);
      if (inside == nullptr) 
          return TCL_ERROR;

      if (strstr(argv[1], "Stress") != 0)
          return builder->addTypedObject<UniaxialMaterial>(tago, new InitStressMaterial(tago, *inside, initial));
      else if (strstr(argv[1], "Strain") != 0)
          return builder->addTypedObject<UniaxialMaterial>(tago, new InitStrainMaterial(tago, *inside, initial));
  }

  else if (strcmp(argv[0], "nDMaterial") == 0 && (
            strstr(argv[1], "InitStress") != 0 || 
            strstr(argv[1], "InitialStress") != 0 || 
            strstr(argv[1], "InitialStrain") != 0 || 
            strstr(argv[1], "InitStrain") != 0)) {
          
      Vector initial(6);
      NDMaterial* inside = builder->getTypedObject<NDMaterial>(tagi);
      if (inside == nullptr) 
        return TCL_ERROR;
      
      inside = inside->getCopy("ThreeDimensional");
      if (!inside || strcmp(inside->getType(), "ThreeDimensional") != 0) {
        opserr << OpenSees::PromptValueError 
                << "InitStressNDMaterial only works with 3D materials\n";
        return TCL_ERROR;
      }

      if (argc == 5) {
        double evol;
        if (Tcl_GetDouble(interp, argv[4], &evol) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "failed to read initial value\n";
          return TCL_ERROR;
        }
        for (int i = 0; i < 3; ++i)
          initial(i) = evol;

      } else {
          for (int i=0; i<6; ++i) {
            if (Tcl_GetDouble(interp, argv[4+i], &initial(i)) != TCL_OK) {
              opserr << OpenSees::PromptValueError << "failed to read initial value\n";
              return TCL_ERROR;
            }
          }
      }

      NDMaterial* wrapper = nullptr;
      if (strstr(argv[1], "Strain") != 0)
        wrapper = new InitStrainNDMaterial(tago, *inside, initial);
      
      int status = builder->addTypedObject<NDMaterial>(tago, wrapper);

      delete inside;
      return status;
  }


  return TCL_ERROR;
}


int
TclCommand_newParallelMaterial(ClientData clientData, Tcl_Interp* interp, int argc, G3_Char ** const argv)
{
  ModelRegistry* builder = static_cast<ModelRegistry*>(clientData);
  assert(builder != nullptr);

  if (argc < 4) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial Parallel tag? tag1? tag2? ...";
    opserr << " <-min min?> <-max max?>" << "\n";
    return TCL_ERROR;
  }

  int tag;
  UniaxialMaterial* theMaterial = nullptr;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial Parallel tag" << "\n";
    return TCL_ERROR;
  }

  int numMaterials = argc-3;
  
  if (numMaterials == 0) {
    opserr << "WARNING no component material(s) provided\n";
    return TCL_ERROR;
  }

  // Create an array to hold pointers to component materials
  UniaxialMaterial **theMats = new UniaxialMaterial *[numMaterials];
  
  // For each material get the tag and ensure it exists in model already
  for (int i=0; i<numMaterials; i++) {
    int tagI;
    if (Tcl_GetInt(interp, argv[i+3], &tagI) != TCL_OK) {
      opserr << "WARNING invalid component tag " << argv[i+3] << "\n";
      return TCL_ERROR;
    }

    UniaxialMaterial *theMat = builder->getTypedObject<UniaxialMaterial>(tagI);
    if (theMat == nullptr) {
      delete [] theMats;
      return TCL_ERROR;
    }
    theMats[i] = theMat;
  }
  
  // Parsing was successful, allocate the material
  theMaterial = new ParallelMaterial(tag, numMaterials, theMats);
  builder->addTaggedObject<UniaxialMaterial>(*theMaterial);
  
  delete [] theMats;
  return TCL_OK;
}


// material type? tag? nd_tag?
int
TclCommand_newPlateFiber(ClientData clientData, Tcl_Interp* interp, int argc, G3_Char ** const argv)
{
  ModelRegistry *builder = static_cast<ModelRegistry*>(clientData);

  if (argc < 4) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: nDMaterial " << argv[1] << " tag? matTag?" << "\n";
    return TCL_ERROR;
  }

  int tag, matTag;
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid nDMaterial tag" << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
    opserr << "WARNING invalid matTag " << argv[3] << "\n";
    return TCL_ERROR;
  }

  NDMaterial *threeDMaterial = builder->getTypedObject<NDMaterial>(matTag);
  if (threeDMaterial == nullptr)
    return TCL_ERROR;
  
  // Check if the material is a 3D material
  threeDMaterial = threeDMaterial->getCopy("ThreeDimensional");
  if (threeDMaterial == nullptr) {
    opserr << "a ThreeDimensional material was expected\n";
    return TCL_ERROR;
  }

  NDMaterial *theMaterial = nullptr;
  if (strcmp(argv[1], "PlateFiberMaterial") == 0 ||
      strcmp(argv[1], "PlateFiber") == 0) {
    theMaterial = new PlateFiberMaterial(tag, *threeDMaterial);
  }
  // else if (strcmp(argv[1], "PlateFiberMaterialThermal") == 0 ||
  //     strcmp(argv[1], "PlateFiberThermal") == 0) {
  //   theMaterial = new PlateFiberMaterialThermal(tag, *threeDMaterial);
  // }
  else if (strcmp(argv[1], "BeamFiberMaterial") == 0 ||
           strcmp(argv[1], "BeamFiber") == 0) {
    theMaterial = new BeamFiberMaterial(tag, *threeDMaterial);
  }
  else if (strcmp(argv[1], "BeamFiberMaterial2d") == 0 ||
           strcmp(argv[1], "BeamFiber2d") == 0) {
    theMaterial = new BeamFiberMaterial2d(tag, *threeDMaterial);
  }
  else if (strcmp(argv[1], "BeamFiberMaterial2dPS") == 0 ||
           strcmp(argv[1], "BeamFiber2dPS") == 0) {
    theMaterial = new BeamFiberMaterial2dPS(tag, *threeDMaterial);
  }

  else {
    opserr << "WARNING invalid material type " << "\n";
    return TCL_ERROR;
  }

  if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK) {
    delete theMaterial;
    delete threeDMaterial;
    return TCL_ERROR;
  }

  delete threeDMaterial;
  return TCL_OK;
}


int 
TclCommand_newPlateRebar(ClientData clientData, Tcl_Interp* interp, int argc, G3_Char ** const argv)
{
  //
  // nDMaterial type? tag? uni_tag? angle?
  //
  ModelRegistry *builder = static_cast<ModelRegistry*>(clientData);
  if (argc < 5) {
    opserr << "WARNING insufficient arguments\n";
    return TCL_ERROR;
  }

  int tag, matTag;
  double angle;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid nDMaterial PlateRebar tag" << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
    opserr << "WARNING invalid matTag " << argv[3] << "\n";
    return TCL_ERROR;
  }

  UniaxialMaterial *theMat = builder->getTypedObject<UniaxialMaterial>(matTag);
  if (theMat == nullptr)
    return TCL_ERROR;

  if (Tcl_GetDouble(interp, argv[4], &angle) != TCL_OK) {
    opserr << "WARNING invalid angle" << "\n";
    return TCL_ERROR;
  }

  if (strcmp(argv[1], "PlateRebarMaterial") == 0 ||
      strcmp(argv[1], "PlateRebar") == 0) {
    if (builder->addTaggedObject<NDMaterial>(*new PlateRebarMaterial(tag, *theMat, angle)) != TCL_OK) {
      return TCL_ERROR;
    }
  }

  else if (strcmp(argv[1], "PlaneStressRebarMaterial") == 0 ||
           strcmp(argv[1], "PlaneStressRebar") == 0) {
    if (builder->addTaggedObject<NDMaterial>(*new PlaneStressRebarMaterial(tag, *theMat, angle)) != TCL_OK) {
      return TCL_ERROR;
    }
  }

  else {
    opserr << "WARNING invalid material type " << "\n";
    return TCL_ERROR;
  }
  return TCL_OK;
}



int
TclCommand_newFatigueMaterial(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  ModelRegistry *builder = static_cast<ModelRegistry*>(clientData);

  if (argc < 4) {
    opserr << OpenSees::PromptValueError << "insufficient arguments\n";
    opserr << "Want: uniaxialMaterial Fatigue tag? matTag?";
    opserr << " <-D_max dmax?> <-e0 e0?> <-m m?>" << "\n";
    opserr << " <-min min?> <-max max?>" << "\n";
    return TCL_ERROR;
  }

  int tag, matTag;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid uniaxialMaterial Fatigue tag" << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid component tag\n";
    return TCL_ERROR;
  }

  double Dmax = 1.0;
  double E0 = 0.191;
  double m = -0.458;
  double epsmin = NEG_INF_STRAIN;
  double epsmax = POS_INF_STRAIN;

  for (int j = 4; j < argc; j++) {
    if (strcmp(argv[j], "-Dmax") == 0) {
      if ((j + 1 >= argc) ||
          (Tcl_GetDouble(interp, argv[++j], &Dmax) != TCL_OK)) {
        opserr << OpenSees::PromptValueError << "invalid -Dmax";
        return TCL_ERROR;
      }
    } else if (strcmp(argv[j], "-E0") == 0) {
      if ((j + 1 >= argc) || (Tcl_GetDouble(interp, argv[++j], &E0) != TCL_OK)) {
        opserr << OpenSees::PromptValueError << "invalid -E0";
        return TCL_ERROR;
      }
    } else if (strcmp(argv[j], "-m") == 0) {
      if ((j + 1 >= argc) ||
          (Tcl_GetDouble(interp, argv[++j], &m) != TCL_OK)) {
        opserr << OpenSees::PromptValueError << "invalid -m";
        return TCL_ERROR;
      }
    } else if (strcmp(argv[j], "-min") == 0) {
      if ((j + 1 >= argc) ||
          (Tcl_GetDouble(interp, argv[++j], &epsmin) != TCL_OK)) {
        opserr << OpenSees::PromptValueError << "invalid -min ";
        return TCL_ERROR;
      }
    } else if (strcmp(argv[j], "-max") == 0) {
      if ((j + 1 >= argc) ||
          (Tcl_GetDouble(interp, argv[++j], &epsmax) != TCL_OK)) {
        opserr << OpenSees::PromptValueError << "invalid -max";
        return TCL_ERROR;
      }
    }
  }

  UniaxialMaterial *theMat = builder->getTypedObject<UniaxialMaterial>(matTag);

  if (theMat == nullptr) {
    opserr << OpenSees::PromptValueError << "component material does not exist\n";
    return TCL_ERROR;
  }

  // Parsing was successful, allocate the material
  UniaxialMaterial *theMaterial =
      new FatigueMaterial(tag, *theMat, Dmax, E0, m, epsmin, epsmax);

  if (builder->addTaggedObject<UniaxialMaterial>(*theMaterial) != TCL_OK) {
    delete theMaterial;
    return TCL_ERROR;
  }
  return TCL_OK;
}


int
TclCommand_addPlaneWrapper(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{

  ModelRegistry *builder = static_cast<ModelRegistry*>(clientData);

  NDMaterial * theMaterial = nullptr;

  if (strcmp(argv[1], "PlaneStressMaterial") == 0 ||
           strcmp(argv[1], "PlaneStress") == 0) {
    if (argc < 4) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial PlaneStress tag? matTag?" << "\n";
      return TCL_ERROR;
    }

    int tag, matTag;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid nDMaterial PlaneStress tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
      opserr << "WARNING invalid matTag" << "\n";
      return TCL_ERROR;
    }

    NDMaterial *threeDMaterial = builder->getTypedObject<NDMaterial>(matTag);
    if (threeDMaterial == nullptr)
      return TCL_ERROR;

    theMaterial = new PlaneStressMaterial(tag, *threeDMaterial);
  }

  // PlaneStrainMaterial
  else if (strcmp(argv[1], "PlaneStrainMaterial") == 0 ||
           strcmp(argv[1], "PlaneStrain") == 0) {
    if (argc < 4) {
      opserr << "WARNING insufficient arguments\n";
      return TCL_ERROR;
    }

    int tag, matTag;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid nDMaterial tag " << argv[2] << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
      opserr << "WARNING invalid matTag " << argv[3] << "\n";
      return TCL_ERROR;
    }

    NDMaterial *threeDMaterial = builder->getTypedObject<NDMaterial>(matTag);
    if (threeDMaterial == nullptr)
      return TCL_ERROR;

    theMaterial = new PlaneStrainMaterial(tag, *threeDMaterial);
  }

  // Done parsing
  if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK ) {
    delete theMaterial;
    return TCL_ERROR;
  }

  return TCL_OK;

}


int
TclCommand_addParallel3DMaterial(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  ModelRegistry *builder = static_cast<ModelRegistry*>(clientData);

  // Accept:
  //   nDMaterial Parallel3D tag? matTag1? matTag2? ... <-weights w1? w2? ...>?
  // or
  //   nDMaterial Parallel3D tag? {matTag1? matTag2? ...} <-weights {w1? w2? ...}>?

  const char* info =
    "nDMaterial Parallel3D tag? matTag1? matTag2? ... <-weights w1? w2? ...>?";

  if (argc < 4) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: " << info << "\n";
    return TCL_ERROR;
  }

  int tag;
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid nDMaterial Parallel3D tag " << argv[2] << "\n";
    return TCL_ERROR;
  }

  int weightLoc = -1;
  for (int i = 3; i < argc; i++) {
    if (strcmp(argv[i], "-weights") == 0) {
      if (weightLoc >= 0) {
        opserr << "WARNING nDMaterial Parallel3D cannot use -weights more than once\n";
        return TCL_ERROR;
      }
      weightLoc = i;
    }
  }

  int matStart = 3;
  int matEnd   = (weightLoc >= 0) ? weightLoc : argc;

  if (matEnd <= matStart) {
    opserr << "WARNING nDMaterial Parallel3D requires at least one component material\n";
    return TCL_ERROR;
  }

  std::vector<NDMaterial*> materials;
  std::vector<double> weights;

  bool freeMatArgv = false;
  bool freeWeightArgv = false;

  TCL_Char** matArgv = nullptr;
  Tcl_Size matArgc = 0;

  TCL_Char** weightArgv = nullptr;
  Tcl_Size weightArgc = 0;

  int result = TCL_ERROR;

  //
  // Material arguments.
  //
  if (matEnd - matStart == 1) {
    // A single Tcl word may be either one material tag or a Tcl list of tags.
    if (Tcl_SplitList(interp, argv[matStart], &matArgc, &matArgv) != TCL_OK) {
      opserr << "WARNING invalid list of component material tags\n";
      goto cleanup;
    }
    freeMatArgv = true;
  }
  else {
    matArgc = matEnd - matStart;
    matArgv = argv + matStart;
  }

  if (matArgc <= 0) {
    opserr << "WARNING nDMaterial Parallel3D received an empty material list\n";
    goto cleanup;
  }

  materials.reserve(static_cast<std::size_t>(matArgc));

  for (Tcl_Size i = 0; i < matArgc; i++) {
    int matTag;

    if (Tcl_GetInt(interp, matArgv[i], &matTag) != TCL_OK) {
      opserr << "WARNING invalid component material tag " << matArgv[i] << "\n";
      goto cleanup;
    }

    NDMaterial *theMat = builder->getTypedObject<NDMaterial>(matTag);
    if (theMat == nullptr) {
      opserr << "WARNING could not find nDMaterial with tag " << matTag << "\n";
      goto cleanup;
    }

    materials.push_back(theMat);
  }

  //
  // Weight arguments.
  //
  if (weightLoc < 0) {
    weights.assign(materials.size(), 1.0);
  }
  else {
    int weightStart = weightLoc + 1;
    int weightEnd   = argc;

    if (weightEnd <= weightStart) {
      opserr << "WARNING nDMaterial Parallel3D -weights requires at least one weight\n";
      goto cleanup;
    }

    if (weightEnd - weightStart == 1) {
      // A single Tcl word may be either one weight or a Tcl list of weights.
      if (Tcl_SplitList(interp, argv[weightStart], &weightArgc, &weightArgv) != TCL_OK) {
        opserr << "WARNING invalid list of Parallel3D weights\n";
        goto cleanup;
      }
      freeWeightArgv = true;
    }
    else {
      weightArgc = weightEnd - weightStart;
      weightArgv = argv + weightStart;
    }

    if (weightArgc <= 0) {
      opserr << "WARNING nDMaterial Parallel3D received an empty weight list\n";
      goto cleanup;
    }

    if (weightArgc != matArgc) {
      opserr << "WARNING nDMaterial Parallel3D got "
             << matArgc << " component materials but "
             << weightArgc << " weights\n";
      goto cleanup;
    }

    weights.reserve(static_cast<std::size_t>(weightArgc));

    for (Tcl_Size i = 0; i < weightArgc; i++) {
      double weight;

      if (Tcl_GetDouble(interp, weightArgv[i], &weight) != TCL_OK) {
        opserr << "WARNING invalid Parallel3D weight " << weightArgv[i] << "\n";
        goto cleanup;
      }

      weights.push_back(weight);
    }
  }

  {
    NDMaterial* theMat = new Parallel3DMaterial(tag, materials, weights);

    if (builder->addTaggedObject<NDMaterial>(*theMat) != TCL_OK) {
      delete theMat;
      opserr << "WARNING could not add nDMaterial Parallel3D " << tag << "\n";
      goto cleanup;
    }
  }

  result = TCL_OK;

cleanup:
  if (freeMatArgv) {
    Tcl_Free((char*)matArgv);
  }

  if (freeWeightArgv) {
    Tcl_Free((char*)weightArgv);
  }

  return result;
}



int 
TclCommand_addOrthotropicWrapper(ClientData clientData, Tcl_Interp* interp, Tcl_Size argc, TCL_Char** const argv)
{
  // nDMaterial Orthotropic $tag $theIsoMat $Ex $Ey $Ez $Gxy $Gyz $Gzx $vxy $vyz $vzx $Asigmaxx $Asigmayy $Asigmazz $Asigmaxyxy $Asigmayzyz $Asigmaxzxz.
  assert(clientData != nullptr);
  ModelRegistry *builder = static_cast<ModelRegistry*>(clientData);
  enum class Positions : int {
    Tag,
    MaterialTag,
    Ex, Ey, Ez,
    Gxy, Gyz, Gzx,
    NuXY, NuYZ, NuZX,
    Axx, Ayy, Azz, 
    Axyxy, Ayzyz, Axzxz,
    EndRequired,
    End
  };

  int tag, matTag;
  struct {
    double Ex, Ey, Ez;
    double Gxy, Gyz, Gzx;
    double NuXY, NuYZ, NuZX;
    double Axx, Ayy, Azz, Axyxy, Ayzyz, Axzxz;
  } data;

  ArgumentTracker<Positions> tracker;
  std::set<int> positions;

  if (argc < 18) {
    opserr << OpenSees::PromptValueError << "insufficient arguments\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "failed to read tag\n";
    return TCL_ERROR;
  }
  tracker.consume(Positions::Tag);


  int i=3;
  while (i < argc) {
    if (strcmp(argv[i], "-material") == 0) {
      if (i + 1 >= argc) {
        opserr << OpenSees::PromptValueError 
               << "Missing material tag."
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      if (Tcl_GetInt(interp, argv[i + 1], &matTag) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "failed to read material tag " << argv[i + 1]
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      tracker.consume(Positions::MaterialTag);
      i += 2;
    }
    else if (strcmp(argv[i], "-Ex") == 0) {
      if (i + 1 >= argc || Tcl_GetDouble(interp, argv[i + 1], &data.Ex) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to read Ex\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::Ex);
      i += 2;
    }
    else if (strcmp(argv[i], "-Ey") == 0) {
      if (i + 1 >= argc || Tcl_GetDouble(interp, argv[i + 1], &data.Ey) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to read Ey\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::Ey);
      i += 2;
    }
    else if (strcmp(argv[i], "-Ez") == 0) {
      if (i + 1 >= argc || Tcl_GetDouble(interp, argv[i + 1], &data.Ez) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to read Ez\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::Ez);
      i += 2;
    }
    else if (strcmp(argv[i], "-Gxy") == 0) {
      if (i + 1 >= argc || Tcl_GetDouble(interp, argv[i + 1], &data.Gxy) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to read Gxy\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::Gxy);
      i += 2;
    }
    else if (strcmp(argv[i], "-Gyz") == 0) {
      if (i + 1 >= argc || Tcl_GetDouble(interp, argv[i + 1], &data.Gyz) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to read Gyz\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::Gyz);
      i += 2;
    }
    else if (strcmp(argv[i], "-Gzx") == 0) {
      if (i + 1 >= argc || Tcl_GetDouble(interp, argv[i + 1], &data.Gzx) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to read Gzx\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::Gzx);
      i += 2;
    }
    else if ((strcasecmp(argv[i], "-NuXY") == 0) || (strcmp(argv[i], "-vxy") == 0)) {
      if (i + 1 >= argc || Tcl_GetDouble(interp, argv[i + 1], &data.NuXY) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to read NuXY\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::NuXY);
      i += 2;
    }
    else if ((strcasecmp(argv[i], "-NuYZ") == 0) || (strcmp(argv[i], "-vyz") == 0)) {
      if (i + 1 >= argc || Tcl_GetDouble(interp, argv[i + 1], &data.NuYZ) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to read NuYZ\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::NuYZ);
      i += 2;
    }
    else if ((strcasecmp(argv[i], "-NuZX") == 0) || (strcmp(argv[i], "-vzx") == 0)) {
      if (i + 1 >= argc || Tcl_GetDouble(interp, argv[i + 1], &data.NuZX) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to read NuZX\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::NuZX);
      i += 2;
    }
    else if (strcmp(argv[i], "-Axx") == 0) {
      if (i + 1 >= argc || Tcl_GetDouble(interp, argv[i + 1], &data.Axx) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to read Axx\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::Axx);
      i += 2;
    }
    else if (strcmp(argv[i], "-Ayy") == 0) {
      if (i + 1 >= argc || Tcl_GetDouble(interp, argv[i + 1], &data.Ayy) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to read Ayy\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::Ayy);
      i += 2;
    }
    else if (strcmp(argv[i], "-Azz") == 0) {
      if (i + 1 >= argc || Tcl_GetDouble(interp, argv[i + 1], &data.Azz) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to read Azz\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::Azz);
      i += 2;
    }
    else if (strcmp(argv[i], "-Axyxy") == 0) {
      if (i + 1 >= argc || Tcl_GetDouble(interp, argv[i + 1], &data.Axyxy) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to read Axyxy\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::Axyxy);
      i += 2;
    }
    else if (strcmp(argv[i], "-Ayzyz") == 0) {
      if (i + 1 >= argc || Tcl_GetDouble(interp, argv[i + 1], &data.Ayzyz) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to read Ayzyz\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::Ayzyz);
      i += 2;
    }
    else if (strcmp(argv[i], "-Axzxz") == 0) {
      if (i + 1 >= argc || Tcl_GetDouble(interp, argv[i + 1], &data.Axzxz) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to read Axzxz\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::Axzxz);
      i += 2;
    }
    else {
      positions.insert(i);
      i++;
    }
  }

  for (int i: positions) {
    if (tracker.current() == Positions::EndRequired) {
      tracker.increment();
    }

    switch (tracker.current()) {
      case Positions::Tag:
        if (Tcl_GetInt(interp, argv[i], &tag) != TCL_OK) {
          opserr << OpenSees::PromptValueError 
                 << "failed to read tag " << argv[i] 
                 << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        }
        tracker.increment();
        break;
      case Positions::MaterialTag:
        if (Tcl_GetInt(interp, argv[i], &matTag) != TCL_OK) {
          opserr << OpenSees::PromptValueError 
                 << "failed to read material tag " << argv[i] 
                 << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        }
        tracker.increment();
        break;
      case Positions::Ex:
        if (Tcl_GetDouble(interp, argv[i], &data.Ex) != TCL_OK) {
          opserr << OpenSees::PromptValueError 
                 << "failed to read Ex " << argv[i] 
                 << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        }
        tracker.increment();
        break;
      case Positions::Ey:
        if (Tcl_GetDouble(interp, argv[i], &data.Ey) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "failed to read Ey " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.increment();
        break;
      case Positions::Ez:
        if (Tcl_GetDouble(interp, argv[i], &data.Ez) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "failed to read Ez " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.increment();
        break;
      case Positions::Gxy:
        if (Tcl_GetDouble(interp, argv[i], &data.Gxy) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "failed to read Gxy " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.increment();
        break;
      case Positions::Gyz:
        if (Tcl_GetDouble(interp, argv[i], &data.Gyz) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "failed to read Gyz " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.increment();
        break;
      case Positions::Gzx:
        if (Tcl_GetDouble(interp, argv[i], &data.Gzx) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "failed to read Gzx " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.increment();
        break;
      case Positions::NuXY:
        if (Tcl_GetDouble(interp, argv[i], &data.NuXY) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "failed to read NuXY " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.increment();
        break;
      case Positions::NuYZ:
        if (Tcl_GetDouble(interp, argv[i], &data.NuYZ) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "failed to read NuYZ " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.increment();
        break;
      case Positions::NuZX:
        if (Tcl_GetDouble(interp, argv[i], &data.NuZX) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "failed to read NuZX\n";
          return TCL_ERROR;
        }
        tracker.increment();
        break;
      case Positions::Axx:
        if (Tcl_GetDouble(interp, argv[i], &data.Axx) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "failed to read Axx\n";
          return TCL_ERROR;
        }
        tracker.increment();
        break;
      case Positions::Ayy:
        if (Tcl_GetDouble(interp, argv[i], &data.Ayy) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "failed to read Ayy\n";
          return TCL_ERROR;
        }
        tracker.increment();
        break;
      case Positions::Azz:
        if (Tcl_GetDouble(interp, argv[i], &data.Azz) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "failed to read Azz\n";
          return TCL_ERROR;
        }
        tracker.increment();
        break;
      case Positions::Axyxy:
        if (Tcl_GetDouble(interp, argv[i], &data.Axyxy) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "failed to read Axyxy\n";
          return TCL_ERROR;
        }
        tracker.increment();
        break;
      case Positions::Ayzyz:
        if (Tcl_GetDouble(interp, argv[i], &data.Ayzyz) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "failed to read Ayzyz\n";
          return TCL_ERROR;
        }
        tracker.increment();
        break;
      case Positions::Axzxz:
        if (Tcl_GetDouble(interp, argv[i], &data.Axzxz) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "failed to read Axzxz\n";
          return TCL_ERROR;
        }
        tracker.increment();
        break;
      case Positions::EndRequired:
      case Positions::End:
        opserr << OpenSees::PromptValueError << "unexpected argument " << argv[i] << "\n";
        return TCL_ERROR;
    }
  }

  if (tracker.current() < Positions::EndRequired) {
    opserr << OpenSees::PromptValueError << "missing required arguments: ";
    while (tracker.current() != Positions::End) {
      switch (tracker.current()) {
        case Positions::MaterialTag:  opserr << "material "; break;
        case Positions::Ex:   opserr << "Ex "; break;
        case Positions::Ey:   opserr << "Ey "; break;
        case Positions::Ez:   opserr << "Ez "; break;
        case Positions::Gxy:  opserr << "Gxy "; break;
        case Positions::Gyz:  opserr << "Gyz "; break;
        case Positions::Gzx:  opserr << "Gzx "; break;
        case Positions::NuXY: opserr << "NuXY "; break;
        case Positions::NuYZ: opserr << "NuYZ "; break;
        case Positions::NuZX: opserr << "NuZX "; break;
        case Positions::Axx:  opserr << "Axx "; break;
        case Positions::Ayy:  opserr << "Ayy "; break;
        case Positions::Azz:  opserr << "Azz "; break;
        case Positions::Axyxy: opserr << "Axyxy "; break;
        case Positions::Ayzyz: opserr << "Ayzyz "; break;
        case Positions::Axzxz: opserr << "Axzxz "; break;
        default:
          // Should not happen
          assert(false);
      }
    }
    opserr << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }


  NDMaterial *threeDMaterial = builder->getTypedObject<NDMaterial>(matTag);
  if (threeDMaterial == nullptr) {
    return TCL_ERROR;
  }

  if (data.Axx <= 0 || data.Ayy <= 0 || data.Azz <= 0) {
    opserr << OpenSees::PromptValueError 
           << "Axx, Ayy, and Azz must be positive\n";
    return TCL_ERROR;
  }
  if (data.Axyxy <= 0 || data.Ayzyz <= 0 || data.Axzxz <= 0) {
    opserr << OpenSees::PromptValueError 
           << "Axyxy, Ayzyz, and Axzxz must be positive\n";
    return TCL_ERROR;
  }

  NDMaterial *theMaterial = new OrthotropicMaterial(tag, *threeDMaterial,
                                                   data.Ex, data.Ey, data.Ez,
                                                   data.Gxy, data.Gyz, data.Gzx,
                                                   data.NuXY, data.NuYZ, data.NuZX,
                                                   data.Axx, data.Ayy, data.Azz,
                                                   data.Axyxy, data.Ayzyz, data.Axzxz);

  if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK ) {
    delete theMaterial;
    return TCL_ERROR;
  }

  return TCL_OK;
}
