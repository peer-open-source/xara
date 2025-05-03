//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
#include <tcl.h>
#include <Logging.h>
#include <Parsing.h>
#include <BasicModelBuilder.h>
#include <NDMaterial.h>
#include <UniaxialMaterial.h>

#include <InitStressMaterial.h>
#include <InitStrainMaterial.h>

#include <InitStrainNDMaterial.h>
#include <InitStressNDMaterial.h>
#include <ContinuumUniaxial.h>

#include <PlaneStressMaterial.h>
#include <PlaneStrainMaterial.h>

#include <ParallelMaterial.h>
#include <FatigueMaterial.h>

int
TclCommand_addWrappingMaterial(ClientData clientData, Tcl_Interp* interp,
                               int argc, TCL_Char** const argv)
{

    BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

    if (argc < 4) {
        opserr << G3_ERROR_PROMPT << " insufficient arguments\n";
        return TCL_ERROR;
    }

    int tago, tagi;
    if (Tcl_GetInt(interp, argv[2], &tago) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "failed to read tag\n";
        return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[3], &tagi) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "failed to read tag\n";
        return TCL_ERROR;
    }



    if (strcmp(argv[1], "ContinuumWrapper") == 0 || strcmp(argv[1], "Continuum") == 0) {
        NDMaterial* inside = builder->getTypedObject<NDMaterial>(tagi);
        if (inside == nullptr) 
            return TCL_ERROR;
        
        return builder->addTypedObject<UniaxialMaterial>(tago, new ContinuumUniaxial(tago, *inside));
    }

    else if (strcmp(argv[0], "uniaxialMaterial") == 0 && (
              strstr(argv[1], "InitStress") != 0 || 
              strstr(argv[1], "InitialStress") != 0 || 
              strstr(argv[1], "InitialStrain") != 0 || 
              strstr(argv[1], "InitStrain") != 0)) {
            
        double initial;
        if (Tcl_GetDouble(interp, argv[4], &initial) != TCL_OK) {
            opserr << G3_ERROR_PROMPT << "failed to read initial value\n";
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
            opserr << G3_ERROR_PROMPT << "InitStressNDMaterial only works with 3D materials\n";
            return TCL_ERROR;
        }

        if (argc == 5) {
          double evol;
          if (Tcl_GetDouble(interp, argv[4], &evol) != TCL_OK) {
              opserr << G3_ERROR_PROMPT << "failed to read initial value\n";
              return TCL_ERROR;
          }
          for (int i = 0; i < 3; ++i)
            initial(i) = evol;

        } else {
            for (int i=0; i<6; ++i) {
                if (Tcl_GetDouble(interp, argv[4+i], &initial(i)) != TCL_OK) {
                    opserr << G3_ERROR_PROMPT << "failed to read initial value\n";
                    return TCL_ERROR;
                }
            }
        }

        if (strstr(argv[1], "Strain") != 0)
          return builder->addTypedObject<NDMaterial>(tago, new InitStrainNDMaterial(tago, *inside, initial));
    }


    return TCL_ERROR;
}


int
TclCommand_newParallelMaterial(ClientData clientData, Tcl_Interp* interp, int argc, G3_Char ** const argv)
{
    if (argc < 4) {
        opserr << "WARNING insufficient arguments\n";
        opserr << "Want: uniaxialMaterial Parallel tag? tag1? tag2? ...";
        opserr << " <-min min?> <-max max?>" << endln;
        return TCL_ERROR;
    }

    int tag;
    UniaxialMaterial* theMaterial = nullptr;
    BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
        opserr << "WARNING invalid uniaxialMaterial Parallel tag" << endln;
        return TCL_ERROR;
    }

    int numMaterials = argc-3;
    
    if (numMaterials == 0) {
        opserr << "WARNING no component material(s) provided\n";
        opserr << "uniaxialMaterial Parallel: " << tag << endln;
        return TCL_ERROR;
    }

    // Create an array to hold pointers to component materials
    UniaxialMaterial **theMats = new UniaxialMaterial *[numMaterials];
    
    // For each material get the tag and ensure it exists in model already
    for (int i=0; i<numMaterials; i++) {
      int tagI;
      if (Tcl_GetInt(interp, argv[i+3], &tagI) != TCL_OK) {
        opserr << "WARNING invalid component tag\n";
        opserr << "uniaxialMaterial Parallel: " << tag << endln;
        return TCL_ERROR;
      }

      UniaxialMaterial *theMat = builder->getTypedObject<UniaxialMaterial>(tagI);
      
      if (theMat == nullptr) {
        delete [] theMats;
        return TCL_ERROR;
      } else
        theMats[i] = theMat;
    }
    
    // Parsing was successful, allocate the material
    theMaterial = new ParallelMaterial(tag, numMaterials, theMats);
    builder->addTaggedObject<UniaxialMaterial>(*theMaterial);
    
    // Deallocate the temporary pointers
    delete [] theMats;
    return TCL_OK;
}


int
TclCommand_newFatigueMaterial(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

  if (argc < 4) {
    opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
    opserr << "Want: uniaxialMaterial Fatigue tag? matTag?";
    opserr << " <-D_max dmax?> <-e0 e0?> <-m m?>" << "\n";
    opserr << " <-min min?> <-max max?>" << "\n";
    return TCL_ERROR;
  }

  int tag, matTag;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "invalid uniaxialMaterial Fatigue tag" << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "invalid component tag\n";
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
        opserr << G3_ERROR_PROMPT << "invalid -Dmax";
        return TCL_ERROR;
      }
    } else if (strcmp(argv[j], "-E0") == 0) {
      if ((j + 1 >= argc) || (Tcl_GetDouble(interp, argv[++j], &E0) != TCL_OK)) {
        opserr << G3_ERROR_PROMPT << "invalid -E0";
        return TCL_ERROR;
      }
    } else if (strcmp(argv[j], "-m") == 0) {
      if ((j + 1 >= argc) ||
          (Tcl_GetDouble(interp, argv[++j], &m) != TCL_OK)) {
        opserr << G3_ERROR_PROMPT << "invalid -m";
        return TCL_ERROR;
      }
    } else if (strcmp(argv[j], "-min") == 0) {
      if ((j + 1 >= argc) ||
          (Tcl_GetDouble(interp, argv[++j], &epsmin) != TCL_OK)) {
        opserr << G3_ERROR_PROMPT << "invalid -min ";
        return TCL_ERROR;
      }
    } else if (strcmp(argv[j], "-max") == 0) {
      if ((j + 1 >= argc) ||
          (Tcl_GetDouble(interp, argv[++j], &epsmax) != TCL_OK)) {
        opserr << G3_ERROR_PROMPT << "invalid -max";
        return TCL_ERROR;
      }
    }
  }

  UniaxialMaterial *theMat = builder->getTypedObject<UniaxialMaterial>(matTag);

  if (theMat == nullptr) {
    opserr << G3_ERROR_PROMPT << "component material does not exist\n";
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
