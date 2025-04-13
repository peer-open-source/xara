#include <tcl.h>
#include <ParallelMaterial.h>
#include <BasicModelBuilder.h>
#include <runtimeAPI.h>


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
