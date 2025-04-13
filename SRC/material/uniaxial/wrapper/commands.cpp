#include <tcl.h>
#include <G3_Logging.h>
#include <BasicModelBuilder.h>
#include "ContinuumUniaxial.h"
#include <NDMaterial.h>


int TclCommand_ContinuumUniaxialMaterial(ClientData clientData, Tcl_Interp* interp,
                                         int argc, const char**const argv)
{

    BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

    if (argc < 4) {
      opserr << G3_ERROR_PROMPT << " insufficient arguments\n";
      opserr << "Want: uniaxialMaterial Continuum tag? ndMatTag?" << endln;
      return TCL_ERROR;
    }

    int tag[2];
    if (Tcl_GetInt(interp, argv[2], &tag[0]) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "failed to read tag\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[3], &tag[1]) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "failed to read tag\n";
      return TCL_ERROR;
    }

    NDMaterial* theMat = builder->getTypedObject<NDMaterial>(tag[1]);
    if (theMat == nullptr) {
      opserr << G3_ERROR_PROMPT << " material does not exist\n";
      return TCL_ERROR;
    }

    UniaxialMaterial* mat = new ContinuumUniaxial(tag[0],*theMat);

    return builder->addTaggedObject(*theMat);

}
