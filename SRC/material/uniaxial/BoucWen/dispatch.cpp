#include <tcl.h>
#include <Parsing.h>
#include <Logging.h>
#include <BasicModelBuilder.h>
#include "BoucWenMaterial.h"

int
TclCommand_newUniaxialBoucWen(ClientData cd, Tcl_Interp* interp, int argc, G3_Char ** const argv)
{
    if (argc < 12) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: uniaxialMaterial BoucWen tag? alpha? ko? n? gamma?"
             << endln << " beta? Ao? deltaA? deltaNu? deltaEta?" << endln;
      return TCL_ERROR;
    }

    int tag;
    double alpha, ko, n, gamma, beta, Ao, deltaA, deltaNu, deltaEta;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid uniaxialMaterial BoucWen tag" << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[3], &alpha) != TCL_OK) {
      opserr << "WARNING invalid alpha\n";
      opserr << "uniaxialMaterial BoucWen: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[4], &ko) != TCL_OK) {
      opserr << "WARNING invalid ko\n";
      opserr << "uniaxialMaterial BoucWen: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[5], &n) != TCL_OK) {
      opserr << "WARNING invalid n\n";
      opserr << "uniaxialMaterial BoucWen: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[6], &gamma) != TCL_OK) {
      opserr << "WARNING invalid gamma\n";
      opserr << "uniaxialMaterial BoucWen: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[7], &beta) != TCL_OK) {
      opserr << "WARNING invalid beta\n";
      opserr << "uniaxialMaterial BoucWen: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[8], &Ao) != TCL_OK) {
      opserr << "WARNING invalid Ao\n";
      opserr << "uniaxialMaterial BoucWen: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[9], &deltaA) != TCL_OK) {
      opserr << "WARNING invalid deltaA\n";
      opserr << "uniaxialMaterial BoucWen: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[10], &deltaNu) != TCL_OK) {
      opserr << "WARNING invalid deltaNu\n";
      opserr << "uniaxialMaterial BoucWen: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[11], &deltaEta) != TCL_OK) {
      opserr << "WARNING invalid deltaEta\n";
      opserr << "uniaxialMaterial BoucWen: " << tag << endln;
      return TCL_ERROR;
    }

    // Check if the user has given a tolerance for the Newton scheme
    double tolerance = 1.0e-8;
    if (argc > 12) {
      if (Tcl_GetDouble(interp, argv[12], &tolerance) != TCL_OK) {
        opserr << "WARNING invalid tolerance\n";
        opserr << "uniaxialMaterial BoucWen: " << tolerance << endln;
        return TCL_ERROR;
      }
    }

    // Check if the user has given a maxNumIter for the Newton scheme
    int maxNumIter = 20;
    if (argc > 13) {
      if (Tcl_GetInt(interp, argv[13], &maxNumIter) != TCL_OK) {
        opserr << "WARNING invalid maxNumIter\n";
        opserr << "uniaxialMaterial BoucWen: " << maxNumIter << endln;
        return TCL_ERROR;
      }
    }

    // Parsing was successful, allocate the material
    UniaxialMaterial* theMaterial =
        new BoucWenMaterial(tag, alpha, ko, n, gamma, beta, Ao, deltaA, deltaNu,
                            deltaEta, tolerance, maxNumIter);

    return ((BasicModelBuilder*)cd)->addTaggedObject<UniaxialMaterial>(*theMaterial);
}
