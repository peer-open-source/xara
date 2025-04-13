
#include <tcl.h>
#include "AxialSp.h"
#include "AxialSpHD.h"
#include <runtimeAPI.h>
#include <BasicModelBuilder.h>

int
TclCommand_AxialSp(ClientData clientData, Tcl_Interp *interp, int argc,
                   TCL_Char ** const argv)
{
  // arguments (necessary)
  int tag;
  double sce;
  double fty;
  double fcy;

  // arguments (optional)
  double bte = 0.0;
  double bty = 0.0;
  double bcy = 0.0;
  double fcr = 0.0;

  //
  UniaxialMaterial *theMaterial = nullptr;

  // error flag
  bool ifNoError = true;

  if (argc < 6 || argc > 10) {
    // uniaxialMaterial AxialSp matTag? sce? fty? fcy?

    opserr << "WARNING invalid number of arguments\n";
    ifNoError = false;
  }

  // argv[2-5]
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid AxialSp tag" << "\n";
    ifNoError = false;
  }

  if (Tcl_GetDouble(interp, argv[3], &sce) != TCL_OK) {
    opserr << "WARNING invalid sce\n";
    opserr << "AxialSp: " << tag << "\n";
    ifNoError = false;
  }

  if (Tcl_GetDouble(interp, argv[4], &fty) != TCL_OK) {
    opserr << "WARNING invalid fty\n";
    opserr << "AxialSp: " << tag << "\n";
    ifNoError = false;
  }

  if (Tcl_GetDouble(interp, argv[5], &fcy) != TCL_OK) {
    opserr << "WARNING invalid fcy\n";
    opserr << "AxialSp: " << tag << "\n";
    ifNoError = false;
  }

  // argv[6~]
  if (argc >= 7) {
    if (Tcl_GetDouble(interp, argv[6], &bte) != TCL_OK) {
      opserr << "WARNING invalid bte\n";
      opserr << "AxialSp: " << tag << "\n";
      ifNoError = false;
    }
  }

  if (argc >= 8) {
    if (Tcl_GetDouble(interp, argv[7], &bty) != TCL_OK) {
      opserr << "WARNING invalid bty\n";
      opserr << "AxialSp: " << tag << "\n";
      ifNoError = false;
    }
  }

  if (argc >= 9) {
    if (Tcl_GetDouble(interp, argv[8], &bcy) != TCL_OK) {
      opserr << "WARNING invalid bcy\n";
      opserr << "AxialSp: " << tag << "\n";
      ifNoError = false;
    }
  }

  if (argc == 10) {
    if (Tcl_GetDouble(interp, argv[9], &fcr) != TCL_OK) {
      opserr << "WARNING invalid fcr\n";
      opserr << "AxialSp: " << tag << "\n";
      ifNoError = false;
    }
  }

  // if error detected
  if (!ifNoError) {
    // input:
    opserr << "Input command: ";
    for (int i = 0; i < argc; i++) {
      opserr << argv[i] << " ";
    }
    opserr << "\n";

    // want:
    opserr << "WANT: AxialSp tag? sce? fty? fcy? <bte?> <bty?> <bcy?> <fcr?>"
           << "\n";

    return TCL_ERROR;
  }

  // Parsing was successful, allocate the material
  theMaterial = new AxialSp(tag, sce, fty, fcy, bte, bty, bcy, fcr);

  // Now add the material to the modelBuilder

  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);
  if (builder->addTaggedObject<UniaxialMaterial>(*theMaterial) != TCL_OK) {
    delete theMaterial;
    return TCL_ERROR;

  } else {
    return TCL_OK;
  }
}

int
TclCommand_AxialSpHD(ClientData clientData, Tcl_Interp *interp, int argc,
                     TCL_Char ** const argv)
{
  // arguments (necessary)
  int tag;
  double sce;
  double fty;
  double fcy;

  // arguments (optional)
  double bte = 1.0;
  double bty = 1.0;
  double bth = 1.0;
  double bcy = 1.0;
  double fcr = 0.0;
  double ath = 1.0;

  //
  UniaxialMaterial *theMaterial = 0;

  // error flag
  bool ifNoError = true;

  if (argc < 6 ||
      argc > 12) { // uniaxialMaterial AxialSpHD matTag? sce? fty? fcy?

    opserr << "WARNING invalid number of arguments\n";
    ifNoError = false;
  }

  // argv[2~5]
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid AxialSpHD tag" << "\n";
    ifNoError = false;
  }

  if (Tcl_GetDouble(interp, argv[3], &sce) != TCL_OK) {
    opserr << "WARNING invalid sce\n";
    opserr << "AxialSpHD: " << tag << "\n";
    ifNoError = false;
  }

  if (Tcl_GetDouble(interp, argv[4], &fty) != TCL_OK) {
    opserr << "WARNING invalid fty\n";
    opserr << "AxialSpHD: " << tag << "\n";
    ifNoError = false;
  }

  if (Tcl_GetDouble(interp, argv[5], &fcy) != TCL_OK) {
    opserr << "WARNING invalid fcy\n";
    opserr << "AxialSpHD: " << tag << "\n";
    ifNoError = false;
  }

  // argv[6~]
  if (argc >= 7) {
    if (Tcl_GetDouble(interp, argv[6], &bte) != TCL_OK) {
      opserr << "WARNING invalid bte\n";
      opserr << "AxialSpHD: " << tag << "\n";
      ifNoError = false;
    }
  }

  if (argc >= 8) {
    if (Tcl_GetDouble(interp, argv[7], &bty) != TCL_OK) {
      opserr << "WARNING invalid bty\n";
      opserr << "AxialSpHD: " << tag << "\n";
      ifNoError = false;
    }
  }

  if (argc >= 9) {
    if (Tcl_GetDouble(interp, argv[8], &bth) != TCL_OK) {
      opserr << "WARNING invalid bth\n";
      opserr << "AxialSpHD: " << tag << "\n";
      ifNoError = false;
    }
  }

  if (argc >= 10) {
    if (Tcl_GetDouble(interp, argv[9], &bcy) != TCL_OK) {
      opserr << "WARNING invalid bcy\n";
      opserr << "AxialSpHD: " << tag << "\n";
      ifNoError = false;
    }
  }

  if (argc >= 11) {
    if (Tcl_GetDouble(interp, argv[10], &fcr) != TCL_OK) {
      opserr << "WARNING invalid fcr\n";
      opserr << "AxialSpHD: " << tag << "\n";
      ifNoError = false;
    }
  }

  if (argc == 12) {
    if (Tcl_GetDouble(interp, argv[11], &ath) != TCL_OK) {
      opserr << "WARNING invalid ath\n";
      opserr << "AxialSpHD: " << tag << "\n";
      ifNoError = false;
    }
  }

  // if error detected
  if (!ifNoError) {

    // want:
    opserr << "WANT: AxialSpHD tag? sce? fty? fcy? <bte?> <bty?> <bth?> <bcy?> "
              "<fcr?> <ath?>"
           << "\n";

    return TCL_ERROR;
  }

  // Parsing was successful, allocate the material
  theMaterial = new AxialSpHD(tag, sce, fty, fcy, bte, bty, bth, bcy, fcr, ath);

  // Now add the material to the modelBuilder
  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);
  if (builder->addTaggedObject<UniaxialMaterial>(*theMaterial) != TCL_OK) {
    delete theMaterial;
    return TCL_ERROR;
  }

  return TCL_OK;
}
