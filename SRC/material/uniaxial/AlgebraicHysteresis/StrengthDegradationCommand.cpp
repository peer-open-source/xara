// Written: MHS
// Created: August 2000
//
// Description: This file contains the parsing routines for the
// TCL strengthDegradation command.

#include <string.h>
#include <tcl.h>
#include <OPS_Globals.h>

#include <BasicModelBuilder.h>

#include <SectionStrengthDegradation.h>
#include <EnergyStrengthDegradation.h>
#include <ConstantStrengthDegradation.h>
#include <DuctilityStrengthDegradation.h>
#include <SectionForceDeformation.h>

#include <elementAPI.h>
#include <packages.h>

extern "C" int OPS_ResetInputNoBuilder(ClientData clientData,
                                       Tcl_Interp *interp, int cArg, int mArg,
                                       TCL_Char ** const argv, Domain *domain);


extern OPS_Routine OPS_SectionStrengthDegradation;
extern OPS_Routine OPS_DuctilityStrengthDegradation;
extern OPS_Routine OPS_EnergyStrengthDegradation;
extern OPS_Routine OPS_ConstantStrengthDegradation;
extern OPS_Routine OPS_ACIStrengthDegradation;
extern OPS_Routine OPS_PetrangeliStrengthDegradation;

static void
printCommand(int argc, TCL_Char ** const argv)
{
  opserr << "Input command: ";
  for (int i = 0; i < argc; i++)
    opserr << argv[i] << " ";
  opserr << endln;
}

int
TclBasicBuilderStrengthDegradationCommand(ClientData clientData,
                                          Tcl_Interp *interp, int argc,
                                          TCL_Char ** const argv, Domain *theDomain)
{
  G3_Runtime *rt = G3_getRuntime(interp);

  // Make sure there is a minimum number of arguments
  if (argc < 2) {
    opserr << "WARNING insufficient number of strengthDegradation arguments\n";
    opserr << "Want: strengthDegradation type? tag? <specific strengthDegradation args>"
           << endln;
    return TCL_ERROR;
  }

    OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, theDomain);

  // Pointer to a strengthDegradation that will be added to the model builder
  StrengthDegradation *theState = 0;

  // Check argv[1] for strengthDegradation type
  if (strcmp(argv[1], "Section") == 0) {
    if (argc < 7) {
      opserr << "WARNING insufficient arguments\n";
      printCommand(argc, argv);
      opserr << "Want: strengthDegradation Section tag? code e1? V2? e2? "
                "<-yield ey?>"
             << endln;
      return TCL_ERROR;
    }

    int tag, code;
    double e1, V2, e2, ey;
    bool isDuctility = false;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid strengthDegradation Section tag" << endln;
      return TCL_ERROR;
    }

    if (strcmp(argv[3], "Mz") == 0)
      code = SECTION_RESPONSE_MZ;
    else if (strcmp(argv[3], "P") == 0)
      code = SECTION_RESPONSE_P;
    else if (strcmp(argv[3], "Vy") == 0)
      code = SECTION_RESPONSE_VY;
    else if (strcmp(argv[3], "My") == 0)
      code = SECTION_RESPONSE_MY;
    else if (strcmp(argv[3], "Vz") == 0)
      code = SECTION_RESPONSE_VZ;
    else if (strcmp(argv[3], "T") == 0)
      code = SECTION_RESPONSE_T;
    else {
      opserr << "WARNING invalid code" << argv[3] << endln;
      opserr << "strengthDegradation Section: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &e1) != TCL_OK) {
      opserr << "WARNING invalid e1\n";
      opserr << "strengthDegradation Section: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &V2) != TCL_OK) {
      opserr << "WARNING invalid V2\n";
      opserr << "strengthDegradation Section: " << tag << endln;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &e2) != TCL_OK) {
      opserr << "WARNING invalid e2\n";
      opserr << "strengthDegradation Section: " << tag << endln;
      return TCL_ERROR;
    }

    if (argc > 8 && strcmp(argv[7], "-yield") == 0) {
      if (Tcl_GetDouble(interp, argv[8], &ey) != TCL_OK) {
        opserr << "WARNING invalid ey\n";
        opserr << "strengthDegradation Section: " << tag << endln;
        return TCL_ERROR;
      }
      isDuctility = true;
    }

    // Parsing was successful, allocate the material
    if (isDuctility)
      theState = new SectionStrengthDegradation(tag, ey, e1, V2, e2, code);
    else
      theState = new SectionStrengthDegradation(tag, e1, V2, e2, code);
  }

  else if (strcmp(argv[1], "Ductility") == 0) {
    void *theDegr = OPS_DuctilityStrengthDegradation(rt, argc, argv);
    if (theDegr != 0)
      theState = (StrengthDegradation *)theDegr;
    else
      return TCL_ERROR;
  }

  else if (strcmp(argv[1], "Energy") == 0) {
    void *theDegr = OPS_EnergyStrengthDegradation(rt, argc, argv);
    if (theDegr != 0)
      theState = (StrengthDegradation *)theDegr;
    else
      return TCL_ERROR;
  }

  else if (strcmp(argv[1], "Constant") == 0) {
    void *theDegr = OPS_ConstantStrengthDegradation(rt, argc, argv);
    if (theDegr != 0)
      theState = (StrengthDegradation *)theDegr;
    else
      return TCL_ERROR;
  }

  else if (strcmp(argv[1], "ACI") == 0) {
    void *theDegr = OPS_ACIStrengthDegradation(rt, argc, argv);
    if (theDegr != 0)
      theState = (StrengthDegradation *)theDegr;
    else
      return TCL_ERROR;
  }

  else if (strcmp(argv[1], "Petrangeli") == 0) {
    void *theDegr = OPS_PetrangeliStrengthDegradation(rt, argc, argv);
    if (theDegr != 0)
      theState = (StrengthDegradation *)theDegr;
    else
      return TCL_ERROR;
  }

  else {
    opserr << "WARNING unknown type of strengthDegradation: " << argv[1];
    opserr << "\nValid types: Section, Energy, Constant, Ductility, ACI, "
              "Petrangeli\n";
    return TCL_ERROR;
  }

  // Ensure we have created the Degradation, out of memory if got here and no
  // strengthDegradation
  if (theState == 0) {
    opserr << "WARNING ran out of memory creating strengthDegradation\n";
    opserr << argv[1] << endln;
    return TCL_ERROR;
  }

  // Now add the material to the modelBuilder
  if (OPS_addStrengthDegradation(theState) == false) {
    opserr << "WARNING could not add strengthDegradation to the domain\n";
    opserr << *theState << endln;
    delete theState; // Avoid memory leak
    return TCL_ERROR;
  }

  return TCL_OK;
}
