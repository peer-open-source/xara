// Written: MHS
// Created: August 2000
//
// Description: This file contains the parsing routines for the
// TCL stiffnessDegradation command.

#include <tcl.h>
#include <OPS_Globals.h>
#include <BasicModelBuilder.h>

extern "C" int OPS_ResetInputNoBuilder(ClientData clientData,
                                       Tcl_Interp *interp, int cArg, int mArg,
                                       TCL_Char ** const argv, Domain *domain);

#include <DuctilityStiffnessDegradation.h>
#include <EnergyStiffnessDegradation.h>
#include <ConstantStiffnessDegradation.h>

#include <string.h>
#include <elementAPI.h>
#include <packages.h>

extern OPS_Routine OPS_DuctilityStiffnessDegradation;
extern OPS_Routine OPS_EnergyStiffnessDegradation;
extern OPS_Routine OPS_ConstantStiffnessDegradation;
extern OPS_Routine OPS_PincheiraStiffnessDegradation;

int
TclBasicBuilderStiffnessDegradationCommand(ClientData clientData,
                                           Tcl_Interp *interp, int argc,
                                           TCL_Char ** const argv, Domain *theDomain)
{
  G3_Runtime *rt = G3_getRuntime(interp);
  // Make sure there is a minimum number of arguments
  if (argc < 2) {
    opserr << "WARNING insufficient number of stiffnessDegradation arguments\n";
    opserr << "Want: stiffnessDegradation type? tag? <specific "
              "stiffnessDegradation args>"
           << endln;
    return TCL_ERROR;
  }

  OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, theDomain);

  // Pointer to a stiffnessDegradation that will be added to the model builder
  StiffnessDegradation *theState = 0;

  // Check argv[1] for stiffnessDegradation type
  if (strcmp(argv[1], "Ductility") == 0) {
    void *theDegr = OPS_DuctilityStiffnessDegradation(rt, argc, argv);
    if (theDegr != 0)
      theState = (StiffnessDegradation *)theDegr;
    else
      return TCL_ERROR;
  }

  else if (strcmp(argv[1], "Energy") == 0) {
    void *theDegr = OPS_EnergyStiffnessDegradation(rt, argc, argv);
    if (theDegr != 0)
      theState = (StiffnessDegradation *)theDegr;
    else
      return TCL_ERROR;
  }

  else if (strcmp(argv[1], "Constant") == 0) {
    void *theDegr = OPS_ConstantStiffnessDegradation(rt, argc, argv);
    if (theDegr != 0)
      theState = (StiffnessDegradation *)theDegr;
    else
      return TCL_ERROR;
  }

  else if (strcmp(argv[1], "Pincheira") == 0) {
    void *theDegr = OPS_PincheiraStiffnessDegradation(rt, argc, argv);
    if (theDegr != 0)
      theState = (StiffnessDegradation *)theDegr;
    else
      return TCL_ERROR;
  }

  else {
    opserr << "WARNING unknown type of stiffnessDegradation: " << argv[1];
    opserr << "\nValid types: Ductility, Energy, Constant\n";
    return TCL_ERROR;
  }

  // Ensure we have created the Degradation, out of memory if got here and no
  // stiffnessDegradation
  if (theState == 0) {
    opserr << "WARNING ran out of memory creating stiffnessDegradation\n";
    opserr << argv[1] << endln;
    return TCL_ERROR;
  }

  // Now add the material to the modelBuilder
  if (OPS_addStiffnessDegradation(theState) == false) {
    opserr << "WARNING could not add stiffnessDegradation to the domain\n";
    opserr << *theState << endln;
    delete theState; // Avoid memory leak
    return TCL_ERROR;
  }

  return TCL_OK;
}
