// Written: MHS
// Created: August 2000
//
// Description: This file contains the parsing routines for the
// TCL unloadingRule command.

#include <OPS_Globals.h>
#include <UnloadingRule.h>
// #include <TakedaUnloadingRule.h>
// #include <EnergyUnloadingRule.h>
// #include <ConstantUnloadingRule.h>

#include <elementAPI.h>
#include <tcl.h>
#include <string.h>

extern OPS_Routine OPS_TakedaUnloadingRule;
extern OPS_Routine OPS_EnergyUnloadingRule;
extern OPS_Routine OPS_ConstantUnloadingRule;
extern OPS_Routine OPS_KarsanUnloadingRule;

extern "C" int OPS_ResetInputNoBuilder(ClientData clientData,
                                       Tcl_Interp *interp, int cArg, int mArg,
                                       TCL_Char ** const argv, Domain *domain);

#include <packages.h>

int
TclBasicBuilderUnloadingRuleCommand(ClientData clientData, Tcl_Interp *interp,
                                    int argc, TCL_Char ** const argv,
                                    Domain *theDomain)
{
  G3_Runtime *rt = G3_getRuntime(interp);

  // Make sure there is a minimum number of arguments
  if (argc < 2) {
    opserr << "WARNING insufficient number of unloadingRule arguments\n";
    opserr << "Want: unloadingRule type? tag? <specific unloadingRule args>"
           << endln;
    return TCL_ERROR;
  }

    OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, theDomain);

  // Pointer to a unloadingRule that will be added to the model builder
  UnloadingRule *theState = 0;

  // Check argv[1] for unloadingRule type
  if (strcmp(argv[1], "Ductility") == 0 || strcmp(argv[1], "Takeda") == 0) {
    void *theDegr = OPS_TakedaUnloadingRule(rt, argc, argv);
    if (theDegr != 0)
      theState = (UnloadingRule *)theDegr;
    else
      return TCL_ERROR;
  }

  else if (strcmp(argv[1], "Energy") == 0) {
    void *theDegr = OPS_EnergyUnloadingRule(rt, argc, argv);
    if (theDegr != 0)
      theState = (UnloadingRule *)theDegr;
    else
      return TCL_ERROR;
  }

  else if (strcmp(argv[1], "Constant") == 0) {
    void *theDegr = OPS_ConstantUnloadingRule(rt, argc, argv);
    if (theDegr != 0)
      theState = (UnloadingRule *)theDegr;
    else
      return TCL_ERROR;
  }

  else if (strcmp(argv[1], "Karsan") == 0) {
    void *theDegr = OPS_KarsanUnloadingRule(rt, argc, argv);
    if (theDegr != 0)
      theState = (UnloadingRule *)theDegr;
    else
      return TCL_ERROR;
  }

  else {
    opserr << "WARNING unknown type of unloadingRule: " << argv[1];
    opserr << "\nValid types: Ductility, Energy, Constant\n";
    return TCL_ERROR;
  }

  // Ensure we have created the Degradation, out of memory if got here and no
  // unloadingRule
  if (theState == 0) {
    opserr << "WARNING ran out of memory creating unloadingRule\n";
    opserr << argv[1] << endln;
    return TCL_ERROR;
  }

  // Now add the material to the modelBuilder
  if (OPS_addUnloadingRule(theState) == false) {
    opserr << "WARNING could not add unloadingRule to the domain\n";
    opserr << *theState << endln;
    delete theState; // Avoid memory leak
    return TCL_ERROR;
  }

  return TCL_OK;
}
