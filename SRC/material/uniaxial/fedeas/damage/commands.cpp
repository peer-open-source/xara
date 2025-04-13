//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Written: cmp
//
#include <stdio.h>
#include <string.h>
#include <Parsing.h>
#include <Logging.h>
#include "DegradingUniaxialWrapper.h"
#include <BasicModelBuilder.h>

#define WRAPPER_CMD "FedeasUniaxialDamage"
// #define WRAPPER_CMD "FedeasDamage"

int
TclCommand_newFedeasUniaxialDamage(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);
  // Pointer to a uniaxial material that will be returned
  DegradingUniaxialWrapper *theMaterial = nullptr;
  UniaxialMaterial *theWrappedMaterial = nullptr;
  int tags[2];

  if (argc < 2) {
    opserr << "WARNING invalid uniaxialMaterial " WRAPPER_CMD " $tag "
              "$wrapTag <-damage $damageTag>"
           << endln;
    return TCL_ERROR;;
  }

  // Get wrapper tag
  if (Tcl_GetInt(interp, argv[2], &tags[0]) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial tag\n";
    // printCommand(argc, argv);
    return TCL_ERROR;;
  }
  // Get base tag
  if (Tcl_GetInt(interp, argv[3], &tags[1]) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial tag\n";
    // printCommand(argc, argv);
    return TCL_ERROR;;
  }

  // Get base material
  theWrappedMaterial = builder->getTypedObject<UniaxialMaterial>(tags[1]);
  if (theWrappedMaterial == nullptr) {
    opserr << "WARNING unable to retrieve uniaxialMaterial with tag" WRAPPER_CMD " tag: "
           << tags[1] << endln;
    return TCL_ERROR;;
  }

  int argn = 4;
  double Ccd = 0.5;
  StateOperator *damage = new StateOperator;
  while (argn < argc) {
    const char *param = argv[argn];

    if ((strcmp(param, "-damage") == 0) || 
        (strcmp(param, "-dmg") == 0)    ||
        (strcmp(param, "-DMG") == 0))   {
      *damage = *(StateOperator*)Tcl_GetAssocData(interp, 
                                                  "fedeas::damage::UniaxialDamage", NULL);
      
      damage->call(damage, interp, ISW_CREATE, argc - argn, &argv[1+argn], 0, 0, 0, 0, 0);
      damage->call(damage, interp, ISW_MALLOC, 0, 0, 0, 0, 0, 0, 0);
      argn++;

    } else if ((strcmp(param, "-couple") == 0) || 
               (strcmp(param, "-ccd") == 0)    ||
               (strcmp(param, "-Ccd") == 0))   {
      Ccd = std::stod(argv[++argn]);
    } else {
      break;
    }
    argn++;
  }

  // Parsing was successful, allocate the material
  theMaterial = new DegradingUniaxialWrapper(tags[0], *theWrappedMaterial, damage);

  theMaterial->setCoupling(Ccd);

  // if (dmgtag){
  //   if (theMaterial->setDamageWrapper(G3_getInterpreter(rt), dmgtag) > 0)
  //     opserr << "#Set damage wrapper '" << dmgtag << "'\n";
  // }

  // return G3_addUniaxialMaterial(rt, theMaterial);
  builder->addTaggedObject<UniaxialMaterial>(*theMaterial);
  return TCL_OK;
}

