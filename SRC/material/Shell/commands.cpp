/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// Written: rms, MHS, cmp
// Created: 07/99
//
// Description: This file contains the function invoked when the user invokes
// the section command in the interpreter.
//
// What: "@(#) TclModelBuilderMaterialCommands.C, revA"
//
#include <tcl.h>
#include <string.h>

#include <BasicModelBuilder.h>
#include <ElasticMembranePlateSection.h>
#include <MembranePlateFiberSection.h>

typedef SectionForceDeformation ShellSection;


int
TclCommand_ShellSection(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** argv)
{
  // Pointer to a section that will be added to the model builder
  SectionForceDeformation* theSection = nullptr;
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

  // Check argv[1] for section type

  if ((strcmp(argv[1], "ElasticShell") == 0) ||
      (strcmp(argv[1], "ElasticMembranePlateSection") == 0) ||
      (strcmp(argv[1], "ElasticPlateSection") == 0)) {

    if (argc < 5) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: section ElasticMembranePlateSection tag? E? nu? h? <rho?>\n";
      return TCL_ERROR;
    }

    double E, nu, h;
    double rho = 0.0;

    int tag;
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid section tag for ElasticShell section.\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
      opserr << "WARNING invalid E\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &nu) != TCL_OK) {
      opserr << "WARNING invalid nu\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &h) != TCL_OK) {
      opserr << "WARNING invalid h\n";
      return TCL_ERROR;
    }

    if (argc > 6 && Tcl_GetDouble(interp, argv[6], &rho) != TCL_OK) {
      opserr << "WARNING invalid rho\n";
      return TCL_ERROR;
    }

    theSection = new ElasticMembranePlateSection(tag, E, nu, h, rho);
  }

  else if (strcmp(argv[1], "PlateFiber") == 0) {
    if (argc < 5) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: section PlateFiber tag? matTag? h? \n";
      return TCL_ERROR;
    }

    int tag, matTag;
    double h;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid section PlateFiber tag\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &matTag) != TCL_OK) {
      opserr << "WARNING invalid matTag\n";
      opserr << "PlateFiber section: " << matTag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &h) != TCL_OK) {
      opserr << "WARNING invalid h\n";
      opserr << "PlateFiber section: " << tag << "\n";
      return TCL_ERROR;
    }

    NDMaterial* theMaterial = builder->getTypedObject<NDMaterial>(matTag);
    if (theMaterial == nullptr)
      return TCL_ERROR;

    theSection = new MembranePlateFiberSection(tag, h, *theMaterial);
  }


  // Now add the material to the modelBuilder
  if (builder->addTaggedObject<ShellSection>(*theSection) < 0) {
    opserr << "WARNING could not add section to the domain\n";
    opserr << *theSection << "\n";
    delete theSection;
    return TCL_ERROR;
  }

  return TCL_OK;

//return TCL_ERROR;
}

