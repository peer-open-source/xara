//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
//
// Written: cmp
//

#include <Parsing.h>
#include <Logging.h>
#include <string.h>
#include <assert.h>

#include <tcl.h>
#include <Vector.h>
#include <VectorND.h>
#include <DummyStream.h>
#include <Logging.h>
#include <Response.h>
#include <NDMaterial.h>
#include <ModelRegistry.h>

static Tcl_CmdProc MaterialTest_setStrainSection;
static Tcl_CmdProc MaterialTest_getStressSection;
static Tcl_CmdProc MaterialTest_getTangent;
static Tcl_CmdProc MaterialTest_getResponse;
static Tcl_CmdProc MaterialTest_Commit;


using namespace OpenSees;

const struct {const char*name; Tcl_CmdProc*func;} MaterialTestCommands[] = {
  {"material::update",     MaterialTest_setStrainSection },
  {"material::stress",     MaterialTest_getStressSection },
  {"material::tangent",    MaterialTest_getTangent       },
  {"material::response",   MaterialTest_getResponse      },
  {"material::commit",     MaterialTest_Commit           },
};

// invoke Material $tag $commands
int
TclCommand_useMaterial(ClientData clientData, Tcl_Interp *interp, 
                       Tcl_Size argc, TCL_Char ** const argv)
{

  assert(clientData != nullptr);
  // TODO: Parse tag properly
  NDMaterial *theMaterial = 
    ((ModelRegistry*)clientData)->getTypedObject<NDMaterial>(std::atoi(argv[2]));

  if (theMaterial == nullptr) {
    opserr << OpenSees::PromptValueError << "no material found with tag '" << argv[2] << "'\n";
    return TCL_ERROR;
  }

  //
  //
  //
  const int ncmd = sizeof(MaterialTestCommands) / sizeof(MaterialTestCommands[0]);
  for (int i = 0; i < ncmd; ++i) {
    Tcl_CreateCommand(interp, MaterialTestCommands[i].name,
                      MaterialTestCommands[i].func,
                      (ClientData)theMaterial, nullptr);
  }

  //
  //
  int status = Tcl_Eval(interp, argv[3]);
  //
  //

  for (int i = 0; i < ncmd; ++i) {
    Tcl_DeleteCommand(interp, MaterialTestCommands[i].name);
  }
//   Tcl_DeleteCommand(interp, "material::update");
//   Tcl_DeleteCommand(interp, "material::stress");
//   Tcl_DeleteCommand(interp, "material::tangent");
//   Tcl_DeleteCommand(interp, "material::commit");
//   Tcl_DeleteCommand(interp, "material::response");

  return status;
}


static int
MaterialTest_setStrainSection(ClientData clientData, Tcl_Interp *interp,
                              Tcl_Size argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  NDMaterial *theMaterial = (NDMaterial*)clientData;

  // check number of arguments in command line
  if (argc < 2) {
    opserr << OpenSees::PromptValueError 
           << "insufficient arguments"
           << "\n";
    return TCL_ERROR;
  }

  // get the sectionID form command line
  // Need to set the data based on argc, otherwise it crashes when setting
  // "data(i-1) = strain"
  // VectorND<12> e{};
  int order = theMaterial->getOrder();
  if (order != argc - 1) {
    opserr << OpenSees::PromptValueError 
           << "expected " << order << " strain values, but got " << argc - 1 << "\n";
    return TCL_ERROR;
  }
  Vector data(order);
  for (int i = 1; i < argc && i < order; ++i) {
    double strain;
    if (Tcl_GetDouble(interp, argv[i], &strain) != TCL_OK) {
      opserr << OpenSees::PromptValueError 
             << "could not read strain: strainSectionTest strain1? "
                "strain2? ... strainN?"
             << "\n";
      return TCL_ERROR;
    }
    data(i - 1) = strain;
  }

  if (theMaterial->setTrialStrain(data) != 0) {
    opserr << OpenSees::PromptValueError 
           << "failed to set trial strain\n";
    return TCL_ERROR;
  }

  return TCL_OK;
}


static int
MaterialTest_Commit(ClientData clientData, Tcl_Interp *interp,
                    Tcl_Size argc, TCL_Char ** const argv)
{
  NDMaterial *theMaterial = (NDMaterial*)clientData;
  const int status = theMaterial->commitState();
  return TCL_OK;
}


static int
MaterialTest_getStressSection(ClientData clientData, Tcl_Interp *interp,
                             Tcl_Size argc, TCL_Char ** const argv)
{
  NDMaterial *theMaterial = (NDMaterial*)clientData;
  const Vector &stress = theMaterial->getStress();
  const int nsr = stress.Size();
  Tcl_Obj* list = Tcl_NewListObj(nsr, nullptr);
  for (int i = 0; i < nsr; ++i) {
    Tcl_ListObjAppendElement(interp, list, Tcl_NewDoubleObj(stress(i)));
  }
  Tcl_SetObjResult(interp, list);
  return TCL_OK;
}


static int
MaterialTest_getTangent(ClientData clientData, Tcl_Interp *interp,
                        Tcl_Size argc, TCL_Char ** const argv)
{
  NDMaterial *theMaterial = (NDMaterial*)clientData;

  const Matrix &tangent = theMaterial->getTangent();
  const int nr = tangent.noRows();
  const int nc = tangent.noCols();
  Tcl_Obj* list = Tcl_NewListObj(nr * nc, nullptr);
  for (int i = 0; i < nr; ++i)
    for (int j = 0; j < nc; j++) {
      Tcl_ListObjAppendElement(interp, list, Tcl_NewDoubleObj(tangent(i, j)));
    }
  Tcl_SetObjResult(interp, list);
  return TCL_OK;
}


static int
MaterialTest_getResponse(ClientData clientData, Tcl_Interp *interp,
                                Tcl_Size argc, TCL_Char ** const argv)
{
  NDMaterial *theMaterial = (NDMaterial*)clientData;
  DummyStream dummy;
  Response *theResponse =
      theMaterial->setResponse(argv + 1, argc - 1, dummy);

  if (theResponse == nullptr) {
    opserr << OpenSees::PromptValueError 
           << "Response returned a null pointer\n";
    return TCL_ERROR;
  }

  if (theResponse->getResponse() < 0) {
    delete theResponse;
    opserr << OpenSees::PromptValueError << "Failed to get response\n";
    return TCL_ERROR;
  }

  Information &eleInfo = theResponse->getInformation();
  const Vector &data = eleInfo.getData();
  const int ni = data.Size();
  Tcl_Obj* list = Tcl_NewListObj(ni, nullptr);

  for (int i = 0; i < ni; ++i) {
    Tcl_ListObjAppendElement(interp, list, Tcl_NewDoubleObj(data(i)));
  }

  Tcl_SetObjResult(interp, list);
  delete theResponse;
  return TCL_OK;
}
