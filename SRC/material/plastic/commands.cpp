
#include <tcl.h>
#include <string.h>
#include <Logging.h>
#include <Parsing.h>
#include <isotropy.h>
#include <BasicModelBuilder.h>
// #include <PlasticMaterial.h>

#include <J2Plasticity.h>
#include <MultiaxialCyclicPlasticity.h>

Tcl_CmdProc TclCommand_setIsotropicParameters;

int
TclCommand_newJ2Material(ClientData clientData,
                         Tcl_Interp* interp,
                         int argc,
                         const char** const argv)
{
  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);

  if ((strcmp(argv[1], "J2Plasticity") == 0) || 
      (strcmp(argv[1], "J2") == 0)) {

    if (argc < 9) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial J2Plasticity tag? K? G? sig0? sigInf? delta? H? <eta?>" << "\n";
      return TCL_ERROR;
    }

    int tag;
    double K, G, sig0, sigInf, delta, H;
    double eta = 0.0;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid J2Plasticity tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &K) != TCL_OK) {
      opserr << "WARNING invalid K\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &G) != TCL_OK) {
      opserr << "WARNING invalid G\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &sig0) != TCL_OK) {
      opserr << "WARNING invalid sig0\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &sigInf) != TCL_OK) {
      opserr << "WARNING invalid sigInf\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[7], &delta) != TCL_OK) {
      opserr << "WARNING invalid delta\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[8], &H) != TCL_OK) {
      opserr << "WARNING invalid H\n";
      return TCL_ERROR;
    }
    if (argc > 9 && Tcl_GetDouble(interp, argv[9], &eta) != TCL_OK) {
      opserr << "WARNING invalid eta\n";
      return TCL_ERROR;
    }

    builder->addTaggedObject<NDMaterial>(*new J2Plasticity(tag, 0, K, G, sig0, sigInf, delta, H, eta));
  }
}

#if 0
int
TclCommand_newPlasticMaterial(ClientData clientData, Tcl_Interp* interp, int argc, const char**const argv)
{
  if ((strcmp(argv[1], "J2Plasticity") == 0) || (strcmp(argv[1], "J2") == 0))
  {
    if (argc < 9) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: nDMaterial J2Plasticity tag? K? G? sig0? sigInf? delta? H? <eta?>" << "\n";
      return TCL_ERROR;
    }

    int tag;
    double K, G, sig0, sigInf, delta, H;
    double eta = 0.0;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid J2Plasticity tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[3], &K) != TCL_OK) {
      opserr << "WARNING invalid K\n";
      opserr << "nDMaterial J2Plasticity: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[4], &G) != TCL_OK) {
      opserr << "WARNING invalid G\n";
      opserr << "nDMaterial J2Plasticity: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[5], &sig0) != TCL_OK) {
      opserr << "WARNING invalid sig0\n";
      opserr << "nDMaterial J2Plasticity: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &sigInf) != TCL_OK) {
      opserr << "WARNING invalid sigInf\n";
      opserr << "nDMaterial J2Plasticity: " << tag << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[7], &delta) != TCL_OK) {
      opserr << "WARNING invalid delta\n";
      opserr << "nDMaterial J2Plasticity: " << tag << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[8], &H) != TCL_OK) {
      opserr << "WARNING invalid H\n";
      opserr << "nDMaterial J2Plasticity: " << tag << "\n";
      return TCL_ERROR;
    }
    if (argc > 9 && Tcl_GetDouble(interp, argv[9], &eta) != TCL_OK) {
      opserr << "WARNING invalid eta\n";
      opserr << "nDMaterial J2Plasticity: " << tag << "\n";
      return TCL_ERROR;
    }

    theMaterial = new J2Plasticity(tag, 0, K, G, sig0, sigInf, delta, H, eta);
  }
}
#endif
