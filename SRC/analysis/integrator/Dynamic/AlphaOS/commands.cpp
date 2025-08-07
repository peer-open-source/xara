#include "Parsing.h"
#include "Logging.h"
#include "BasicAnalysisBuilder.h"

#include "elementAPI.h"
#include "AlphaOSGeneralized.h"
#include "AlphaOSGeneralized_TP.h"

int
TclCommand_createAlphaOSGeneralized(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);

  argc -= 2;

  if (argc != 1 && argc != 2 && argc != 4 && argc != 5) {
    opserr << "WARNING - incorrect number of args want AlphaOSGeneralized "
              "$rhoInf <-updateElemDisp>\n";
    opserr << "          or AlphaOSGeneralized $alphaI $alphaF $beta $gamma "
              "<-updateElemDisp>\n";
    return TCL_ERROR;
  }

  bool updElemDisp = false;
  double dData[4];
  int numData;
  if (argc < 3)
    numData = 1;
  else
    numData = 4;

  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING - invalid args want AlphaOSGeneralized $alpha "
              "<-updateElemDisp>\n";
    opserr << "          or AlphaOSGeneralized $alphaI $alphaF $beta $gamma "
              "<-updateElemDisp>\n";
    return TCL_ERROR;
  }

  if (argc == 2 || argc == 5) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-updateElemDisp") == 0)
      updElemDisp = true;
  }

  TransientIntegrator *theIntegrator = nullptr;
  if (argc < 3)
    theIntegrator = new AlphaOSGeneralized(dData[0], updElemDisp);
  else
    theIntegrator = new AlphaOSGeneralized(dData[0], dData[1], dData[2],
                                           dData[3], updElemDisp);

  builder->set(*theIntegrator);
  return TCL_OK;
}


int
TclCommand_createAlphaOSGeneralized_TP(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);

  argc -= 2;

  if (argc != 1 && argc != 2 && argc != 4 && argc != 5) {
    opserr << "WARNING - incorrect number of args want AlphaOSGeneralized_TP "
              "$rhoInf <-updateElemDisp>\n";
    opserr << "          or AlphaOSGeneralized_TP $alphaI $alphaF $beta $gamma "
              "<-updateElemDisp>\n";
    return TCL_ERROR;
  }

  bool updElemDisp = false;
  double dData[4];
  int numData;
  if (argc < 3)
    numData = 1;
  else
    numData = 4;

  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING - invalid args want AlphaOSGeneralized_TP $alpha "
              "<-updateElemDisp>\n";
    opserr << "          or AlphaOSGeneralized_TP $alphaI $alphaF $beta $gamma "
              "<-updateElemDisp>\n";
    return TCL_ERROR;
  }

  if (argc == 2 || argc == 5) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-updateElemDisp") == 0)
      updElemDisp = true;
  }

  TransientIntegrator *theIntegrator = nullptr;
  if (argc < 3)
    theIntegrator = new AlphaOSGeneralized_TP(dData[0], updElemDisp);
  else
    theIntegrator = new AlphaOSGeneralized_TP(dData[0], dData[1], dData[2],
                                              dData[3], updElemDisp);

  builder->set(*theIntegrator);
  return TCL_OK;
}




int
TclCommand_createAlphaOS(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv) 
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);

  argc -= 2;

  if (argc < 1 || argc > 4) {
    opserr << "WARNING - incorrect number of args want AlphaOS $alpha "
              "<-updateElemDisp>\n";
    opserr << "          or AlphaOS $alpha $beta $gamma <-updateElemDisp>\n";
    return TCL_ERROR;
  }

  bool updElemDisp = false;
  double dData[3];
  int numData;
  if (argc < 3)
    numData = 1;
  else
    numData = 3;

  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING - invalid args want AlphaOS $alpha <-updateElemDisp>\n";
    opserr << "          or AlphaOS $alpha $beta $gamma <-updateElemDisp>\n";
    return TCL_ERROR;
  }

  if (argc == 2 || argc == 4) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-updateElemDisp") == 0)
      updElemDisp = true;
  }

  TransientIntegrator *theIntegrator = nullptr;
  if (argc < 3)
    theIntegrator = new AlphaOS(dData[0], updElemDisp);
  else
    theIntegrator = new AlphaOS(dData[0], dData[1], dData[2], updElemDisp);


  builder->set(*theIntegrator);
  return TCL_OK;
}

 
int
TclCommand_createAlphaOS_TP(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);
  TransientIntegrator *theIntegrator = nullptr;

  argc -= 2;

  if (argc < 1 || argc > 4) {
    opserr << "WARNING - incorrect number of args want AlphaOS_TP $alpha "
              "<-updateElemDisp>\n";
    opserr << "          or AlphaOS_TP $alpha $beta $gamma <-updateElemDisp>\n";
    return TCL_ERROR;
  }

  bool updElemDisp = false;
  double dData[3];
  int numData;
  if (argc < 3)
    numData = 1;
  else
    numData = 3;

  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr
        << "WARNING - invalid args want AlphaOS_TP $alpha <-updateElemDisp>\n";
    opserr << "          or AlphaOS_TP $alpha $beta $gamma <-updateElemDisp>\n";
    return TCL_ERROR;
  }

  if (argc == 2 || argc == 4) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-updateElemDisp") == 0)
      updElemDisp = true;
  }

  if (argc < 3)
    theIntegrator = new AlphaOS_TP(dData[0], updElemDisp);
  else
    theIntegrator = new AlphaOS_TP(dData[0], dData[1], dData[2], updElemDisp);


  builder->set(*theIntegrator);
  return TCL_OK;
}


