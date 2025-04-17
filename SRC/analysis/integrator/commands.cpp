//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
#include <Parsing.h>
#include <Logging.h>
#include <Node.h>
#include <Domain.h>
#include <BasicAnalysisBuilder.h>
#include <Dynamic/Houbolt.h>

#if 0
int
TclCommand_createHoubolt(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);
  return new Houbolt(); 
}
#endif

#include <Dynamic/AlphaOSGeneralized.h>
#include <Dynamic/AlphaOSGeneralized_TP.h>
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




#include <Dynamic/AlphaOS.h>
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




#include <analysis/integrator/Dynamic/AlphaOS_TP.h>
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




#include <analysis/integrator/Static/ArcLength1.h>
int
TclCommand_createArcLength1(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv) 
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);

  double arcLength;
  double alpha;
  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "WARNING integrator ArcLength arcLength alpha \n";
    return TCL_ERROR;
  }

  int numdata = 1;
  if (OPS_GetDoubleInput(&numdata, &arcLength) < 0) {
    opserr << "WARNING integrator ArcLength failed to read arc length\n";
    return TCL_ERROR;
  }
  if (OPS_GetDoubleInput(&numdata, &alpha) < 0) {
    opserr << "WARNING integrator ArcLength failed to read alpha\n";
    return TCL_ERROR;
  }

  builder->set(*new ArcLength1(arcLength, alpha));
  return TCL_OK;
}


#include <analysis/integrator/Static/ArcLength.h>
int
TclCommand_createArcLength(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv) 
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);

  double arcLength;
  double alpha;
  if (OPS_GetNumRemainingInputArgs() < 2) {
    opserr << "WARNING integrator ArcLength arcLength alpha \n";
    return TCL_ERROR;
  }

  int numdata = 1;
  if (OPS_GetDoubleInput(&numdata, &arcLength) < 0) {
    opserr << "WARNING integrator ArcLength failed to read arc lenght\n";
    return TCL_ERROR;
  }
  if (OPS_GetDoubleInput(&numdata, &alpha) < 0) {
    opserr << "WARNING integrator ArcLength failed to read alpha\n";
    return TCL_ERROR;
  }
  builder->set(*new ArcLength(arcLength, alpha));
  return TCL_OK;
}




#include <analysis/integrator/Dynamic/BackwardEuler.h>
int
TclCommand_createBackwardEuler(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv) 
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);

  int optn = 0;
  if (OPS_GetNumRemainingInputArgs() > 0) {
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &optn) < 0) {
      opserr << "WARNING integrator BackwardEuler <option> - undefined option "
                "specified\n";
      return TCL_ERROR;
    }
  }
  builder->set(* new BackwardEuler(optn));
  return TCL_OK;
}



#if 0
#include <analysis/integrator/CollocationHSFixedNumIter.h>
int
TclCommand_createCollocationHSFixedNumIter(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv) 
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);

  TransientIntegrator *theIntegrator = 0;

  argc -= 2;

  if (argc != 1 && argc != 3 && argc != 5) {
    opserr << "WARNING - incorrect number of args want "
              "CollocationHSFixedNumIter $theta <-polyOrder $O>\n";
    opserr << "          or CollocationHSFixedNumIter $theta $beta $gamma "
              "<-polyOrder $O>\n";
    return TCL_ERROR;
  }

  double dData[3];
  int polyOrder = 2;
  int numData = 0;

  // count number of numeric parameters
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-polyOrder") == 0) {
      break;
    }
    numData++;
  }
  // reset to read from beginning
  OPS_ResetCurrentInputArg(2);

  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING - invalid args want CollocationHSFixedNumIter $theta "
              "<-polyOrder $O>\n";
    opserr << "          or CollocationHSFixedNumIter $theta $beta $gamma "
              "<-polyOrder $O>\n";
    return TCL_ERROR;
  }

  if (numData + 2 == argc) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-polyOrder") == 0) {
      int numData2 = 1;
      if (OPS_GetInt(&numData2, &polyOrder) != 0) {
        opserr << "WARNING - invalid polyOrder want CollocationHSFixedNumIter "
                  "$rhoInf <-polyOrder $O>\n";
        opserr << "          or CollocationHSFixedNumIter $alphaI $alphaF "
                  "$beta $gamma <-polyOrder $O>\n";
      }
    }
  }

  if (numData == 1)
    theIntegrator = new CollocationHSFixedNumIter(dData[0], polyOrder);
  else if (numData == 3)
    theIntegrator =
        new CollocationHSFixedNumIter(dData[0], dData[1], dData[2], polyOrder);

  builder->set(*theIntegrator);
  return TCL_OK;
}





#include <analysis/integrator/CollocationHSIncrLimit.h>
int
TclCommand_createCollocationHSIncrLimit(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);

  argc -= 2;

  if (argc != 2 && argc != 4 && argc != 6) {
    opserr << "WARNING - incorrect number of args want CollocationHSIncrLimit "
              "$theta $limit <-normType $T>\n";
    opserr << "          or CollocationHSIncrLimit $theta $beta $gamma $limit "
              "<-normType $T>\n";
    return TCL_ERROR;
  }

  double dData[4];
  int normType = 2;
  int numData = 0;

  // count number of numeric parameters
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-normType") == 0) {
      break;
    }
    numData++;
  }
  // reset to read from beginning
  OPS_ResetCurrentInputArg(2);

  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING - invalid args want CollocationHSIncrLimit $theta "
              "$limit <-normType $T>\n";
    opserr << "          or CollocationHSIncrLimit $theta $beta $gamma $limit "
              "<-normType $T>\n";
    return TCL_ERROR;
  }

  if (numData + 2 == argc) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-normType") == 0) {
      int numData2 = 1;
      if (OPS_GetInt(&numData2, &normType) != 0) {
        opserr << "WARNING - invalid normType want CollocationHSIncrLimit "
                  "$theta $limit <-normType $T>\n";
        opserr << "          or CollocationHSIncrLimit $theta $beta $gamma "
                  "$limit <-normType $T>\n";
      }
    }
  }

  TransientIntegrator *theIntegrator = 0;
  if (numData == 2)
    theIntegrator = new CollocationHSIncrLimit(dData[0], dData[1], normType);
  else if (numData == 4)
    theIntegrator = new CollocationHSIncrLimit(dData[0], dData[1], dData[2],
                                               dData[3], normType);


  builder->set(*theIntegrator);
  return TCL_OK;
}


#include <analysis/integrator/CollocationHSIncrReduct.h>

int
TclCommand_createCollocationHSIncrReduct(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);

  TransientIntegrator *theIntegrator = 0;

  argc -= 2;

  if (argc != 2 && argc != 4) {
    opserr << "WARNING - incorrect number of args want CollocationHSIncrReduct "
              "$theta $reduct\n";
    opserr
        << "          or CollocationHSIncrReduct $theta $beta $gamma $reduct\n";
    return TCL_ERROR;
  }

  double dData[4];
  if (OPS_GetDouble(&argc, dData) != 0) {
    opserr << "WARNING - invalid args want CollocationHSIncrReduct $theta "
              "$reduct\n";
    opserr
        << "          or CollocationHSIncrReduct $theta $beta $gamma $reduct\n";
    return TCL_ERROR;
  }

  if (argc == 2)
    theIntegrator = new CollocationHSIncrReduct(dData[0], dData[1]);
  else
    theIntegrator =
        new CollocationHSIncrReduct(dData[0], dData[1], dData[2], dData[3]);

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating CollocationHSIncrReduct "
              "integrator\n";

  builder->set(*theIntegrator);
  return TCL_OK;
}




#include <analysis/integrator/Collocation.h>
int
TclCommand_createCollocation(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);
  TransientIntegrator *theIntegrator = 0;

  argc -= 2;

  if (argc != 1 && argc != 3) {
    opserr << "WARNING - incorrect number of args want Collocation $theta\n";
    opserr << "          or Collocation $theta $beta $gamma\n";
    return TCL_ERROR;
  }

  double dData[3];
  if (OPS_GetDouble(&argc, dData) != 0) {
    opserr << "WARNING - invalid args want Collocation $theta\n";
    opserr << "          or Collocation $theta $beta $gamma\n";
    return TCL_ERROR;
  }

  if (argc == 1)
    theIntegrator = new Collocation(dData[0]);
  else
    theIntegrator = new Collocation(dData[0], dData[1], dData[2]);

  builder->set(*theIntegrator);
  return TCL_OK;
}
#endif




#include <analysis/integrator/Static/DisplacementControl.h>
int
TclCommand_createDisplacementControlIntegrator(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv) 
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);

  if (OPS_GetNumRemainingInputArgs() < 3) {
    opserr << "insufficient arguments for DisplacementControl\n";
    return TCL_ERROR;
  }

  // node, dof
  int iData[2];
  int numData = 2;
  if (OPS_GetIntInput(&numData, &iData[0]) < 0) {
    opserr << "WARNING failed to read node tag and ndf\n";
    return TCL_ERROR;
  }

  double incr;
  numData = 1;
  if (OPS_GetDoubleInput(&numData, &incr) < 0) {
    opserr << "WARNING failed to read incr\n";
    return TCL_ERROR;
  }

  // numIter,dumin,dumax
  int numIter = 1;
  int formTangent = 0;
  double data[2] = {incr, incr};
  if (OPS_GetNumRemainingInputArgs() > 2) {
    numData = 1;
    if (OPS_GetIntInput(&numData, &numIter) < 0) {
      opserr << "WARNING failed to read numIter\n";
      return TCL_ERROR;
    }
    numData = 2;
    if (OPS_GetDoubleInput(&numData, &data[0]) < 0) {
      opserr << "WARNING failed to read dumin and dumax\n";
      return TCL_ERROR;
    }
  }

  if (OPS_GetNumRemainingInputArgs() == 1) {
    std::string type = OPS_GetString();
    if (type == "-initial" || type == "-Initial") {
      formTangent = 1;
    }
  }

  // check node
  Domain *theDomain = builder->getDomain();
  Node *theNode = theDomain->getNode(iData[0]);
  if (theNode == nullptr) {
    opserr << "WARNING integrator DisplacementControl node dof dU : Node does "
              "not exist\n";
    return TCL_ERROR;
  }

  int numDOF = theNode->getNumberDOF();
  if (iData[1] <= 0 || iData[1] > numDOF) {
    opserr << "WARNING integrator DisplacementControl node dof dU : invalid "
              "dof given\n";
    return TCL_ERROR;
  }

  builder->set(*new DisplacementControl(iData[0], iData[1] - 1, incr, theDomain,
                                 numIter, data[0], data[1], formTangent));
  return TCL_OK;
}



#include <analysis/integrator/Dynamic/GeneralizedAlpha.h>
int
TclCommand_createGeneralizedAlpha(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{

  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);
  TransientIntegrator *theIntegrator = nullptr;

  argc -= 2;

  if (argc != 2 && argc != 4) {
    opserr << "WARNING - incorrect number of args want GeneralizedAlpha "
              "$alphaM $alphaF <$gamma $beta>\n";
    return TCL_ERROR;
  }

  double dData[4];
  if (OPS_GetDouble(&argc, dData) != 0) {
    opserr << "WARNING - invalid args want GeneralizedAlpha $alphaM $alphaF "
              "<$gamma $beta>\n";
    return TCL_ERROR;
  }

  if (argc == 2)
    theIntegrator = new GeneralizedAlpha(dData[0], dData[1]);
  else
    theIntegrator =
        new GeneralizedAlpha(dData[0], dData[1], dData[2], dData[3]);

  builder->set(*theIntegrator);
  return TCL_OK;
}




#if 0
#include <analysis/integrator/GimmeMCK.h>
int
TclCommand_createGimmeMCK(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);
  TransientIntegrator *theIntegrator = 0;

  argc -= 2;

  if (argc < 3) {
    opserr
        << "WARNING - incorrect number of args want GimmeMCK $m $c $k <$ki>\n";
    return TCL_ERROR;
  }

  int numdata = 3;
  double ddata[3];
  if (OPS_GetDouble(&numdata, ddata) != 0) {
    opserr << "WARNING - invalid args want GimmeMCK $m $c $k <$ki>\n";
    return TCL_ERROR;
  }
  numdata = 1;
  double ki = 0.0;
  if (argc > 3) {
    if (OPS_GetDouble(&numdata, &ki) != 0) {
      opserr << "WARNING - invalid args want GimmeMCK $m $c $k <$ki>\n";
      return TCL_ERROR;
    }
  }

  theIntegrator = new GimmeMCK(ddata[0], ddata[1], ddata[2], ki);


  builder->set(*theIntegrator);
  return TCL_OK;
}


#include <analysis/integrator/Dynamic/HHTExplicit.h>
int
TclCommand_createHHTExplicit(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);
  TransientIntegrator *theIntegrator = 0;

  argc -= 2;

  if (argc < 1 || argc > 3) {
    opserr << "WARNING - incorrect number of args want HHTExplicit $alpha "
              "<-updateElemDisp>\n";
    opserr << "          or HHTExplicit $alpha $gamma <-updateElemDisp>\n";
    return TCL_ERROR;
  }

  bool updElemDisp = false;
  double dData[2];
  int numData = 0;

  // count number of numeric parameters
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-updateElemDisp") == 0) {
      break;
    }
    numData++;
  }
  // reset to read from beginning
  OPS_ResetCurrentInputArg(2);

  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr
        << "WARNING - invalid args want HHTExplicit $alpha <-updateElemDisp>\n";
    opserr << "          or HHTExplicit $alpha $gamma <-updateElemDisp>\n";
    return TCL_ERROR;
  }

  if (numData + 1 == argc) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-updateElemDisp") == 0)
      updElemDisp = true;
  }

  if (numData == 1)
    theIntegrator = new HHTExplicit(dData[0], updElemDisp);
  else if (numData == 2)
    theIntegrator = new HHTExplicit(dData[0], dData[1], updElemDisp);

  builder->set(*theIntegrator);
  return TCL_OK;
}




#include <analysis/integrator/Dynamic/HHTExplicit_TP.h>
int
TclCommand_createHHTExplicit_TP(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);

  argc -= 2;

  if (argc < 1 || argc > 2) {
    opserr << "WARNING - incorrect number of args want HHTExplicit_TP $alpha\n";
    opserr << "          or HHTExplicit_TP $alpha $gamma\n";
    return TCL_ERROR;
  }

  double dData[2];
  if (OPS_GetDouble(&argc, dData) != 0) {
    opserr << "WARNING - invalid args want HHTExplicit_TP $alpha\n";
    opserr << "          or HHTExplicit_TP $alpha $gamma\n";
    return TCL_ERROR;
  }

  TransientIntegrator *theIntegrator = 0;
  if (argc == 1)
    theIntegrator = new HHTExplicit_TP(dData[0]);
  else if (argc == 2)
    theIntegrator = new HHTExplicit_TP(dData[0], dData[1]);

  builder->set(*theIntegrator);
  return TCL_OK;
}





#include <analysis/integrator/Dynamic/HHTGeneralizedExplicit.h>
int
TclCommand_createHHTGeneralizedExplicit(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);
  TransientIntegrator *theIntegrator = nullptr;

  argc -= 2;

  if (argc < 2 || argc > 5) {
    opserr << "WARNING - incorrect number of args want HHTGeneralizedExplicit "
              "$rhoB $alphaF <-updateElemDisp>\n";
    opserr << "          or HHTGeneralizedExplicit $alphaI $alphaF $beta "
              "$gamma <-updateElemDisp>\n";
    return TCL_ERROR;
  }

  bool updElemDisp = false;
  double dData[4];
  int numData;
  if (argc < 4)
    numData = 2;
  else
    numData = 4;

  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING - invalid args want HHTGeneralizedExplicit $rhoB "
              "$alphaF <-updateElemDisp>\n";
    opserr << "          or HHTGeneralizedExplicit $alphaI $alphaF $beta "
              "$gamma <-updateElemDisp>\n";
    return TCL_ERROR;
  }

  if (argc == 3 || argc == 5) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-updateElemDisp") == 0)
      updElemDisp = true;
  }

  if (argc < 4)
    theIntegrator = new HHTGeneralizedExplicit(dData[0], dData[1], updElemDisp);
  else
    theIntegrator = new HHTGeneralizedExplicit(dData[0], dData[1], dData[2],
                                               dData[3], updElemDisp);

  builder->set(*theIntegrator);
  return TCL_OK;
}




#include <analysis/integrator/Dynamic/HHTGeneralizedExplicit_TP.h>
int
TclCommand_createHHTGeneralizedExplicit_TP(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  TransientIntegrator *theIntegrator = 0;

  argc -= 2;

  if (argc != 2 && argc != 4) {
    opserr << "WARNING - incorrect number of args want "
              "HHTGeneralizedExplicit_TP $rhoB $alphaF\n";
    opserr << "          or HHTGeneralizedExplicit_TP $alphaI $alphaF $beta "
              "$gamma\n";
    return TCL_ERROR;
  }

  double dData[4];
  if (OPS_GetDouble(&argc, dData) != 0) {
    opserr << "WARNING - invalid args want HHTGeneralizedExplicit_TP $rhoB "
              "$alphaF\n";
    opserr << "          or HHTGeneralizedExplicit_TP $alphaI $alphaF $beta "
              "$gamma\n";
    return TCL_ERROR;
  }

  if (argc == 2)
    theIntegrator = new HHTGeneralizedExplicit_TP(dData[0], dData[1]);
  else if (argc == 4)
    theIntegrator =
        new HHTGeneralizedExplicit_TP(dData[0], dData[1], dData[2], dData[3]);

  builder->set(*theIntegrator);
  return TCL_OK;
}




#include <analysis/integrator/Dynamic/HHTGeneralized.h>
int
TclCommand_createHHTGeneralized(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  TransientIntegrator *theIntegrator = 0;

  argc -= 2;

  if (argc != 1 && argc != 4) {
    opserr
        << "WARNING - incorrect number of args want HHTGeneralized $rhoInf\n";
    opserr << "          or HHTGeneralized $alphaI $alphaF $beta $gamma\n";
    return TCL_ERROR;
  }

  double dData[4];
  if (OPS_GetDouble(&argc, dData) != 0) {
    opserr << "WARNING - invalid args want HHTGeneralized $rhoInf\n";
    opserr << "          or HHTGeneralized $alphaI $alphaF $beta $gamma\n";
    return TCL_ERROR;
  }

  if (argc == 1)
    theIntegrator = new HHTGeneralized(dData[0]);
  else
    theIntegrator = new HHTGeneralized(dData[0], dData[1], dData[2], dData[3]);

  builder->set(*theIntegrator);
  return TCL_OK;
}




#include <analysis/integrator/HHTGeneralized_TP.h>
int
TclCommand_createHHTGeneralized_TP(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  TransientIntegrator *theIntegrator = 0;

  argc -= 2;

  if (argc != 1 && argc != 4) {
    opserr << "WARNING - incorrect number of args want HHTGeneralized_TP "
              "$rhoInf\n";
    opserr << "          or HHTGeneralized_TP $alphaI $alphaF $beta $gamma\n";
    return TCL_ERROR;
  }

  double dData[4];
  if (OPS_GetDouble(&argc, dData) != 0) {
    opserr << "WARNING - invalid args want HHTGeneralized_TP $rhoInf\n";
    opserr << "          or HHTGeneralized_TP $alphaI $alphaF $beta $gamma\n";
    return TCL_ERROR;
  }

  if (argc == 1)
    theIntegrator = new HHTGeneralized_TP(dData[0]);
  else
    theIntegrator =
        new HHTGeneralized_TP(dData[0], dData[1], dData[2], dData[3]);

  builder->set(*theIntegrator);
  return TCL_OK;
}




#include <analysis/integrator/Dynamic/HHTHSFixedNumIter.h>
int
TclCommand_createHHTHSFixedNumIter(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  TransientIntegrator *theIntegrator = 0;

  argc -= 2;

  if (argc != 1 && argc != 3 && argc != 4 && argc != 6) {
    opserr << "WARNING - incorrect number of args want HHTHSFixedNumIter "
              "$rhoInf <-polyOrder $O>\n";
    opserr << "          or HHTHSFixedNumIter $alphaI $alphaF $beta $gamma "
              "<-polyOrder $O>\n";
    return TCL_ERROR;
  }

  double dData[4];
  int polyOrder = 2;
  bool updDomFlag = true;
  int numData;
  if (argc < 4)
    numData = 1;
  else
    numData = 4;

  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING - invalid args want HHTHSFixedNumIter $rhoInf "
              "<-polyOrder $O>\n";
    opserr << "          or HHTHSFixedNumIter $alphaI $alphaF $beta $gamma "
              "<-polyOrder $O>\n";
    return TCL_ERROR;
  }

  if (argc == 3 || argc == 6) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-polyOrder") == 0) {
      numData = 1;
      if (OPS_GetInt(&numData, &polyOrder) != 0) {
        opserr << "WARNING - invalid polyOrder want HHTHSFixedNumIter $rhoInf "
                  "<-polyOrder $O>\n";
        opserr << "          or HHTHSFixedNumIter $alphaI $alphaF $beta $gamma "
                  "<-polyOrder $O>\n";
      }
    }
  }

  if (argc < 4)
    theIntegrator = new HHTHSFixedNumIter(dData[0], polyOrder, updDomFlag);
  else
    theIntegrator = new HHTHSFixedNumIter(dData[0], dData[1], dData[2],
                                          dData[3], polyOrder, updDomFlag);

  builder->set(*theIntegrator);
  return TCL_OK;
}




#include <analysis/integrator/Dynamic/HHTHSFixedNumIter_TP.h>
int
TclCommand_createHHTHSFixedNumIter_TP(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  TransientIntegrator *theIntegrator = 0;

  argc -= 2;

  if (argc != 1 && argc != 3 && argc != 4 && argc != 6) {
    opserr << "WARNING - incorrect number of args want HHTHSFixedNumIter_TP "
              "$rhoInf <-polyOrder $O>\n";
    opserr << "          or HHTHSFixedNumIter_TP $alphaI $alphaF $beta $gamma "
              "<-polyOrder $O>\n";
    return TCL_ERROR;
  }

  double dData[4];
  int polyOrder = 2;
  bool updDomFlag = true;
  int numData;
  if (argc < 4)
    numData = 1;
  else
    numData = 4;

  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING - invalid args want HHTHSFixedNumIter_TP $rhoInf "
              "<-polyOrder $O>\n";
    opserr << "          or HHTHSFixedNumIter_TP $alphaI $alphaF $beta $gamma "
              "<-polyOrder $O>\n";
    return TCL_ERROR;
  }

  if (argc == 3 || argc == 6) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-polyOrder") == 0) {
      numData = 1;
      if (OPS_GetInt(&numData, &polyOrder) != 0) {
        opserr << "WARNING - invalid polyOrder want HHTHSFixedNumIter_TP "
                  "$rhoInf <-polyOrder $O>\n";
        opserr << "          or HHTHSFixedNumIter_TP $alphaI $alphaF $beta "
                  "$gamma <-polyOrder $O>\n";
      }
    }
  }

  if (argc < 4)
    theIntegrator = new HHTHSFixedNumIter_TP(dData[0], polyOrder, updDomFlag);
  else
    theIntegrator = new HHTHSFixedNumIter_TP(dData[0], dData[1], dData[2],
                                             dData[3], polyOrder, updDomFlag);

  builder->set(*theIntegrator);
  return TCL_OK;
}




#include <analysis/integrator/Dynamic/HHTHSIncrLimit.h>
int
TclCommand_createHHTHSIncrLimit(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  TransientIntegrator *theIntegrator = 0;

  argc -= 2;

  if (argc != 2 && argc != 4 && argc != 5 && argc != 7) {
    opserr << "WARNING - incorrect number of args want HHTHSIncrLimit $rhoInf "
              "$limit <-normType $T>\n";
    opserr << "          or HHTHSIncrLimit $alphaI $alphaF $beta $gamma $limit "
              "<-normType $T>\n";
    return TCL_ERROR;
  }

  double dData[5];
  int normType = 2;
  int numData;
  if (argc < 5)
    numData = 2;
  else
    numData = 5;

  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING - invalid args want HHTHSIncrLimit $rhoInf $limit "
              "<-normType $T>\n";
    opserr << "          or HHTHSIncrLimit $alphaI $alphaF $beta $gamma $limit "
              "<-normType $T>\n";
    return TCL_ERROR;
  }

  if (argc == 4 || argc == 7) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-normType") == 0) {
      numData = 1;
      if (OPS_GetInt(&numData, &normType) != 0) {
        opserr << "WARNING - invalid normType want HHTHSIncrLimit $rhoInf "
                  "$limit <-normType $T>\n";
        opserr << "          or HHTHSIncrLimit $alphaI $alphaF $beta $gamma "
                  "$limit <-normType $T>\n";
      }
    }
  }

  if (argc < 5)
    theIntegrator = new HHTHSIncrLimit(dData[0], dData[1], normType);
  else
    theIntegrator = new HHTHSIncrLimit(dData[0], dData[1], dData[2], dData[3],
                                       dData[4], normType);

  builder->set(*theIntegrator);
  return TCL_OK;
}




#include <analysis/integrator/Dynamic/HHTHSIncrLimit_TP.h>
int
TclCommand_createHHTHSIncrLimit_TP(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  TransientIntegrator *theIntegrator = 0;

  argc -= 2;

  if (argc != 2 && argc != 4 && argc != 5 && argc != 7) {
    opserr << "WARNING - incorrect number of args want HHTHSIncrLimit_TP "
              "$rhoInf $limit <-normType $T>\n";
    opserr << "          or HHTHSIncrLimit_TP $alphaI $alphaF $beta $gamma "
              "$limit <-normType $T>\n";
    return TCL_ERROR;
  }

  double dData[5];
  int normType = 2;
  int numData;
  if (argc < 5)
    numData = 2;
  else
    numData = 5;

  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING - invalid args want HHTHSIncrLimit_TP $rhoInf $limit "
              "<-normType $T>\n";
    opserr << "          or HHTHSIncrLimit_TP $alphaI $alphaF $beta $gamma "
              "$limit <-normType $T>\n";
    return TCL_ERROR;
  }

  if (argc == 4 || argc == 7) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-normType") == 0) {
      numData = 1;
      if (OPS_GetInt(&numData, &normType) != 0) {
        opserr << "WARNING - invalid normType want HHTHSIncrLimit_TP $rhoInf "
                  "$limit <-normType $T>\n";
        opserr << "          or HHTHSIncrLimit_TP $alphaI $alphaF $beta $gamma "
                  "$limit <-normType $T>\n";
      }
    }
  }

  if (argc < 5)
    theIntegrator = new HHTHSIncrLimit_TP(dData[0], dData[1], normType);
  else
    theIntegrator = new HHTHSIncrLimit_TP(dData[0], dData[1], dData[2],
                                          dData[3], dData[4], normType);

  builder->set(*theIntegrator);
  return TCL_OK;
}




#include <analysis/integrator/Dynamic/HHTHSIncrReduct.h>
int
TclCommand_createHHTHSIncrReduct(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  TransientIntegrator *theIntegrator = 0;

  argc -= 2;

  if (argc != 2 && argc != 5) {
    opserr << "WARNING - incorrect number of args want HHTHSIncrReduct $rhoInf "
              "$reduct\n";
    opserr << "          or HHTHSIncrReduct $alphaI $alphaF $beta $gamma "
              "$reduct\n";
    return TCL_ERROR;
  }

  double dData[5];
  if (OPS_GetDouble(&argc, dData) != 0) {
    opserr << "WARNING - invalid args want HHTHSIncrReduct $rhoInf $reduct\n";
    opserr << "          or HHTHSIncrReduct $alphaI $alphaF $beta $gamma "
              "$reduct\n";
    return TCL_ERROR;
  }

  if (argc == 2)
    theIntegrator = new HHTHSIncrReduct(dData[0], dData[1]);
  else
    theIntegrator =
        new HHTHSIncrReduct(dData[0], dData[1], dData[2], dData[3], dData[4]);

  builder->set(*theIntegrator);
  return TCL_OK;
}




#include <analysis/integrator/Dynamic/HHTHSIncrReduct_TP.h>
int
TclCommand_createHHTHSIncrReduct_TP(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  TransientIntegrator *theIntegrator = 0;

  argc -= 2;

  if (argc != 2 && argc != 5) {
    opserr << "WARNING - incorrect number of args want HHTHSIncrReduct_TP "
              "$rhoInf $reduct\n";
    opserr << "          or HHTHSIncrReduct_TP $alphaI $alphaF $beta $gamma "
              "$reduct\n";
    return TCL_ERROR;
  }

  double dData[5];
  if (OPS_GetDouble(&argc, dData) != 0) {
    opserr
        << "WARNING - invalid args want HHTHSIncrReduct_TP $rhoInf $reduct\n";
    opserr << "          or HHTHSIncrReduct_TP $alphaI $alphaF $beta $gamma "
              "$reduct\n";
    return TCL_ERROR;
  }

  if (argc == 2)
    theIntegrator = new HHTHSIncrReduct_TP(dData[0], dData[1]);
  else
    theIntegrator = new HHTHSIncrReduct_TP(dData[0], dData[1], dData[2],
                                           dData[3], dData[4]);

  builder->set(*theIntegrator);
  return TCL_OK;
}
#endif

//
//
//
#include <analysis/integrator/Dynamic/HHT.h>

int
TclCommand_createHHT(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv) 
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);

  argc -= 2;

  if (argc != 1 && argc != 3) {
    opserr << "WARNING - incorrect number of args want HHT $alpha <$gamma "
              "$beta>\n";
    return TCL_ERROR;
  }

  double dData[3];
  if (OPS_GetDouble(&argc, dData) != 0) {
    opserr << "WARNING - invalid args want HHT $alpha <$gamma $beta>\n";
    return TCL_ERROR;
  }

  TransientIntegrator *theIntegrator = nullptr;
  if (argc == 1)
    theIntegrator = new HHT(dData[0]);
  else
    theIntegrator = new HHT(dData[0], dData[1], dData[2]);


  builder->set(*theIntegrator);
  return TCL_OK;
}




#include <analysis/integrator/Dynamic/HHT_TP.h>
int
TclCommand_createHHT_TP(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv) 
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);
  TransientIntegrator *theIntegrator = 0;

  argc -= 2;

  if (argc != 1 && argc != 3) {
    opserr << "WARNING - incorrect number of args want HHT_TP $alpha <$gamma "
              "$beta>\n";
    return TCL_ERROR;
  }

  double dData[3];
  if (OPS_GetDouble(&argc, dData) != 0) {
    opserr << "WARNING - invalid args want HHT_TP $alpha <$gamma $beta>\n";
    return TCL_ERROR;
  }

  if (argc == 1)
    theIntegrator = new HHT_TP(dData[0]);
  else
    theIntegrator = new HHT_TP(dData[0], dData[1], dData[2]);

  builder->set(*theIntegrator);
  return TCL_OK;
}



#if 0
#include <analysis/integrator/HSConstraint.h>
int
TclCommand_createHSConstraint() {
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 1) {
    opserr << "WARNING integrator HSConstraint <arcLength> <psi_u> <psi_f> "
              "<u_ref> \n";
    return TCL_ERROR;
  }
  if (numdata > 4)
    numdata = 4;

  double data[4];
  if (OPS_GetDoubleInput(&numdata, data) < 0) {
    opserr << "WARNING integrator HSConstraint invalid double inputs\n";
    return TCL_ERROR;
  }
  double arcLength = data[0];
  double psi_u = data[1];
  double psi_f = data[2];
  double u_ref = data[3];

  switch (numdata) {
  case 1:
    return new HSConstraint(arcLength);
  case 2:
    return new HSConstraint(arcLength, psi_u);
  case 3:
    return new HSConstraint(arcLength, psi_u, psi_f);
  case 4:
    return new HSConstraint(arcLength, psi_u, psi_f, u_ref);
  }

  return TCL_ERROR;
}
#endif




#include <analysis/integrator/Dynamic/KRAlphaExplicit.h>
int
TclCommand_createKRAlphaExplicit(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);
  TransientIntegrator *theIntegrator = 0;

  argc -= 2;

  if (argc != 1 && argc != 2) {
    opserr << "WARNING - incorrect number of args want KRAlphaExplicit $rhoInf "
              "<-updateElemDisp>\n";
    return TCL_ERROR;
  }

  bool updElemDisp = false;
  double rhoInf;
  int numData = 1;
  if (OPS_GetDouble(&numData, &rhoInf) != 0) {
    opserr << "WARNING - invalid args want KRAlphaExplicit $rhoInf "
              "<-updateElemDisp>\n";
    return TCL_ERROR;
  }

  if (argc == 2) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-updateElemDisp") == 0)
      updElemDisp = true;
  }

  theIntegrator = new KRAlphaExplicit(rhoInf, updElemDisp);

  builder->set(*theIntegrator);
  return TCL_OK;
}




#include <analysis/integrator/Dynamic/KRAlphaExplicit_TP.h>
int
TclCommand_createKRAlphaExplicit_TP(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);
  TransientIntegrator *theIntegrator = 0;

  argc -= 2;

  if (argc != 1) {
    opserr << "WARNING - incorrect number of args want KRAlphaExplicit_TP "
              "$rhoInf\n";
    return TCL_ERROR;
  }

  double rhoInf;
  if (OPS_GetDouble(&argc, &rhoInf) != 0) {
    opserr << "WARNING - invalid args want KRAlphaExplicit_TP $rhoInf\n";
    return TCL_ERROR;
  }

  theIntegrator = new KRAlphaExplicit_TP(rhoInf);

  builder->set(*theIntegrator);
  return TCL_OK;
}




#include <analysis/integrator/Static/LoadControl.h>
int
TclCommand_createLoadControlIntegrator(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);

  if (OPS_GetNumRemainingInputArgs() < 1) {
    opserr << "LoadControl - insufficient arguments\n";
    return TCL_ERROR;
  }

  double lambda;
  int numData = 1;
  if (OPS_GetDoubleInput(&numData, &lambda) < 0) {
    opserr << "WARNING LoadControl - failed to read double lambda\n";
    return TCL_ERROR;
  }

  int numIter = 1;
  double mLambda[2] = {lambda, lambda};
  if (OPS_GetNumRemainingInputArgs() > 2) {
    if (OPS_GetIntInput(&numData, &numIter) < 0) {
      opserr << "WARNING LoadControl - failed to read int numIter\n";
      return TCL_ERROR;
    }
    numData = 2;
    if (OPS_GetDoubleInput(&numData, &mLambda[0]) < 0) {
      opserr << "WARNING LoadControl - failed to read double min and max\n";
      return TCL_ERROR;
    }
  }

  builder->set(*new LoadControl(lambda, numIter, mLambda[0], mLambda[1]));
  return TCL_OK;
}


#include <analysis/integrator/Static/MinUnbalDispNorm.h>
int
TclCommand_createMinUnbalDispNorm(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);

  double lambda11, minlambda, maxlambda;
  int numIter;
  if (OPS_GetNumRemainingInputArgs() < 1) {
    opserr << "WARNING integrator MinUnbalDispNorm lambda11 <Jd minLambda1j "
              "maxLambda1j>\n";
    return TCL_ERROR;
  }

  int numdata = 1;
  if (OPS_GetDoubleInput(&numdata, &lambda11) < 0) {
    opserr << "WARNING integrator MinUnbalDispNorm invalid lambda11\n";
    return TCL_ERROR;
  }

  if (OPS_GetNumRemainingInputArgs() >= 3) {
    if (OPS_GetIntInput(&numdata, &numIter) < 0) {
      opserr << "WARNING integrator MinUnbalDispNorm invalid numIter\n";
      return TCL_ERROR;
    }
    if (OPS_GetDoubleInput(&numdata, &minlambda) < 0) {
      opserr << "WARNING integrator MinUnbalDispNorm invalid minlambda\n";
      return TCL_ERROR;
    }
    if (OPS_GetDoubleInput(&numdata, &maxlambda) < 0) {
      opserr << "WARNING integrator MinUnbalDispNorm invalid maxlambda\n";
      return TCL_ERROR;
    }
  } else {
    minlambda = lambda11;
    maxlambda = lambda11;
    numIter = 1;
  }

  int signFirstStepMethod = MinUnbalDispNorm::SIGN_LAST_STEP;
  if (OPS_GetNumRemainingInputArgs() > 0) {
    const char *flag = OPS_GetString();
    if ((strcmp(flag, "-determinant") == 0) || (strcmp(flag, "-det") == 0)) {
      signFirstStepMethod = MinUnbalDispNorm::CHANGE_DETERMINANT;
    }
  }


  builder->set(*new MinUnbalDispNorm(lambda11, numIter, minlambda, maxlambda,
                                     signFirstStepMethod));
  return TCL_OK;
}



#include <analysis/integrator/Dynamic/Newmark1.h>
int
TclCommand_createNewmark1(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);

  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata != 2 && numdata != 6) {
    opserr << "WARNING integrator Newmark1 gamma beta <alphaM> <betaKcurrent> "
              "<betaKi> <betaKlastCommitted>\n";
    return TCL_ERROR;
  }

  double data[6] = {0, 0, 0, 0, 0, 0};
  if (OPS_GetDoubleInput(&numdata, data) < 0) {
    opserr << "WARNING integrator Newmark1 invalid double inputs\n";
    return TCL_ERROR;
  }

  double gamma = data[0];
  double beta = data[1];
  double alphaM = data[2], betaK = data[3], betaKi = data[4], betaKc = data[5];

  if (numdata == 2)
    builder->set(* new Newmark1(gamma, beta));
  else
    builder->set(*new Newmark1(gamma, beta, alphaM, betaK, betaKi, betaKc));
  return TCL_OK;
}




#include <analysis/integrator/Dynamic/NewmarkExplicit.h>
int
TclCommand_createNewmarkExplicit(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);

  if (argc != 3) {
    opserr
        << "WARNING - incorrect number of args want NewmarkExplicit $gamma\n";
    return TCL_ERROR;
  }

  double gamma;
  if (OPS_GetDouble(&argc, &gamma) != 0) {
    opserr << "WARNING - invalid args want NewmarkExplicit $gamma\n";
    return TCL_ERROR;
  }

  TransientIntegrator *theIntegrator = new NewmarkExplicit(gamma);

  builder->set(*theIntegrator);
  return TCL_OK;
}



#if 0
#include <analysis/integrator/Dynamic/NewmarkHSFixedNumIter.h>
int
TclCommand_createNewmarkHSFixedNumIter(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  TransientIntegrator *theIntegrator = 0;

  argc -= 2;

  if (argc != 2 && argc != 4) {
    opserr << "WARNING - incorrect number of args want NewmarkHSFixedNumIter "
              "$gamma $beta <-polyOrder $O>\n";
    return TCL_ERROR;
  }

  double dData[2];
  int polyOrder = 2;
  bool updDomFlag = true;
  int numData = 2;

  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING - invalid args want NewmarkHSFixedNumIter $gamma $beta "
              "<-polyOrder $O>\n";
    return TCL_ERROR;
  }

  if (argc == 4) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-polyOrder") == 0) {
      numData = 1;
      if (OPS_GetInt(&numData, &polyOrder) != 0) {
        opserr << "WARNING - invalid polyOrder want NewmarkHSFixedNumIter "
                  "$gamma $beta <-polyOrder $O>\n";
      }
    }
  }

  theIntegrator =
      new NewmarkHSFixedNumIter(dData[0], dData[1], polyOrder, updDomFlag);


  builder->set(*theIntegrator);
  return TCL_OK;
}




#include <analysis/integrator/Dynamic/NewmarkHSIncrLimit.h>
int
TclCommand_createNewmarkHSIncrLimit(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  TransientIntegrator *theIntegrator = 0;

  argc -= 2;

  if (argc != 3 && argc != 5) {
    opserr << "WARNING - incorrect number of args want NewmarkHSIncrLimit "
              "$gamma $beta $limit <-normType $T>\n";
    return TCL_ERROR;
  }

  double dData[3];
  int normType = 2;
  int numData = 3;

  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING - invalid args want NewmarkHSIncrLimit $gamma $beta "
              "$limit <-normType $T>\n";
    return TCL_ERROR;
  }

  if (argc == 5) {
    const char *argvLoc = OPS_GetString();
    if (strcmp(argvLoc, "-normType") == 0) {
      numData = 1;
      if (OPS_GetInt(&numData, &normType) != 0) {
        opserr << "WARNING - invalid normType want NewmarkHSIncrLimit $gamma "
                  "$beta $limit <-normType $T>\n";
      }
    }
  }

  theIntegrator =
      new NewmarkHSIncrLimit(dData[0], dData[1], dData[2], normType);


  builder->set(*theIntegrator);
  return TCL_OK;
}




#include <analysis/integrator/Dynamic/NewmarkHSIncrReduct.h>
int
TclCommand_createNewmarkHSIncrReduct(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  TransientIntegrator *theIntegrator = 0;

  argc -= 2;

  if (argc != 3) {
    opserr << "WARNING - incorrect number of args want NewmarkHSIncrReduct "
              "$gamma $beta $reduct\n";
    return TCL_ERROR;
  }

  double dData[3];
  if (OPS_GetDouble(&argc, dData) != 0) {
    opserr << "WARNING - invalid args want NewmarkHSIncrReduct $gamma $beta "
              "$reduct\n";
    return TCL_ERROR;
  }

  theIntegrator = new NewmarkHSIncrReduct(dData[0], dData[1], dData[2]);


  builder->set(*theIntegrator);
  return TCL_OK;
}
#endif



#include <analysis/integrator/Dynamic/Newmark.h>
int
TclCommand_createNewmark(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);

  TransientIntegrator *theIntegrator = nullptr;

  argc -= 2;

  if (argc != 2 && argc != 4) {
    opserr << "WARNING - incorrect number of args want Newmark $gamma $beta "
              "<-form $typeUnknown>\n";
    return TCL_ERROR;
  }

  int dispFlag = 1;
  double dData[2];
  int numData = 2;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING - invalid args want Newmark $gamma $beta <-form "
              "$typeUnknown>\n";
    return TCL_ERROR;
  }

  if (argc == 2)
    theIntegrator = new Newmark(dData[0], dData[1]);
  else {
    //    char nextString[10];
    const char *nextString = OPS_GetString();
    //    OPS_GetString(nextString, 10);
    if (strcmp(nextString, "-form") == 0) {
      //      OPS_GetString(nextString, 10);
      nextString = OPS_GetString();
      if ((nextString[0] == 'D') || (nextString[0] == 'd'))
        dispFlag = 1;
      else if ((nextString[0] == 'A') || (nextString[0] == 'a'))
        dispFlag = 3;
      else if ((nextString[0] == 'V') || (nextString[0] == 'v'))
        dispFlag = 2;
    }
    theIntegrator = new Newmark(dData[0], dData[1], dispFlag);
  }

  builder->set(*theIntegrator);
  return TCL_OK;
}



#if 0
#include <analysis/integrator/StagedLoadControl.h>
int
TclCommand_createStagedLoadControlIntegrator()
{
  if (OPS_GetNumRemainingInputArgs() < 1) {
    opserr << "insufficient arguments\n";
    return TCL_ERROR;
  }

  double lambda;
  int numData = 1;
  if (OPS_GetDoubleInput(&numData, &lambda) < 0) {
    opserr << "WARNING failed to read double lambda\n";
    return TCL_ERROR;
  }

  int numIter = 1;
  double mLambda[2] = {lambda, lambda};
  if (OPS_GetNumRemainingInputArgs() > 2) {
    if (OPS_GetIntInput(&numData, &numIter) < 0) {
      opserr << "WARNING failed to read int numIter\n";
      return TCL_ERROR;
    }
    numData = 2;
    if (OPS_GetDoubleInput(&numData, &mLambda[0]) < 0) {
      opserr << "WARNING failed to read double min and max\n";
      return TCL_ERROR;
    }
  }

  return new StagedLoadControl(lambda, numIter, mLambda[0], mLambda[1]);
}
#endif



#include <analysis/integrator/Dynamic/WilsonTheta.h>

int
TclCommand_createWilsonTheta(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
  BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);

  argc -= 2;

  if (argc != 1) {
    opserr << "WARNING - incorrect number of args want WilsonTheta $theta\n";
    return TCL_ERROR;
  }

  double theta;
  if (OPS_GetDouble(&argc, &theta) != 0) {
    opserr << "WARNING - invalid args want WilsonTheta $theta\n";
    return TCL_ERROR;
  }

  TransientIntegrator *theIntegrator = new WilsonTheta(theta);


  builder->set(*theIntegrator);
  return TCL_OK;
}
