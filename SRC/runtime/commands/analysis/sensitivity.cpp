//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
#include <tcl.h>
#include <Logging.h>
#include <Parsing.h>
#include <Integrator.h>
#include <StaticIntegrator.h>
#include <Domain.h>
#include <LoadPattern.h>
#include <Parameter.h>
#include <ParameterIter.h>
#include <TransientIntegrator.h>
#include <BasicAnalysisBuilder.h>


int 
computeGradients(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
    BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);
    Integrator* theIntegrator = nullptr;

    if(builder->getStaticIntegrator() != nullptr) {
        theIntegrator = builder->getStaticIntegrator();

    } else if(builder->getTransientIntegrator() != nullptr) {
        theIntegrator = builder->getTransientIntegrator();
    }

    if (theIntegrator == nullptr) {
        opserr << OpenSees::PromptValueError << "No integrator is created\n";
        return TCL_ERROR;
    }

    if (theIntegrator->computeSensitivities() < 0) {
      opserr << OpenSees::PromptValueError << "failed to compute sensitivities\n";
      return TCL_ERROR;
    }
    
    return TCL_OK;
}


int
TclCommand_sensLambda(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  BasicAnalysisBuilder* builder = static_cast<BasicAnalysisBuilder*>(clientData);

  if (argc < 3) {
    opserr << OpenSees::PromptValueError << "insufficient arguments\n";
    return TCL_ERROR;
  }

  int pattern, paramTag;
  if (Tcl_GetInt(interp, argv[1], &pattern) != TCL_OK) {
    opserr << "ERROR reading load pattern tag\n";
    return TCL_ERROR;
  }

  LoadPattern *thePattern = builder->getDomain()->getLoadPattern(pattern);
  if (thePattern == nullptr) {
    opserr << "ERROR load pattern with tag " << pattern
           << " not found in domain\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2], &paramTag) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "sensLambda patternTag?  paramTag? - could not read "
              "paramTag? ";
    return TCL_ERROR;
  }

  Parameter *theParam = builder->getDomain()->getParameter(paramTag);
  if (theParam == nullptr) {
    opserr << OpenSees::PromptValueError 
           << "sensLambda: parameter " << paramTag << " not found" << "\n";
    return TCL_ERROR;
  }

  int gradIndex = theParam->getGradIndex();
  double value = thePattern->getLoadFactorSensitivity(gradIndex);

  Tcl_SetObjResult(interp, Tcl_NewDoubleObj(value));

  return TCL_OK;
}


int
TclCommand_sensitivityAlgorithm(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
    BasicAnalysisBuilder *builder = static_cast<BasicAnalysisBuilder*>(clientData);


    Integrator* theIntegrator = nullptr;

    if (builder->getStaticIntegrator() != nullptr) {
        theIntegrator = builder->getStaticIntegrator();

    } else if(builder->getTransientIntegrator() != nullptr) {
        theIntegrator = builder->getTransientIntegrator();
    }


    if (argc < 2) {
        opserr << "ERROR: Wrong number of parameters to sensitivity algorithm." << "\n";
        return TCL_ERROR;
    }

    if (theIntegrator == nullptr) {
        opserr << "The integrator needs to be instantiated before " << "\n"
               << " setting  sensitivity algorithm." << "\n";
        return -1;
    }


    // 1: compute at each step (default); 
    // 2: compute by command; 

    int analysisTypeTag = 1;
    if (strcmp(argv[1],"-computeAtEachStep") == 0)
        analysisTypeTag = 1;

    else if (strcmp(argv[1],"-computeByCommand") == 0)
        analysisTypeTag = 2;

    else {
        opserr << "Unknown sensitivity algorithm option: " << argv[1] << "\n";
        return TCL_ERROR;
    }

    theIntegrator->setComputeType(analysisTypeTag);
    theIntegrator->activateSensitivityKey();
    
    return TCL_OK;
}

