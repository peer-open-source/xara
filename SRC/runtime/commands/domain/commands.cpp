//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
//
// Description: This file contains the functions that will be called by
// the interpreter when the appropriate command name is specified.
//
#include <Logging.h>
#include <Parsing.h>
#include <elementAPI.h>
#include <classTags.h>
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tgmath.h>
#include <assert.h>
#include <set>
#include <vector>
#include <algorithm>
//
#include "commands.h"
// Domain
#include <Domain.h>
#include <DOF_Group.h>
#include <Matrix.h>
#include <Node.h>
#include <Element.h>
#include <ElementIter.h>
#include <LoadPattern.h>
#include <ElementalLoad.h>
#include <ElementalLoadIter.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <MP_Constraint.h>
#include <MP_ConstraintIter.h>
//
#include <Parameter.h>
#include <ParameterIter.h>
#include <InitialStateParameter.h>
#include <ElementStateParameter.h>
#include <Pressure_Constraint.h>
// Analysis
#include <AnalysisModel.h>
#include <EquiSolnAlgo.h>
#include <StaticIntegrator.h>
#include <LinearSOE.h>
#include <EigenSOE.h>
// Other
#include <Recorder.h>
#include <DummyStream.h>
#include <Information.h>
#include <Response.h>
#include <packages.h>
//
// Global variables
//
class ModelBuilder;
ModelBuilder          *theBuilder         = nullptr;
//
// Forward declarations
//
const char *getInterpPWD(Tcl_Interp *interp);
extern "C" int OPS_ResetInputNoBuilder(ClientData clientData,
                                       Tcl_Interp *interp, int cArg, int mArg,
                                       TCL_Char ** const argv, Domain *domain);

Tcl_CmdProc TclCommand_record;
Tcl_CmdProc TclCommand_setLoadConst;
Tcl_CmdProc TclCommand_setCreep;

namespace {
static struct {
  const char *name;
  int (*func)(ClientData, Tcl_Interp *, Tcl_Size, TCL_Char ** const);
} domainCommands[] = {
  {"loadConst",           &TclCommand_setLoadConst},
  {"recorder",            &TclAddRecorder},
  {"region",              &TclCommand_addMeshRegion},

  {"printGID",            &printModelGID},

  {"setTime",             &TclCommand_setTime},
  {"getTime",             &TclCommand_getTime},
  {"setCreep",            &TclCommand_setCreep},

  // DAMPING
  {"rayleigh",            &rayleighDamping},
  
  {"getLoadFactor",       &getLoadFactor},

  //
  {"basicDeformation",    &basicDeformation},
  {"basicForce",          &basicForce},
  {"basicStiffness",      &basicStiffness},


  {"nodeDOFs",            &nodeDOFs},
  {"nodeCoord",           &nodeCoord},
  {"nodeMass",            &nodeMass},
  {"nodeVel",             &nodeVel},
  {"nodeDisp",            &nodeDisp},
  {"nodeAccel",           &nodeAccel},
  {"nodeResponse",        &nodeResponse},
  {"nodePressure",        &nodePressure},
  {"nodeBounds",          &nodeBounds},
  {"findNodeWithID",      &findID},
  {"nodeUnbalance",       &nodeUnbalance},
  {"nodeEigenvector",     &nodeEigenvector},
  {"nodeReaction",        &nodeReaction},

  {"reactions",           &calculateNodalReactions},

  {"setNodeVel",          &setNodeVel},
  {"setNodeDisp",         &setNodeDisp},
  {"setNodeAccel",        &setNodeAccel},
  {"setNodeCoord",        &setNodeCoord},
  {"nodeRotation",        &nodeRotation},
  {"getNodeTags",         &getNodeTags},



  {"getParamTags",        &getParamTags},
  {"getParamValue",       &getParamValue},
  {"parameter",           &TclCommand_parameter},
  {"addToParameter",      &TclCommand_parameter},
  {"updateParameter",     &TclCommand_parameter},
  {"setParameter",        &TclCommand_setParameter},


  {"getEleLoadTags",      &getEleLoadTags},
  {"getEleLoadData",      &getEleLoadData},
  {"getEleLoadClassTags", &getEleLoadClassTags},


  {"sectionForce",        &sectionForce},
  {"sectionTag",          &sectionTag},
  {"sectionDisplacement", &sectionDisplacement},
  {"sectionDeformation",  &sectionDeformation},
  {"sectionStiffness",    &sectionStiffness},
  {"sectionFlexibility",  &sectionFlexibility},
  {"sectionLocation",     &sectionLocation},
  {"sectionWeight",       &sectionWeight},

  {"recorderValue",       &OPS_recorderValue},
  {"record",              &TclCommand_record},

  {"InitialStateAnalysis", &InitialStateAnalysis},
};
}

// TODO: reimplement defaultUnits and setParameter
// int defaultUnits(ClientData, Tcl_Interp *, int, TCL_Char ** const argv);
// int setParameter(ClientData, Tcl_Interp *, int, TCL_Char **);
int
G3_AddTclDomainCommands(Tcl_Interp *interp, Domain* the_domain)
{

  ClientData domain = (ClientData)the_domain;


  {
    using namespace OpenSees::DomainCommands;
    // Domain
    Tcl_CreateObjCommand(interp, "fixedNodes",          &fixedNodes,          domain, nullptr);
    Tcl_CreateObjCommand(interp, "fixedDOFs",           &fixedDOFs,           domain, nullptr);
    Tcl_CreateObjCommand(interp, "constrainedNodes",    &constrainedNodes,    domain, nullptr);
    Tcl_CreateObjCommand(interp, "constrainedDOFs",     &constrainedDOFs,     domain, nullptr);
    Tcl_CreateObjCommand(interp, "domainChange",        &domainChange,        domain, nullptr);
    Tcl_CreateObjCommand(interp, "remove",              &removeObject,        domain, nullptr);
    Tcl_CreateCommand(interp,    "retainedNodes",       &retainedNodes,       domain, nullptr);
    Tcl_CreateCommand(interp,    "retainedDOFs",        &retainedDOFs,        domain, nullptr);
    // Elements
    Tcl_CreateCommand(interp, "localForce",          &localForce,    domain, nullptr);
    Tcl_CreateCommand(interp, "eleType",             &eleType,       domain, nullptr);
    Tcl_CreateCommand(interp, "eleNodes",            &eleNodes,            domain, nullptr);
    Tcl_CreateCommand(interp, "getEleTags",          &getEleTags,          domain, nullptr);
    Tcl_CreateCommand(interp, "getNumElements",      &getNumElements,      domain, nullptr);
    Tcl_CreateCommand(interp, "getEleClassTags",     &getEleClassTags,     domain, nullptr);
    Tcl_CreateCommand(interp, "eleForce",            &eleForce,            domain, nullptr);
    Tcl_CreateCommand(interp, "eleResponse",         &eleResponse,         domain, nullptr);
    Tcl_CreateCommand(interp, "eleDynamicalForce",   &eleDynamicalForce,   domain, nullptr);
    Tcl_CreateCommand(interp, "updateElementDomain", &updateElementDomain, nullptr, nullptr);
    // damping
    Tcl_CreateCommand(interp, "setElementRayleighDampingFactors", &addElementRayleigh, domain, nullptr);
    Tcl_CreateCommand(interp, "setElementRayleighFactors",        &addElementRayleigh, domain, nullptr);
    // Modal
    Tcl_CreateCommand(interp, "modalProperties",     &modalProperties,     domain, nullptr);
  }

  for (int i = 0; i < sizeof(domainCommands) / sizeof(domainCommands[0]); ++i) {
    Tcl_CreateCommand(interp, domainCommands[i].name,
                      domainCommands[i].func, domain, nullptr);
  }

  // sensitivity
  Tcl_CreateCommand(interp, "computeGradients",      &computeGradients, (ClientData)domain, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "sensitivityAlgorithm",  &TclCommand_sensitivityAlgorithm, (ClientData)domain, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "sensNodeDisp",          &sensNodeDisp, (ClientData)domain, (Tcl_CmdDeleteProc *)NULL);
//Tcl_CreateCommand(interp, "sensLambda",            &sensLambda, (ClientData)domain, (Tcl_CmdDeleteProc *)NULL); // Abbas
  Tcl_CreateCommand(interp, "sensNodeVel",           &sensNodeVel, (ClientData)domain, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "sensNodeAccel",         &sensNodeAccel, (ClientData)domain, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "sensSectionForce",      &sensSectionForce, (ClientData)domain, (Tcl_CmdDeleteProc *)NULL);
  Tcl_CreateCommand(interp, "sensNodePressure",      &sensNodePressure, (ClientData)domain, (Tcl_CmdDeleteProc *)NULL);


//   TODO: cmp, moved definition to packages/optimization; need to link in optionally

  // Tcl_CreateCommand(interp, "sdfResponse",      &sdfResponse, nullptr, nullptr);
  // Tcl_CreateCommand(interp, "database", &addDatabase, nullptr, nullptr);

  // wipeAnalysis(0, interp, 0, 0);
  return TCL_OK;
}



int
getLoadFactor(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc,
              TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain* domain = (Domain*)clientData; 

  if (argc < 2) {
    opserr << OpenSees::PromptValueError 
           << "no load pattern supplied -- getLoadFactor\n";
    return TCL_ERROR;
  }

  int pattern;
  if (Tcl_GetInt(interp, argv[1], &pattern) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "reading load pattern tag -- getLoadFactor\n";
    return TCL_ERROR;
  }

  LoadPattern *the_pattern = domain->getLoadPattern(pattern);
  if (the_pattern == nullptr) {
    opserr << OpenSees::PromptValueError << "load pattern with tag " << pattern
           << " not found in domain -- getLoadFactor\n";
    return TCL_ERROR;
  }

  double factor = the_pattern->getLoadFactor();
  Tcl_SetObjResult(interp, Tcl_NewDoubleObj(factor));

  return TCL_OK;
}



// added by C.McGann, U.Washington
int
InitialStateAnalysis(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc,
                     TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc < 2) {
    opserr << "WARNING: Incorrect number of arguments for InitialStateAnalysis "
              "command"
           << "\n";
    return TCL_ERROR;
  }

  if (strcmp(argv[1], "on") == 0) {
    opserr << "InitialStateAnalysis ON" << "\n";

    // set global variable to true
    // FMK changes for parallel:
    // ops_InitialStateAnalysis = true;

    Parameter *theP = new InitialStateParameter(true);
    the_domain->addParameter(theP);
    delete theP;

    return TCL_OK;

  } else if (strcmp(argv[1], "off") == 0) {
    opserr << "InitialStateAnalysis OFF" << "\n";

    // call revert to start to zero the displacements
    the_domain->revertToStart();

    // set global variable to false
    // FMK changes for parallel
    // ops_InitialStateAnalysis = false;
    Parameter *theP = new InitialStateParameter(false);
    the_domain->addParameter(theP);
    delete theP;

    return TCL_OK;

  } else {
    opserr << "WARNING: Incorrect arguments - want InitialStateAnalysis on, or "
              "InitialStateAnalysis off"
           << "\n";

    return TCL_ERROR;
  }
}

int
rayleighDamping(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc,
                TCL_Char ** const argv)
{
  //
  // rayleigh alphaM? betaK? betaK0? betaKc?
  //
  if (argc < 3) {
    opserr << OpenSees::PromptValueError
           << "not enough arguments to command\n";
    return TCL_ERROR;
  }

  double alphaM, betaK, betaK0=0.0, betaKc=0.0;
  if (Tcl_GetDouble(interp, argv[1], &alphaM) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "could not read alphaM? \n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[2], &betaK) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "could not read betaK? \n";
    return TCL_ERROR;
  }

  if (argc > 3 && Tcl_GetDouble(interp, argv[3], &betaK0) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "could not read betaK0? \n";
    return TCL_ERROR;
  }

  if (argc > 4 && Tcl_GetDouble(interp, argv[4], &betaKc) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "could not read betaKc? \n";
    return TCL_ERROR;
  }

  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;
  the_domain->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);
  return TCL_OK;
}


int
getEleLoadClassTags(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc,
                    TCL_Char ** const argv)
{
  //
  // getEleLoadClassTags <patternTag?>
  //
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc == 1) {
    LoadPattern *thePattern;
    LoadPatternIter &thePatterns = the_domain->getLoadPatterns();

    char buffer[20];

    while ((thePattern = thePatterns()) != nullptr) {
      ElementalLoadIter theEleLoads = thePattern->getElementalLoads();
      ElementalLoad *theLoad;

      while ((theLoad = theEleLoads()) != nullptr) {
        sprintf(buffer, "%d ", theLoad->getClassTag());
        Tcl_AppendResult(interp, buffer, NULL);
      }
    }
  }
  else if (argc == 2) {
    int patternTag;

    if (Tcl_GetInt(interp, argv[1], &patternTag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "failed to read patternTag\n";
      return TCL_ERROR;
    }

    LoadPattern *thePattern = the_domain->getLoadPattern(patternTag);
    if (thePattern == nullptr) {
      opserr << OpenSees::PromptValueError 
             << "load pattern with tag " << patternTag
             << " not found in domain"
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }

    ElementalLoadIter theEleLoads = thePattern->getElementalLoads();

    char buffer[20];

    ElementalLoad *theLoad;
    while ((theLoad = theEleLoads()) != nullptr) {
      sprintf(buffer, "%d ", theLoad->getClassTag());
      Tcl_AppendResult(interp, buffer, NULL);
    }

  } else {
    opserr << OpenSees::PromptValueError << "unexpected arguments\n" 
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  return TCL_OK;
}

int
getEleLoadTags(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc,
               TCL_Char ** const argv)
{
  //
  // getEleLoadTags <patternTag?>
  //
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc == 1) {
    LoadPattern *thePattern;
    LoadPatternIter &thePatterns = the_domain->getLoadPatterns();

    Tcl_Obj *result = Tcl_NewListObj(0, nullptr);

    while ((thePattern = thePatterns()) != nullptr) {
      ElementalLoadIter theEleLoads = thePattern->getElementalLoads();
      ElementalLoad *theLoad;

      while ((theLoad = theEleLoads()) != nullptr) {
        Tcl_ListObjAppendElement(interp, result, Tcl_NewIntObj(theLoad->getElementTag()));
      }
    }

    Tcl_SetObjResult(interp, result);

  } else if (argc == 2) {
    int patternTag;

    if (Tcl_GetInt(interp, argv[1], &patternTag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "failed to read patternTag \n";
      return TCL_ERROR;
    }

    LoadPattern *thePattern = the_domain->getLoadPattern(patternTag);
    if (thePattern == nullptr) {
      opserr << OpenSees::PromptValueError << "load pattern with tag " << patternTag
             << " not found in domain\n";
      return TCL_ERROR;
    }

    ElementalLoadIter theEleLoads = thePattern->getElementalLoads();
    ElementalLoad *theLoad;

    char buffer[20];

    while ((theLoad = theEleLoads()) != nullptr) {
      sprintf(buffer, "%d ", theLoad->getElementTag());
      Tcl_AppendResult(interp, buffer, NULL);
    }

  } else {
    opserr << OpenSees::PromptValueError << "unexpectd arguments\n" << "\n";
    return TCL_ERROR;
  }

  return TCL_OK;
}

int
getEleLoadData(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc,
               TCL_Char ** const argv)
{
  // getLoadData <patternTag?>
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc == 1) {
    LoadPattern *thePattern;
    LoadPatternIter &thePatterns = the_domain->getLoadPatterns();

    char buffer[40];
    int typeEL;

    while ((thePattern = thePatterns()) != nullptr) {
      ElementalLoadIter &theEleLoads = thePattern->getElementalLoads();
      ElementalLoad *theLoad;

      while ((theLoad = theEleLoads()) != nullptr) {
        const Vector &eleLoadData = theLoad->getData(typeEL, 1.0);

        int eleLoadDataSize = eleLoadData.Size();
        for (int i = 0; i < eleLoadDataSize; ++i) {
          sprintf(buffer, "%35.20f ", eleLoadData(i));
          Tcl_AppendResult(interp, buffer, NULL);
        }
      }
    }
  } 
  
  else if (argc == 2) {
    int patternTag;

    if (Tcl_GetInt(interp, argv[1], &patternTag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "failed to read patternTag\n";
      return TCL_ERROR;
    }

    LoadPattern *thePattern = the_domain->getLoadPattern(patternTag);
    if (thePattern == nullptr) {
      opserr << OpenSees::PromptValueError << "load pattern with tag " << patternTag
             << " not found in domain\n";
      return TCL_ERROR;
    }

    ElementalLoadIter theEleLoads = thePattern->getElementalLoads();
    ElementalLoad *theLoad;

    int typeEL;
    char buffer[40];

    while ((theLoad = theEleLoads()) != nullptr) {
      const Vector &eleLoadData = theLoad->getData(typeEL, 1.0);

      int eleLoadDataSize = eleLoadData.Size();
      for (int i = 0; i < eleLoadDataSize; ++i) {
        sprintf(buffer, "%35.20f ", eleLoadData(i));
        Tcl_AppendResult(interp, buffer, NULL);
      }
    }

  } else {
    opserr << OpenSees::PromptValueError 
           << "want - getEleLoadTags <patternTag?>" << "\n";
    return TCL_ERROR;
  }

  return TCL_OK;
}
