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
#include <StaticPattern.h>
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
// Forward declarations
//
const char *getInterpPWD(Tcl_Interp *interp);
extern "C" int OPS_ResetInputNoBuilder(ClientData,
                                       Tcl_Interp *, int cArg, int mArg,
                                       TCL_Char ** const argv, Domain *);

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
  {"rayleigh",            &TclCommand_rayleighDamping},
  
  {"getLoadFactor",       &XaraCmd_getLoadFactor},

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
  {"setNodeCoord",        &XaraCmd_setNodeCoord},
  {"setNodePressure",     &setNodePressure},

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
  // Element response
  {"localForce",          &OpenSees::DomainCommands::localForce         },
  {"eleType",             &OpenSees::DomainCommands::eleType,           },
  {"eleNodes",            &OpenSees::DomainCommands::eleNodes,          },
  {"getEleTags",          &OpenSees::DomainCommands::getEleTags,        },
  {"getNumElements",      &OpenSees::DomainCommands::getNumElements,    },
  {"getEleClassTags",     &OpenSees::DomainCommands::getEleClassTags,   },
  {"eleForce",            &OpenSees::DomainCommands::eleForce,          },
  {"eleResponse",         &OpenSees::DomainCommands::eleResponse,       },
  {"eleDynamicalForce",   &OpenSees::DomainCommands::eleDynamicalForce, },
  {"updateElementDomain", &OpenSees::DomainCommands::updateElementDomain},

  {"sectionForce",        &sectionForce},
  {"sectionTag",          &sectionTag},
  {"sectionDisplacement", &sectionDisplacement},
  {"sectionDeformation",  &sectionDeformation},
  {"sectionStiffness",    &sectionStiffness},
  {"sectionFlexibility",  &sectionFlexibility},
  {"sectionLocation",     &sectionLocation},
  {"sectionWeight",       &sectionWeight},

  // Recorders
  {"recorderValue",       &OPS_recorderValue},
  {"record",              &TclCommand_record},

  {"InitialStateAnalysis", &InitialStateAnalysis},
};
}


static int 
TclCommand_BadDomainCommand(ClientData, Tcl_Interp* interp, 
                            Tcl_Size argc, 
                            TCL_Char** const argv)
{
  opserr << OpenSees::PromptModelError
         << "The model does not exist or has been destroyed."
         << OpenSees::SignalMessageEnd;
  return TCL_ERROR;
}


int
RemoveTclDomainCommands(Tcl_Interp* interp)
{
  for (size_t i = 0; i < sizeof(domainCommands) / sizeof(domainCommands[0]); ++i) {
    // Tcl_DeleteCommand(interp, domainCommands[i].name);
    Tcl_CreateCommand(interp, domainCommands[i].name,  &TclCommand_BadDomainCommand, nullptr, nullptr);
  }
  return 0;
}

int
AddTclDomainCommands(Tcl_Interp *interp, Domain* the_domain)
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
    // damping
    Tcl_CreateCommand(interp, "setElementRayleighDampingFactors", &addElementRayleigh, domain, nullptr);
    Tcl_CreateCommand(interp, "setElementRayleighFactors",        &addElementRayleigh, domain, nullptr);
    // Modal
    Tcl_CreateCommand(interp, "modalProperties",     &modalProperties,     domain, nullptr);
  }

  for (size_t i = 0; i < sizeof(domainCommands) / sizeof(domainCommands[0]); ++i) {
    Tcl_CreateCommand(interp, domainCommands[i].name,
                      domainCommands[i].func, 
                      domain, nullptr);
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


#include <StaticPattern.h>

int
XaraCmd_getLoadFactor(ClientData clientData, 
                      Tcl_Interp *interp, 
                      ArgSize argc,
                      TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain* domain = (Domain*)clientData; 

  if (argc < 2) {
    opserr << OpenSees::PromptValueError 
           << "no load pattern supplied"
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  Xara::Tag tag;
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "reading load pattern tag\n";
    return TCL_ERROR;
  }

  LoadPattern *the_pattern = domain->getLoadPattern(tag);
  if (the_pattern == nullptr) {
    opserr << OpenSees::PromptValueError 
           << "load pattern with tag " << tag
           << " not found in domain\n";
    return TCL_ERROR;
  }

  if (the_pattern->getClassTag() != LoadPattern::PATTERN_TAG_StaticPattern) {
    opserr << OpenSees::PromptValueError 
           << "load pattern with tag " << tag
           << " is not a StaticPattern\n";
    return TCL_ERROR;
  }

  StaticPattern* theStaticPattern = static_cast<StaticPattern*>(the_pattern);

  double factor = theStaticPattern->getLoadFactor();
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
    opslog << "InitialStateAnalysis ON" << "\n";

    // set global variable to true
    // FMK changes for parallel:
    // ops_InitialStateAnalysis = true;

    Parameter *theP = new InitialStateParameter(true);
    the_domain->addParameter(theP);
    delete theP;

    return TCL_OK;
  }
  else if (strcmp(argv[1], "off") == 0) {
    opslog << "InitialStateAnalysis OFF" << "\n";

    // call revert to start to zero the displacements
    the_domain->revertToStart();

    // set global variable to false
    // FMK changes for parallel
    // ops_InitialStateAnalysis = false;
    Parameter *theP = new InitialStateParameter(false);
    the_domain->addParameter(theP);
    delete theP;

    return TCL_OK;
  }
  else {
    opserr << "WARNING: Incorrect arguments - want InitialStateAnalysis on, or "
              "InitialStateAnalysis off"
           << "\n";

    return TCL_ERROR;
  }
}

int
TclCommand_rayleighDamping(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc,
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

  Tcl_Obj *result = Tcl_NewListObj(0, nullptr);

  if (argc == 1) {
    LoadPattern *thePattern;
    LoadPatternIter &thePatterns = the_domain->getLoadPatterns();

    while ((thePattern = thePatterns()) != nullptr) {
      if (thePattern->getClassTag() != LoadPattern::PATTERN_TAG_StaticPattern)
        continue;
      StaticPattern *theStaticPattern = static_cast<StaticPattern *>(thePattern);
      ElementalLoadIter theEleLoads = theStaticPattern->getElementalLoads();
      ElementalLoad *theLoad;

      while ((theLoad = theEleLoads()) != nullptr) {
        Tcl_ListObjAppendElement(interp, result, Tcl_NewIntObj(theLoad->getClassTag()));
      }
    }
  }
  else {
    // Either of:
    // 1) getEleLoadClassTags <tag>
    // 2) getEleLoadClassTags -pattern <tag>
    int patternTag;
    int arg_pattern = 1;
    if (argc == 3) {
      if (strcmp(argv[1], "-pattern") != 0) {
        opserr << OpenSees::PromptValueError 
               << "Unexpected argument " 
               << argv[1] 
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      arg_pattern = 2;
    }

    if (Tcl_GetInt(interp, argv[arg_pattern], &patternTag) != TCL_OK) {
      opserr << OpenSees::PromptValueError 
             << "Failed to read patternTag"
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }

    LoadPattern *thePattern = the_domain->getLoadPattern(patternTag);
    if (thePattern == nullptr) {
      opserr << OpenSees::PromptValueError 
             << "Failed to find load pattern with tag " 
             << patternTag
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }
    if (thePattern->getClassTag() != LoadPattern::PATTERN_TAG_StaticPattern) {
      opserr << OpenSees::PromptValueError 
             << "Load pattern with tag " << patternTag
             << " is not a StaticPattern"
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }

    StaticPattern *theStaticPattern = static_cast<StaticPattern *>(thePattern);

    ElementalLoadIter theEleLoads = theStaticPattern->getElementalLoads();

    ElementalLoad *theLoad;
    while ((theLoad = theEleLoads()) != nullptr) {
      Tcl_ListObjAppendElement(interp, result, Tcl_NewIntObj(theLoad->getClassTag()));
    }
  }

  Tcl_SetObjResult(interp, result);
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

  Tcl_Obj *result = Tcl_NewListObj(0, nullptr);

  if (argc == 1) {
    LoadPattern *thePattern;
    LoadPatternIter &thePatterns = the_domain->getLoadPatterns();

    while ((thePattern = thePatterns()) != nullptr) {
      if (thePattern->getClassTag() != LoadPattern::PATTERN_TAG_StaticPattern)
        continue;

      StaticPattern *theStaticPattern = static_cast<StaticPattern *>(thePattern);
      ElementalLoadIter theEleLoads = theStaticPattern->getElementalLoads();
      ElementalLoad *theLoad;

      while ((theLoad = theEleLoads()) != nullptr) {
        Tcl_ListObjAppendElement(interp, result, Tcl_NewIntObj(theLoad->getElementTag()));
      }
    }
  }
  else if (argc == 2) {
    int patternTag;

    if (Tcl_GetInt(interp, argv[1], &patternTag) != TCL_OK) {
      opserr << OpenSees::PromptValueError 
             << "failed to read patternTag"
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }

    LoadPattern *thePattern = the_domain->getLoadPattern(patternTag);
    if (thePattern == nullptr) {
      opserr << OpenSees::PromptValueError 
             << "Load pattern with tag " << patternTag
             << " not found in domain\n";
      return TCL_ERROR;
    }
    if (thePattern->getClassTag() != LoadPattern::PATTERN_TAG_StaticPattern) {
      opserr << OpenSees::PromptValueError 
             << "Load pattern with tag " << patternTag
             << " is not a StaticPattern\n";
      return TCL_ERROR;
    }

    StaticPattern *theStaticPattern = static_cast<StaticPattern *>(thePattern);

    ElementalLoadIter theEleLoads = theStaticPattern->getElementalLoads();
    ElementalLoad *theLoad;

    while ((theLoad = theEleLoads()) != nullptr) {
      Tcl_ListObjAppendElement(interp, result, Tcl_NewIntObj(theLoad->getElementTag()));
    }
  }
  else {
    opserr << OpenSees::PromptValueError << "unexpected arguments\n" << "\n";
    return TCL_ERROR;
  }

  Tcl_SetObjResult(interp, result);

  return TCL_OK;
}

int
getEleLoadData(ClientData clientData, 
               Tcl_Interp *interp, Tcl_Size argc,
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
      if (thePattern->getClassTag() != LoadPattern::PATTERN_TAG_StaticPattern)
        continue;
      StaticPattern *theStaticPattern = static_cast<StaticPattern *>(thePattern);
      ElementalLoadIter &theEleLoads = theStaticPattern->getElementalLoads();
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
      opserr << OpenSees::PromptValueError 
             << "Load pattern with tag " << patternTag
             << " not found in domain\n";
      return TCL_ERROR;
    }
    if (thePattern->getClassTag() != LoadPattern::PATTERN_TAG_StaticPattern) {
      opserr << OpenSees::PromptValueError 
             << "Load pattern with tag " << patternTag
             << " is not a StaticPattern\n";
      return TCL_ERROR;
    }

    StaticPattern *theStaticPattern = static_cast<StaticPattern *>(thePattern);
    ElementalLoadIter theEleLoads = theStaticPattern->getElementalLoads();
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

