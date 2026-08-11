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
// Description: This file implements commands that configure a 
// `ModelBuider`, including "model"
//
// Author: cmp
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _MSC_VER 
#  include <string.h>
#  define strcasecmp _stricmp
#else
#  include <strings.h>
#endif
#include <Logging.h>
#include <Parsing.h>
#include <runtimeAPI.h>
#include <Domain.h>
#include <FE_Datastore.h>

#include <ModelRegistry.h>
#include <modeling/commands.h>
#include <runtime/interpreter/Interpreter.h>


using namespace OpenSees;

bool builtModel = false;

FE_Datastore *theDatabase = nullptr;

extern int XaraInit_AnalysisCommands(Tcl_Interp *, ModelRegistry&);
extern int XaraInit_DomainCommands(Tcl_Interp *, Domain*);
extern int RemoveTclDomainCommands(Tcl_Interp* interp);

// 
int
XaraCmd_model(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char *argv[])
{
  G3_Runtime *rt = G3_getRuntime(interp);
  Domain *theNewDomain = (Domain*)clientData;

  ModelRegistry *theNewBuilder = nullptr;
  Rotations::Parameters rotationType = Rotations::Parameters::None;
  {
    const char* rotation_name = getenv("XARA_ROTATE");
    if (rotation_name != nullptr) {
      if (strcmp(rotation_name, "none") == 0) {
        rotationType = Rotations::Parameters::None;
      }
      else if (strcmp(rotation_name, "iter") == 0) {
        rotationType = Rotations::Parameters::Iter;
      }
      else if (strcmp(rotation_name, "incr") == 0) {
        rotationType = Rotations::Parameters::Incr;
      }
      else if (strcmp(rotation_name, "init") == 0) {
        rotationType = Rotations::Parameters::Init;
      }
      else {
        opserr << OpenSees::PromptValueError 
               << "invalid rotation type in environment variable XARA_ROTATE: '" 
               << rotation_name 
               << "'\n";
        return TCL_ERROR;
      }
    }
  }

  //
  //
  //
  bool isNewModel = (clientData == nullptr);
  if (clientData == nullptr) {
    theNewDomain = new Domain();

    // TODO: remove ops_TheActiveDomain
    ops_TheActiveDomain = theNewDomain;

    Tcl_CreateCommand(interp, "model", &XaraCmd_model, theNewDomain, nullptr);

    XaraInit_DomainCommands(interp, theNewDomain);
  }


  // make sure at least one other argument to contain model builder type given
  if (argc < 2) {
    opserr << OpenSees::PromptValueError << "need to specify a model type, valid types:\n";
    opserr << "\tBasicBuilder\n";
    return TCL_ERROR;
  }

  // check argv[1] for type of ModelBuilder and create the object
  if ((strcmp(argv[1], "basic") == 0)        ||
      (strcmp(argv[1], "Basic") == 0)        ||
      (strcmp(argv[1], "-ndm") == 0)         ||
      (strcmp(argv[1], "BasicBuilder") == 0) ||
      (strcmp(argv[1], "basicBuilder") == 0)) {

    if (argc < 3) {
      opserr << OpenSees::PromptValueError 
             << "incorrect number of arguments\n";
      return TCL_ERROR;
    }
    int ndm = 0;
    int ndf = 0;

    int posArg = 1; // track positional argument
    int argPos = 2;
    while (argPos < argc) {
      if (strcmp(argv[argPos], "-ndm") == 0 ||
          strcmp(argv[argPos], "-NDM") == 0) {
        argPos++;
        if (argPos < argc) {
          if (Tcl_GetInt(interp, argv[argPos], &ndm) != TCL_OK) {
            opserr << OpenSees::PromptValueError 
                   << "error reading ndm, got '" << argv[argPos] << "'\n";
            return TCL_ERROR;
          }
        }
        argPos++;
        posArg++;

      } else if (strcmp(argv[argPos], "-ndf") == 0 ||
                 strcmp(argv[argPos], "-NDF") == 0) {
        argPos++;
        if (argPos < argc)
          if (Tcl_GetInt(interp, argv[argPos], &ndf) != TCL_OK) {
            opserr << OpenSees::PromptValueError 
                   << "invalid parameter ndf"
                   << OpenSees::SignalMessageEnd;
            return TCL_ERROR;
          }
        argPos++;
        posArg++;
      }
      else if (strcmp(argv[argPos], "-rotation") == 0) {
        argPos++;
        if (argPos >= argc) {
          opserr << OpenSees::PromptValueError 
                 << "missing rotation type after -rotation"
                 << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        }
        if (strcasecmp(argv[argPos], "none") == 0) {
          rotationType = Rotations::Parameters::None;
        }
        else if (strcasecmp(argv[argPos], "iter") == 0) {
          rotationType = Rotations::Parameters::Iter;
        }
        else if (strcasecmp(argv[argPos], "incr") == 0) {
          rotationType = Rotations::Parameters::Incr;
        }
        else if (strcasecmp(argv[argPos], "init") == 0) {
          rotationType = Rotations::Parameters::Init;
        }
        else {
          opserr << OpenSees::PromptValueError << "invalid rotation type '" << argv[argPos] << "'\n";
          return TCL_ERROR;
        }
      } 
      else if (posArg == 1) {
        if (Tcl_GetInt(interp, argv[argPos], &ndm) != TCL_OK) {
          opserr << OpenSees::PromptValueError 
                 << "invalid parameter ndm"
                 << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        }
        argPos++;
        posArg++;
      }
      else if (posArg == 2) {
        if (Tcl_GetInt(interp, argv[argPos], &ndf) != TCL_OK) {
          opserr << OpenSees::PromptValueError 
                 << "error reading ndf: " 
                 << argv[argPos] 
                 << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        }
        argPos++;
        posArg++;

      } else {
        // no matches, advance to next argument
        argPos++;
      }
    }

    // check that ndm was specified
    if (ndm == 0) {
      opserr << OpenSees::PromptValueError
             << "missing required argument ndm\n";
      return TCL_ERROR;
    }

    // check for ndf, if not assume one
    if (ndf == 0) {
      if (ndm == 1)
        ndf = 1;
      else if (ndm == 2)
        ndf = 3;
      else if (ndm == 3)
        ndf = 6;
      else {
        opserr << OpenSees::PromptValueError << "specified ndm, " << ndm << ", will not work\n";
        opserr << "        with any elements in BasicBuilder\n";
        return TCL_ERROR;
      }
    }

    // TODO: remove this
    int G3_setDomain(G3_Runtime*, Domain*);
    G3_setDomain(rt, theNewDomain);
    // create the model builder
#if 1
    if (!isNewModel) {
      theNewBuilder = G3_getModelBuilder(rt); //static_cast<ModelRegistry*>(clientData);
      theNewBuilder->setDimension(ndm, ndf);
    }
    else
#endif
    theNewBuilder = new ModelRegistry(*theNewDomain, ndm, ndf, rotationType);

    //
    // Add model commands
    //
#ifdef MODEL_CHANNELS
    theNewBuilder->getParallelContext().setup(interp);
#endif 


    static int ncmd = sizeof(ModelBuilderCommands)/sizeof(decltype(ModelBuilderCommands[0])); // CommandTableEntry);

    Tcl_CreateCommand(interp, "wipe", XaraCmd_wipe, (ClientData)theNewBuilder, nullptr);

    for (int i = 0; i < ncmd; i++)
      Tcl_CreateCommand(interp, 
          ModelBuilderCommands[i].name,
          ModelBuilderCommands[i].func,
          (ClientData) theNewBuilder, nullptr);

    Tcl_SetAssocData(interp, "OPS::theTclBuilder", NULL, (ClientData)theNewBuilder);
    Tcl_SetAssocData(interp, "OPS::theBasicModelBuilder", NULL, (ClientData)theNewBuilder);
    Tcl_SetAssocData(interp, "OPS::theTclDomain", NULL, (ClientData)theNewDomain);

    G3_setModelBuilder(rt, theNewBuilder);

    const char* analysis_option;
    if (!(analysis_option = Tcl_GetVar(interp,"opensees::pragma::analysis",TCL_GLOBAL_ONLY)) ||
         (strcmp(analysis_option, "off") != 0)) {
      XaraInit_AnalysisCommands(interp, *theNewBuilder);
    }
  }
  else {
    opserr << OpenSees::PromptValueError 
           << "unknown model builder type '" << argv[1] << "' not supported"
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  return TCL_OK;
}


int
XaraCmd_wipe(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char *argv[])
{
  Tcl_Eval(interp, "_clearAnalysis");

  ModelRegistry *builder = static_cast<ModelRegistry*>(clientData);

  if (theDatabase != nullptr)
    delete theDatabase;

  if (builder != nullptr) {
    Domain* theDomain = builder->getDomain();
    theDomain->clearAll();
    ops_TheActiveDomain = nullptr;
    delete theDomain;
    delete builder;

    RemoveTclDomainCommands(interp);
    static int ncmd = sizeof(ModelBuilderCommands)/sizeof(decltype(ModelBuilderCommands[0]));
    for (int i = 0; i < ncmd; i++)
      Tcl_DeleteCommand(interp, ModelBuilderCommands[i].name);
    builtModel = false;
  }
  Tcl_CreateCommand(interp, "model", &XaraCmd_model, nullptr, nullptr);
  Tcl_CreateCommand(interp, "wipe",  &XaraCmd_wipe,    nullptr, nullptr);

  ops_Dt = 0.0;

  theDatabase = nullptr;

  // the domain deletes the record objects,
  // just have to delete the private array
  return TCL_OK;
}

// command invoked to invoke buildFE_Model() on the ModelBuilder
int
XaraCmd_build(ClientData context, Tcl_Interp *interp, ArgSize argc, TCL_Char *argv[])
{
  return TCL_OK;
}

