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
// Written: cmp
// Created: Spring 2023
//
#include <string>
#include <unordered_map>

#include <Parsing.h>
#include <Logging.h>

Tcl_CmdProc TclCommand_useMaterial;
Tcl_CmdProc TclCommand_useUniaxialMaterial;
Tcl_CmdProc TclCommand_useCrossSection;
Tcl_CmdProc TclCommand_usePlaneStress;

const std::unordered_map<std::string, Tcl_CmdProc*> invoke_commands 
{
  {"UniaxialMaterial",    &TclCommand_useUniaxialMaterial       },

  {"FrameSection",        &TclCommand_useCrossSection           },
  {"section",             &TclCommand_useCrossSection           },

  {"PlaneStress",         &TclCommand_usePlaneStress            },
  {"Material",            &TclCommand_useMaterial               },
  {"TriaxialMaterial",    &TclCommand_useMaterial               },
  {"MultiaxialMaterial",  &TclCommand_useMaterial               },
};


int
XaraCmd_invoke(ClientData context, Tcl_Interp* interp, ArgSize argc, char const** const argv)
{
  // check number of arguments in command line
  if (argc < 4) {
    opserr << OpenSees::PromptValueError 
           << "bad arguments - want: using <obj-type> <obj-tag> {<operations>...}"
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  auto tcl_cmd = invoke_commands.find(std::string(argv[1]));

  if (tcl_cmd != invoke_commands.end()) {

    return (*tcl_cmd->second)(context, interp, argc, &argv[0]);

  } else {
    return TCL_ERROR;
  }

  return TCL_OK;

}

