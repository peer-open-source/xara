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
// This file contains functions that are required by Tcl to load the
// Xara library.
//
// The following macros are defined by the build system:
// - XARA_VERSION: The version of the Xara library.
// - XARA_COMMIT_HASH (optional): The git commit hash of the current build.
// - XARA_PARALLEL_MODE: String, either MP or RT
//
#ifndef XARA_VERSION
#  define XARA_VERSION "0.0.0"
#endif
//
#include <Parsing.h>
#include <runtimeAPI.h>
#include "runtime/G3_Runtime.h"
#include <logging/Logging.h>
#include <handler/OPS_Stream.h>
#include <StandardStream.h>
#include "commands/strings.cpp"
#include <stdio.h>
#include <stdlib.h>

// Determine when stdout is a TTY
#ifdef _WIN32
#  include <io.h>
#  define isatty _isatty
#  define STDERR_FILENO _fileno(stderr)
#else
#  include <unistd.h>               
#endif

// interpreter/runtime.cpp
extern int  XaraInit_InterpreterCommands(Tcl_Interp *interp);
extern void XaraInit_SequentialCommands(Tcl_Interp* interp);
extern int  XaraInit_UtilityCommands(Tcl_Interp*);

//
// Tcl Command that returns the current OpenSees version
//
static int
XaraCmd_version(ClientData clientData, Tcl_Interp *interp, ArgSize argc, TCL_Char **argv)
{
  char buffer[40];

#ifdef XARA_COMMIT_HASH
  sprintf(buffer, "%s (commit %s)", XARA_VERSION, XARA_COMMIT_HASH);
#else
  sprintf(buffer, "%s", XARA_VERSION);
#endif
  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}


//
// Called when the library is loaded as a Tcl extension.
//
extern "C" int 
#ifdef _WIN32
__declspec(dllexport)
#endif
Openseesrt_Init(Tcl_Interp *interp)
{
  if (Tcl_InitStubs(interp, TCL_VERSION, 0) == NULL)
    return TCL_ERROR;

  if (Tcl_PkgProvide(interp, "OpenSeesRT", XARA_VERSION) == TCL_ERROR)
    return TCL_ERROR;

  // Create a runtime instance, and store it with the interpreter
  G3_Runtime *rt = new G3_Runtime{interp};
  Tcl_SetAssocData(interp, "G3_Runtime", NULL, (ClientData)rt);

  // Initialize OpenSees
  XaraInit_InterpreterCommands(interp);
  XaraInit_SequentialCommands(interp); // Add sequential API
  XaraInit_UtilityCommands(interp);    // Add utility commands (linspace, range, etc.)

  char* verbosity = getenv("XARA_VERBOSITY"); // Was OPENSEESRT_VERBOSITY
  if (verbosity != nullptr) {
    if (strcmp(verbosity, "DEBUG") == 0) {
      G3_SetStreamLevel(G3_LevelDebug, true);
    }
  }

  // Prevent coloring output when stderr is not a TTY
  if (isatty(STDERR_FILENO))
    G3_SetStreamColor(nullptr, G3_LevelWarn, 1);


  // Set variables with package information
  Tcl_SetVar(interp, "opensees::copyright", copyright,      TCL_LEAVE_ERR_MSG);
  Tcl_SetVar(interp, "opensees::license",   license,        TCL_LEAVE_ERR_MSG);
  Tcl_SetVar(interp, "opensees::banner",    unicode_banner, TCL_LEAVE_ERR_MSG);
  Tcl_CreateCommand(interp, "version",      XaraCmd_version,  nullptr, nullptr);
  return TCL_OK;
}

