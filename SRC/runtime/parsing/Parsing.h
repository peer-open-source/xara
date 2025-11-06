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
#pragma once

#include <string.h>
#ifndef TCL_Char
typedef const char TCL_Char;
#endif
#ifndef G3_Char
typedef const char G3_Char;
#endif

// As suggested by https://core.tcl-lang.org/tcl/wiki?name=Migrating+C+extensions+to+Tcl+9
#ifndef TCL_SIZE_MAX
typedef int Tcl_Size;
#define TCL_SIZE_MAX INT_MAX
#endif

// #define OPS_API
#if defined(OPS_API)
# if !defined(_TCL)
  struct Tcl_Interp {};
  typedef const char Tcl_Obj;
  typedef void* ClientData;
  enum {
    TCL_OK    = 0,
    TCL_ERROR = 1
  };
  extern "C" {
  typedef int (Tcl_CmdProc)(ClientData, Tcl_Interp *, Tcl_Size, const TCL_Char *const []);
  int Tcl_GetDouble(Tcl_Interp *interp, const char *arg, double *value);
  int Tcl_GetInt(Tcl_Interp *interp, const char *arg, int *value);
  int Tcl_Free(char *ptr);
  int Tcl_SplitList(Tcl_Interp *interp, const char *list, int *argcPtr, TCL_Char ***argvPtr);
  char* Tcl_GetString(const char *arg);


  void* Tcl_GetAssocData(Tcl_Interp *interp, const char *name, void *deleteProc);
  int Tcl_SetAssocData(Tcl_Interp *interp, const char *name, void *deleteProc, void *clientData);
  int Tcl_Eval(Tcl_Interp *interp, const char *script);
  }
# endif

  namespace OpenSees {
  template<Tcl_CmdProc *Proc>
  int InvokeTclCommand(void* context) {
    Tcl_Size argc = 0;
    TCL_Char **argv = nullptr;
    Tcl_Interp interp {};
    return Proc(context, &interp, argc, argv);
  }
  } // namespace OpenSees
#else
# include "InputAPI.h"
#endif


class Parameter;
class Domain;
namespace OpenSees {
namespace Parsing {
int
GetDoubleParam(Tcl_Interp*, Domain&, const char* arg, double* value, Parameter*&);
}
}
