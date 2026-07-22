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
#pragma once
// #include <tcl.h>
// typedef const char** const ArgList;

// namespace Xara {
// typedef Tcl_Interp Interpreter;
// }
#include <tcl.h>
#include <StringStream.h>

namespace Xara {


#if 0
class Interpreter {
public:
  enum class Status : int {
    OK = TCL_OK, 
    ERROR = TCL_ERROR
  };

  constexpr static Status Error = Status::ERROR;
  constexpr static Status Ok = Status::OK;

  typedef const char * Arg;
  typedef Tcl_CmdProc Cmd;

  Interpreter(Tcl_Interp* interp): m_interp(interp) {}

  inline Status read(Arg arg, int& value) const {
    if (Tcl_GetInt(m_interp, arg, &value) != TCL_OK)
      return Error;
    return Ok;
  }

  inline Status read(Arg arg, double& value) const {
    if (Tcl_GetDouble(m_interp, arg, &value) != TCL_OK)
      return Error;
    return Ok;
  }

  inline Status define(Arg arg, Cmd *cmd, void* ctx) const {
    if (Tcl_CreateCommand(m_interp, arg, cmd, ctx, nullptr) == nullptr)
      return Error;
    return Ok;
  }

  inline Status remove(Arg arg) const {
    if (Tcl_DeleteCommand(m_interp, arg) != TCL_OK)
      return Error;
    return Ok;
  }

private:
  Tcl_Interp* const m_interp;
  StringStream m_stream;
};
#endif
}
