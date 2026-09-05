//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2026, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
// Description: This file declares "Xara Command Policies" (XCP) which are used
// to enable backwards-compatible command behaviors.
//
#include <Parsing.h>
#include <unordered_map>

namespace Xara {

enum class XCP {
  XCP0001 = 1, // Override 0 mass with material density
  XCP0002 = 2, // Switch order of thickness and material for SSPquad

  XCP0003 = 3, // Use old determinant switching logic for static integrators

  XCP0004 = 4, // Use inconsistent tangent in PDelta
};


bool CheckPolicy(Tcl_Interp* interp, XCP policy);
}


