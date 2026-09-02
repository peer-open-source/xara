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
#include "XCP.h"
#include <OpenSeesVersion.h>

namespace Xara {

bool
CheckPolicy(Tcl_Interp* interp, XCP policy)
{
  // 1) 
  OpenSeesVersion version = GetCompatibilityVersion(interp);
  if (version < OpenSeesVersion::X1) {
    return true;
  }
  else
    return false;
}

} // namespace Xara
