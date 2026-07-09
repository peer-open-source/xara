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

#include <array>
#include <Parsing.h>
class Element;
class SectionForceDeformation;

Element* CreatePlateQ4(TCL_Char *name,
                       int tag,
                       const std::array<int,4>& nodes,
                       SectionForceDeformation& section);
