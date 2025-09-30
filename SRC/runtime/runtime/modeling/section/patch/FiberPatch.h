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
// File: Patch.h
// Written by Remo M. de Souza
// December 1998
//
#ifndef Patch_h
#define Patch_h
namespace OpenSees {

class FiberCell;

class FiberPatch {
public:
  virtual int getMaterialID() const = 0;
  virtual int getNumCells() const   = 0;
  virtual FiberCell** getCells() const   = 0;
};
} // namespace OpenSees
#endif
