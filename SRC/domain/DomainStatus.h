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
enum class DomainStatus: int {
    Success = 0,
    //
    ElementFailedToConverge = -1,
    SectionFailedToConverge = -2,
    MaterialFailedToConverge = -3,
    MaterialWrapperSingular = -5,
    MaterialWrapperFailedToConverge = -6,
    ElementSingular = -4,
};