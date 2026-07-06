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
namespace OpenSees {
enum class DomainStatus: int {
    Success = 0,
    //
    ElementFailedToConverge = -100,
    SectionFailedToConverge = -200,
    MaterialFailedToConverge = -300,
    MaterialWrapperSingular = -500,
    MaterialWrapperFailedToConverge = -600,
    ElementSingular = -400,
};
}
