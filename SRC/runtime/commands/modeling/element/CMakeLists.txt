#===----------------------------------------------------------------------===//
#
#                                   xara
#                              https://xara.so
#
#===----------------------------------------------------------------------===//
#
# Copyright (c) 2025, OpenSees/Xara Developers
# All rights reserved.  No warranty, explicit or implicit, is provided.
#
# This source code is licensed under the BSD 2-Clause License.
# See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
#
#===----------------------------------------------------------------------===//

target_sources(XaraCmds_Modeling
PRIVATE
    frames/prism.cpp
    frames/frames.cpp
    frames/exact.cpp
    # frames/mixed.cpp
    frames/force.cpp
    shells.cpp
    brick.cpp
    plane.cpp
    truss.cpp
    # plate.cpp
    tests.cpp
)