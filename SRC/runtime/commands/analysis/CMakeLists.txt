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
target_link_libraries(XaraCmds_Analysis PRIVATE 
  Xara_CommandInterface
  OPS_Analysis
  OPS_SysOfEqn
  OPS_Parallel
)

target_include_directories(XaraCmds_Analysis PRIVATE 
  # library
  ${OPS_SRC_DIR}/analysis/algorithm/equiSolnAlgo
  # interpreter
  ${OPS_SRC_DIR}/runtime/interpreter
)

target_sources(XaraCmds_Analysis
  PRIVATE
    algorithm.cpp
    integrator.cpp
    transient.cpp
    analysis.cpp
    numberer.cpp
    ctest.cpp
    solver.cpp
    solver.hpp
    sensitivity.cpp
)

add_subdirectory(modal)
