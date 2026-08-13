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
  Xara_InterpreterCore
  Xara_CommandInterface
  OPS_Analysis
  OPS_Logging
)

target_link_libraries(XaraCmds_AnalysisRT PRIVATE 
  Xara_InterpreterCore
  OPS_SysOfEqn
  OPS_Logging
)

target_link_libraries(XaraCmds_AnalysisMP PRIVATE 
  Xara_InterpreterCore
  OPS_SysOfEqn
  # OPS_SysOfEqnMP
  OPS_Parallel
  OPS_Logging
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
    sensitivity.cpp
)


target_sources(XaraCmds_AnalysisRT
  PRIVATE
    solver.cpp
    solver.hpp
)

target_sources(XaraCmds_AnalysisMP
  PRIVATE
    solver.cpp
    solver.hpp
)

if (XARA_ENABLE_MPI)
  target_compile_definitions(XaraCmds_AnalysisMP PRIVATE XARA_HAVE_PARALLEL_SOLVERS)
endif()

add_subdirectory(modal)
