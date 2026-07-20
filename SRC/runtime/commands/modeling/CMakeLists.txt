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

target_link_libraries(XaraCmds_Modeling PRIVATE Xara_InterpreterCore)

target_sources(XaraCmds_Modeling 
PRIVATE
    element.cpp
    response.cpp
    geomTransf.cpp

    model.cpp
    nodes.cpp
    constraint.cpp
    
    section.cpp
    section/fiber.cpp
    "section/shell.cpp"
    "section/plane.cpp"
    "section/frame.cpp"
    "section/truss.cpp"

    "printing.cpp"

    "utilities/blockND.cpp"
    "utilities/Block2D.cpp"
    "utilities/Block3D.cpp"
)


add_subdirectory(material)
add_subdirectory(element)
add_subdirectory(invoking)
