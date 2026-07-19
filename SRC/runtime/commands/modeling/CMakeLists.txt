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

target_link_libraries(OPS_Runtime PRIVATE Xara_InterpreterCore)

target_sources(OPS_Runtime PRIVATE
    element.cpp
    response.cpp
    geomTransf.cpp
    section/fiber.cpp
)


add_subdirectory(material)
add_subdirectory(element)
add_subdirectory(invoking)
