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

target_sources(XaraCmds_Modeling PRIVATE 
  material.cpp
  nDMaterial.cpp
  uniaxial.cpp
  fedeas.cpp
  legacy.cpp
  wrapper.cpp
  plastic.cpp
  elastic.cpp
  concrete.cpp
  concrete_asd.cpp
  isotropy.cpp
  boucwen.cpp
)