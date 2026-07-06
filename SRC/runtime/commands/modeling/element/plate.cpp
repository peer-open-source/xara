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

#include <Parsing.h>
#include <element/Plate/PlateQ4.hpp>

#ifdef _MSC_VER 
#  include <string.h>
#  define strcasecmp _stricmp
#else
#  include <strings.h>
#endif

using namespace OpenSees;

Element*
CreatePlateQ4(TCL_Char * name,
              int tag,
              const std::array<int,4>& nodes, 
              SectionForceDeformation& section
            )
{
  if (strcasecmp(name, "PlateQ4/F") == 0) {
    return new PlateQ4_Uniform(tag, nodes, section);
  }
  else if (strcasecmp(name, "PlateQ4/U") == 0) {
    // Uniformly reduced integration (URI)
    return new PlateQ4_T<Membrane::Reduced,PlateReduced,DrillShellMITC4Penalty>(tag, nodes, section);
  }
  else if (strcasecmp(name, "PlateQ4/S") == 0) {
    return new PlateQ4_SRI(tag, nodes, section);
  }
  else if (strcasecmp(name, "PlateQ4/L01") == 0) {
    return new PlateQ4_MITC4(tag, nodes, section);
  }
  else if (strcasecmp(name, "PlateQ4/L02") == 0) {
    return new PlateQ4_AGQI_MITC4(tag, nodes, section);
  }
  else if (strcasecmp(name, "PlateQ4/E5") == 0) {
    using PlateType = PlateQ4_T<Membrane::EAS<Membrane::EnhancedQuadMembraneInterpolation>,
                                PlateMITC4,
                                DrillShellMITC4Penalty>;
    return new PlateType(tag, nodes, section);
  }
  return nullptr;
}