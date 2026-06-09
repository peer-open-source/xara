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
#define ALLOW_IMPLICIT_MATRIX
#include "frames.hpp"
#include <for_int.tpp>
#include <ExactFrame3d.h>
// #define XARA_ExactFrame02

#ifdef XARA_ExactFrame02
#include <ExactFrame02.h>
#endif
// #include <ExactFrame03.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
class CrdTransf;

#define MAX_NEN 6

Element*
CreateExactFrame(int tag,
                 int ndf,
                 const std::vector<int>& nodev,
                 std::vector<FrameSection*>& sections,
                 BeamIntegration& beamIntegr,
                 CrdTransf& theTransf,
                 const FrameOptions& options,
                 double mass
                )
{
  // bool use_mass = false;

  Element* element = nullptr;

  if (!options.shear_flag) {
      opserr << OpenSees::PromptValueError 
              << "ExactFrame3d requires shear formulation"
              << OpenSees::SignalMessageEnd;
      return nullptr;
  }
  int exact_version = 0;
  if (getenv("ExactFrame")) {
    exact_version = atoi(getenv("ExactFrame"));
  }
  if ((options.rotation_type != Rotations::Parameters::Iter) && 
      (options.rotation_type != Rotations::Parameters::None) &&
      (options.rotation_type != Rotations::Parameters::Init)) {
    exact_version = 2;
  }

  if (sections.size() < nodev.size()-1)
    for (unsigned i = 0; i < nodev.size()-1; ++i)
      sections.push_back(sections[0]);

  unsigned nen = nodev.size();
  static_loop<2,MAX_NEN>([&](auto nn) constexpr {
    if (nn.value == nen) {
      std::array<int, nn.value> nodes;
      std::copy_n(nodev.begin(), nn.value, nodes.begin());
      static_loop<0,2>([&](auto nwm) constexpr {
        if (nwm.value+6 == ndf) {
          // element = new ExactFrame3d<nn.value, nwm.value>(tag, nodes, sections.data(), theTransf);
          if (!exact_version) {
            element = new ExactFrame3d<nn.value, nwm.value>(tag, nodes, sections.data(), theTransf);
            return;
          }
#ifdef XARA_ExactFrame02
          else if (exact_version == 2) {
            element = new ExactFrame02<nn.value, nwm.value>(tag, nodes, sections.data(), theTransf);
            return;
          }
#endif
        }
      });
    }
  });

  if (element == nullptr) {
      opserr << OpenSees::PromptValueError 
              << "invalid number of dofs for ExactFrame; got " << ndf 
              << OpenSees::SignalMessageEnd;
      return nullptr;
  }

  return element;
}