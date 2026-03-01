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
#include <MixedFrame3d02.h>


Element*
CreateMixedFrame(int tag,
                 int ndf,
                 std::array<int,2>& nodes,
                 std::vector<FrameSection*>& sections,
                 BeamIntegration& beamIntegr,
                 FrameTransformBuilder& tb,
                 Options& options,
                 double mass, int max_iter, double tol)
{
  bool use_mass = false;

  Element* element = nullptr;

  static_loop<0, 3>([&](auto nwm) constexpr {
    if (nwm.value + 6 == ndf) {
      if (!options.shear_flag) {
        static_loop<2,MAX_NIP>([&](auto nip) constexpr {
          if (nip.value == sections.size())
            element = new MixedFrame3d02<nip.value, 4+nwm.value*2, nwm.value, 0>(tag, 
                                          nodes,
                                          sections,
                                          beamIntegr, tb,
                                          mass, options.mass_flag, use_mass,
                                          max_iter, tol
                                          );
          });
      }
      else {
        static_loop<2,MAX_NIP>([&](auto nip) constexpr {
          if (nip.value == sections.size())
            element = new MixedFrame3d02<nip.value, 6+nwm.value*2, nwm.value, 1>(tag, 
                                          nodes,
                                          sections,
                                          beamIntegr, tb,
                                          mass, options.mass_flag, use_mass,
                                          max_iter, tol
                                          );
        });
      }
    }
  });

  return element;
}