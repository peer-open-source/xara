
#define ALLOW_IMPLICIT_MATRIX
#include "frames.hpp"
#include <for_int.tpp>
#include <ForceFrame3d.h>


Element*
CreateForceFrame(int tag,
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

  if (sections.size() > MAX_NIP) {
    opserr << OpenSees::PromptValueError
           << "number of sections (" << sections.size() << ") exceeds maximum (" << MAX_NIP << ")\n";
    return nullptr;
  }

  static_loop<0, 3>([&](auto nwm) constexpr {
    if (nwm.value + 6 == ndf) {
      if (!options.shear_flag) {
        static_loop<2,MAX_NIP>([&](auto nip) constexpr {
          if (nip.value == sections.size())
            element = new ForceFrame3d<nip.value, 4+nwm.value*2, nwm.value, 0>(tag, 
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
            element = new ForceFrame3d<nip.value, 6+nwm.value*2, nwm.value, 1>(tag, 
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