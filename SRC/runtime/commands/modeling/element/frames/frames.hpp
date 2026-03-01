
// Maximum number of integration points. 
// Large values can noticeably impact compile time.
#ifdef _DEBUG
 #define MAX_NIP 8
#elif defined(XARA_RELEASE)
 #define MAX_NIP 20
#else
 #define MAX_NIP 15
#endif
// #define MAX_NIP 6
#include <array>
#include <vector>

#include <FrameSection.h>
class Element;
class BeamIntegration;
namespace OpenSees {
class FrameTransformBuilder;
}

using namespace OpenSees;
struct Options {
  int mass_flag;
  int use_mass;
  int shear_flag;
  int geom_flag;
};


Element*
CreateMixedFrame(int tag,
                 int ndf,
                 std::array<int,2>& nodes,
                 std::vector<FrameSection*>& sections,
                 BeamIntegration& beamIntegr,
                 FrameTransformBuilder& tb,
                 Options& options,
                 double mass, int max_iter, double tol);


Element*
CreateForceFrame(int tag,
                 int ndf,
                 std::array<int,2>& nodes,
                 std::vector<FrameSection*>& sections,
                 BeamIntegration& beamIntegr,
                 FrameTransformBuilder& tb,
                 Options& options,
                 double mass, int max_iter, double tol);