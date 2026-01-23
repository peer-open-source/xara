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
// Written: cmp, mhs, rms, fmk
//
// Created: Feb 2023
//
#define ALLOW_IMPLICIT_MATRIX
// Standard library
#include <string>
#include <array>
#include <algorithm>
#include <vector>
#include <utility>
#include <string.h>
#include <assert.h>
#include <math.h>
#ifdef _MSC_VER 
#  include <string.h>
#  define strcasecmp _stricmp
#else
#  include <strings.h>
#endif
#define strcmp strcasecmp

// Maximum number of integration points. 
// Large values can noticeably impact compile time.
#ifdef _DEBUG
 #define MAX_NIP 8
#elif defined(XARA_RELEASE)
 #define MAX_NIP 30
#else
 #define MAX_NIP 15
#endif

// Parsing
#include <Logging.h>
#include <Parsing.h>
#include <ArgumentTracker.h>

// Model
#include <Node.h>
#include <Domain.h>
#include <ModelRegistry.h>

// Sections
#include <FrameSection.h>
#include <ElasticSection2d.h>
#include <ElasticSection3d.h>

// Elements
#include "ElasticBeam2d.h"
#include "ElasticBeam2d.h"
#include "ElasticBeam3d.h"
#include "ElasticBeam3d.h"
#include "PrismFrame2d.h"
#include "PrismFrame2d.h"
#include "PrismFrame3d.h"
#include "PrismFrame3d.h"

#include <CubicFrame3d.h>
#include <ForceFrame3d.h>
#include <ForceDeltaFrame3d.h>
#include <EulerFrame3d.h>
#include <EulerDeltaFrame3d.h>
#include <ExactFrame3d.h>
#include <ShearFrame3d.h>

#include <DispBeamColumn2d.h>
#include <DispBeamColumn2dThermal.h>
#include <DispBeamColumn3d.h>
#include <DispBeamColumnAsym3d.h>
#include <DispBeamColumn3dThermal.h>
#include <DispBeamColumnNL2d.h>

#include <ElasticForceBeamColumn2d.h>
#include <ElasticForceBeamColumn3d.h>
#include <ElasticForceBeamColumnWarping2d.h>

#include <ForceBeamColumn2d.h>
#include <ForceBeamColumn2d.h>
#include <ForceBeamColumn2dThermal.h>
#include <ForceBeamColumn3d.h>
#include <ForceBeamColumnCBDI2d.h>
#include <ForceBeamColumnCBDI3d.h>
#include <ForceBeamColumnWarping2d.h>
#include <TimoshenkoBeamColumn2d.h>

// Quadrature
#include <BeamIntegration.h>
#include <LobattoBeamIntegration.h>
#include <LegendreBeamIntegration.h>
#include <HingeEndpointBeamIntegration.h>
#include <HingeMidpointBeamIntegration.h>
#include <HingeRadauBeamIntegration.h>
#include <HingeRadauTwoBeamIntegration.h>

#include <transform/FrameTransformBuilder.hpp>

using namespace OpenSees;

struct Options {
  int mass_flag;
  int use_mass;
  int shear_flag;
  int geom_flag;
};

enum class FrameClass {
  ForceFrame3d,
  ForceDeltaFrame3d,
  EulerFrame3d,
  EulerDeltaFrame3d,
  ExactFrame3d,
  CubicFrame3d,
  // Legacy
  DispBeamColumn2d,
  DispBeamColumn3d,
  DispBeamColumnNL2d,
  DispBeamColumn2dThermal,
  DispBeamColumn3dThermal,
  DispBeamColumnAsym3d,
  ElasticForceBeamColumn2d,
  ElasticForceBeamColumn3d,
  ElasticForceBeamColumnWarping2d,
  ForceBeamColumn2d,
  ForceBeamColumn3d,
  ForceBeamColumnCBDI2d,
  ForceBeamColumnCBDI3d,
  ForceBeamColumnWarping2d,
  ForceBeamColumn2dThermal,
  TimoshenkoBeamColumn2d,
  Unknown
};


extern BeamIntegration*     GetBeamIntegration(TCL_Char* type, int);
extern BeamIntegrationRule* GetHingeStencil(int argc, TCL_Char ** const argv);

static inline int
CheckTransformation(Domain& domain, int iNode, int jNode, CrdTransf& transform)
{
  Node* ni = domain.getNode(iNode);
  Node* nj = domain.getNode(jNode);
  if (ni == nullptr || nj == nullptr) {
    opserr << OpenSees::PromptValueError << "nodes not found with tags "
           << iNode << " and " << jNode
           << OpenSees::SignalMessageEnd;
  }

  if (transform.initialize(ni, nj) != 0) {
    if (transform.getInitialLength() <= 0.0) {
      opserr << OpenSees::PromptValueError 
            << "element has zero or negative initial length "
            << transform.getInitialLength()
            << "; check for duplicate nodes"
            << OpenSees::SignalMessageEnd;
    }
    else {
      opserr << OpenSees::PromptValueError 
            << "transformation with tag " << transform.getTag()
            << " could not be initialized with nodes "
            << iNode << " and " << jNode
            << "; check orientation"
            << OpenSees::SignalMessageEnd;
    }
    return TCL_ERROR;
  }
  return TCL_OK;
}



template <int ndm, typename Transform, typename Section>
static Element*
CreateFrame(ModelRegistry& builder, 
            const char* name,
            int tag,
            std::vector<int>& nodev,
            int transfTag,
            const std::vector<int>& section_tags,
            BeamIntegration& beamIntegr,
            double mass, int max_iter, double tol,
            const std::array<double,2>& shear_center,
            Options& options) 
{
  std::vector<Section*> sections;

  // Finalize sections
  assert(section_tags.size() != 0);
  for (int tag : section_tags) {
    Section *section = builder.getTypedObject<Section>(tag);
    if (section == nullptr)
      return nullptr;
    sections.push_back(section);
  }
  int nIP = sections.size();

  SectionForceDeformation** secptrs = (SectionForceDeformation**)(sections.data());

  if (options.shear_flag == -1) {
    options.shear_flag = 0;
    const ID& resultants = sections[0]->getType();
    for (int i=0; i< sections[0]->getOrder(); i++)
      if (resultants(i) == FrameStress::Vy)
        options.shear_flag = 1;
  }

  // Finalize the coordinate transform
  CrdTransf* theTransf = builder.getTypedObject<CrdTransf>(transfTag);
  if (theTransf == nullptr) {
    opserr << OpenSees::PromptValueError 
           << "transformation not found with tag " 
           << transfTag 
           << OpenSees::SignalMessageEnd;
    return nullptr;
  }

  //
  // Create the element
  //
  Element  *theElement = nullptr;

  int iNode = nodev[0],
      jNode = nodev[1];
  bool use_mass = options.use_mass;

  if constexpr (ndm == 2) {

    if (strcmp(name, "elasticForceBeamColumn") == 0)
      theElement = new ElasticForceBeamColumn2d(tag, iNode, jNode, nIP, secptrs, beamIntegr, *theTransf, mass);
    else if (strcmp(name, "timoshenkoBeamColumn") == 0)
      theElement =
          new TimoshenkoBeamColumn2d(tag, iNode, jNode, nIP, secptrs, beamIntegr, *theTransf, mass);
    else if (strcmp(name, "dispBeamColumn") == 0)
      theElement =
          new DispBeamColumn2d(tag, iNode, jNode, nIP, secptrs, beamIntegr, *theTransf, mass, options.mass_flag);
    else if (strcmp(name, "dispBeamColumnNL") == 0)
      theElement =
          new DispBeamColumnNL2d(tag, iNode, jNode, nIP, secptrs,
                                 beamIntegr, *theTransf, mass);

    else if (strcmp(name, "dispBeamColumnThermal") == 0)
      theElement = new DispBeamColumn2dThermal(tag, iNode, jNode, nIP, secptrs, beamIntegr, *theTransf, mass);


    else if (strcmp(name, "dispBeamColumnWithSensitivity") == 0)
      theElement = new DispBeamColumn2d(tag, iNode, jNode, nIP, secptrs, beamIntegr, *theTransf, mass);


    // Force formulations
    else if (strcmp(name, "forceBeamColumnCBDI") == 0)
      theElement = new ForceBeamColumnCBDI2d(tag, iNode, jNode, nIP, secptrs,
                                             beamIntegr, *theTransf, mass, false);

    else if (strcmp(name, "forceBeamColumnCSBDI") == 0)
      return  new ForceBeamColumnCBDI2d(tag, iNode, jNode, nIP, secptrs,
                                              beamIntegr, *theTransf, mass, true);

    else if (strcmp(name, "forceBeamColumnWarping") == 0)
      return
          new ForceBeamColumnWarping2d(tag, iNode, jNode, nIP,
                                       secptrs, beamIntegr, *theTransf);

    else if (strcmp(name, "forceBeamColumnThermal") == 0)
      return new ForceBeamColumn2dThermal(tag, iNode, jNode, nIP, secptrs, beamIntegr, *theTransf, mass);

    else if (strcmp(name, "elasticForceBeamColumnWarping") == 0)
      return new ElasticForceBeamColumnWarping2d(
          tag, iNode, jNode, nIP, secptrs, beamIntegr, *theTransf);
    else
      theElement =
          new ForceBeamColumn2d(tag, iNode, jNode, nIP, secptrs,
                                beamIntegr, *theTransf, mass, max_iter, tol);


  } else {
    //
    // ndm == 3
    //

    if (CheckTransformation(*builder.getDomain(), nodev[0], nodev[nodev.size()-1], *theTransf) != TCL_OK)
      return nullptr;
    
    // Create 3d frame elements
    if (strstr(name, "Frame") != nullptr) {
      if (strstr(name, "Exact") == nullptr) {

        std::array<int, 2> nodes {nodev[0], nodev[1]};

        FrameTransformBuilder* tb = builder.getTypedObject<FrameTransformBuilder>(transfTag);

        if (!tb) {
          opserr << OpenSees::PromptValueError << "invalid transform\n";
          return nullptr;
        }


        if (strcmp(name, "EulerFrame") == 0) {
            theElement = new EulerFrame3d(tag, nodes, nIP, 
                                          sections.data(),
                                          beamIntegr, 
                                          *tb, 
                                          mass, options.mass_flag);
        }

        else if (strcmp(name, "CubicFrame") == 0) {
          if (options.shear_flag)
            theElement = new CubicFrame3d<true,0>(tag, nodes, 
                                          sections,
                                          beamIntegr, 
                                          *theTransf, // TODO: Use FrameTransformBuilder
                                          mass);
          else
            theElement = new CubicFrame3d<false,0>(tag, nodes, 
                                          sections,
                                          beamIntegr, 
                                          *theTransf, // TODO: Use FrameTransformBuilder
                                          mass);
        } 

        else if (strcmp(name, "DisplFrame") == 0) {
          theElement =  new EulerDeltaFrame3d(tag, nodes, sections,
                                              beamIntegr, *theTransf, 
                                              mass, 
                                              options.mass_flag, 
                                              use_mass);
        }

        else if (strcmp(name, "ShearFrame") == 0) {
          if (!options.shear_flag) {
            opserr << OpenSees::PromptValueError 
                  << "ShearFrame3d requires shear"
                  << OpenSees::SignalMessageEnd;
            return nullptr;
          }

          int ndf = builder.getNDF();
          if (sections.size() < nodev.size()-1)
            for (unsigned i = 0; i < nodev.size()-1; ++i)
              sections.push_back(sections[0]);
          
          unsigned nen = nodev.size();
          static_loop<2,6>([&](auto nn) constexpr {
            if (nn.value == nen) {
              std::array<int, nn.value> nodes;
              std::copy_n(nodev.begin(), nn.value, nodes.begin());
              static_loop<0,4>([&](auto nwm) constexpr {
                if (nwm.value+6 == ndf)
                  theElement = new ShearFrame3d<nn.value, nwm.value>(tag, nodes, sections.data(), *theTransf);
              });
            }
          });
          if (theElement == nullptr) {
            opserr << OpenSees::PromptValueError 
                  << "invalid number of dofs for ShearFrame; got " << ndf 
                  << OpenSees::SignalMessageEnd;
            return nullptr;
          }
        }
        else if ((strstr(name, "Force") != 0) ||
                (strcmp(name, "MixedFrame") == 0)) {
          if (strcmp(name, "ForceDeltaFrame") == 0 || options.geom_flag) {
            if (sections.size() > MAX_NIP) {
              opserr << OpenSees::PromptValueError 
                     << "too many sections for ForceDeltaFrame3d: " << static_cast<int>(sections.size())
                     << OpenSees::SignalMessageEnd;
              return nullptr;
            }
            if (!options.shear_flag)
              static_loop<2,MAX_NIP>([&](auto nip) constexpr {
                if (nip.value == sections.size())
                  theElement = new ForceDeltaFrame3d<nip.value, 4>(tag, nodes, sections,
                                                beamIntegr, *tb, 
                                                mass, 
                                                options.mass_flag, 
                                                use_mass,
                                                max_iter, tol,
                                                options.shear_flag
                                                );
              });
            else
              static_loop<2,MAX_NIP>([&](auto nip) constexpr {
                if (nip.value == sections.size())
                  theElement = new ForceDeltaFrame3d<nip.value, 6>(tag, nodes, sections,
                                                beamIntegr, *tb, 
                                                mass, 
                                                options.mass_flag, 
                                                use_mass,
                                                max_iter, tol,
                                                options.shear_flag
                                                );
              });
          } else {
            int ndf = builder.getNDF();
            if (sections.size() > MAX_NIP) {
              opserr << OpenSees::PromptValueError 
                     << "too many sections for element: " 
                     << static_cast<int>(sections.size())
                     << OpenSees::SignalMessageEnd;
              return nullptr;
            }
            static_loop<0, 3>([&](auto nwm) constexpr {
              if (nwm.value + 6 == ndf) {
                if (!options.shear_flag) {
                  static_loop<2,MAX_NIP>([&](auto nip) constexpr {
                    if (nip.value == sections.size())
                      theElement = new ForceFrame3d<nip.value, 4+nwm.value*2, nwm.value, 0>(tag, 
                                                    nodes,
                                                    sections,
                                                    beamIntegr, *tb,
                                                    mass, options.mass_flag, use_mass,
                                                    max_iter, tol
                                                    );
                    });
                }
                else {
                  static_loop<2,MAX_NIP>([&](auto nip) constexpr {
                    if (nip.value == sections.size())
                      theElement = new ForceFrame3d<nip.value, 6+nwm.value*2, nwm.value, 1>(tag, 
                                                    nodes,
                                                    sections,
                                                    beamIntegr, *tb,
                                                    mass, options.mass_flag, use_mass,
                                                    max_iter, tol
                                                    );
                  });
                }
              }
            });
          }
        }
      }

      else if (strcmp(name, "ExactFrame") == 0) {
        if (!options.shear_flag) {
          opserr << OpenSees::PromptValueError 
                 << "ExactFrame3d requires shear formulation"
                 << OpenSees::SignalMessageEnd;
          return nullptr;
        }

        int ndf = builder.getNDF();
        if (sections.size() < nodev.size()-1)
          for (unsigned i = 0; i < nodev.size()-1; ++i)
            sections.push_back(sections[0]);
        
        unsigned nen = nodev.size();
        static_loop<2,6>([&](auto nn) constexpr {
          if (nn.value == nen) {
            std::array<int, nn.value> nodes;
            std::copy_n(nodev.begin(), nn.value, nodes.begin());
            static_loop<0,4>([&](auto nwm) constexpr {
              if (nwm.value+6 == ndf)
                theElement = new ExactFrame3d<nn.value, nwm.value>(tag, nodes, sections.data(), *theTransf);
            });
          }
        });
        if (theElement == nullptr) {
          opserr << OpenSees::PromptValueError 
                 << "invalid number of dofs for ExactFrame; got " << ndf 
                 << OpenSees::SignalMessageEnd;
          return nullptr;
        }
      }
    }

    else if (strcmp(name, "elasticForceBeamColumn") == 0)
      theElement = new ElasticForceBeamColumn3d(tag, iNode, jNode, nIP, secptrs, 
                                                beamIntegr, *theTransf, mass);

    else if (strcasecmp(name, "dispBeamColumn") == 0)
      theElement = new DispBeamColumn3d(tag, iNode, jNode, nIP, secptrs,
                                        beamIntegr, *theTransf, 
                                        mass, options.mass_flag);

    else if (strcasecmp(name, "dispBeamColumnWithSensitivity") == 0)
      theElement = new DispBeamColumn3d(
          tag, iNode, jNode, nIP, secptrs, beamIntegr, *theTransf, mass);

    else if (strcasecmp(name, "dispBeamColumnThermal") == 0)
      theElement = new DispBeamColumn3dThermal(
          tag, iNode, jNode, nIP, secptrs, beamIntegr, *theTransf, mass);

    else if (strcasecmp(name, "forceBeamColumnCBDI") == 0)
      theElement = new ForceBeamColumnCBDI3d(tag, iNode, jNode, nIP, secptrs,
                                             beamIntegr, *theTransf, 
                                             mass, false, max_iter, tol);
    else if (strcasecmp(name, "dispBeamColumnAsym") == 0)
      theElement = new DispBeamColumnAsym3d(tag, iNode, jNode, nIP, secptrs,
                                            beamIntegr, *theTransf, 
                                            shear_center[0], shear_center[1],
                                            mass, options.mass_flag);
    else
      theElement = new ForceBeamColumn3d(tag, iNode, jNode, nIP, secptrs,
                                         beamIntegr, *theTransf, mass, max_iter, tol);
  }
  return theElement;
}


// 0       1    2 3  4
// element beam 1 $i $j 0 1 2
//
// Legacy:
//  a) Gauss Inline
//     0       1     2    3    4    5    6
//                                  0    1
//     element $type $tag $ndi $ndj $trn "Gauss arg1 arg2 ..." 
//             <-mass $mass> <-iter $iter $tol>
//
//  b) Gauss Predefined
//                                  0    1
//     element(type, tag, ndi, ndj, trn, itag, 
//              iter=(10, 1e-12), mass=0.0)
//
//  c) "Original/Obsolete": Gauss Tabulated, Prismatic
//                                  0    1    2
//     element $type $tag $ndi $ndj $nip $sec $trn 
//             <-mass $mass> <-iter $iter $tol> <-integration $ityp>
//
//  d) Gauss Tabulated, Non-Prismatic
//                                  0    1         ... (2 + nIP)
//     element $type $tag $ndi $ndj $nip -sections ... $trn 
//             <-mass $massDens> <-cMass> <-integration $ityp>
//
//   or
//                                   0
//     element $type $tag $ndi $ndj  $trn 
//             -sections {...}
//             <-mass $massDens> <-cMass> <-integration $ityp>
//
//
// Integration may be specitied as either 
//   i  ) a single name, 
//   ii ) a pattern spec, or 
//   iii) a tag for a pattern
// 
// if a list of sections is given with -sections, or nIP is provided, then 
// we must have an integration with the form (i)
//
// 1) Parse common keyword args: 
//      "-mass" $mass, 
//      "-cMass"/"-lMass"
//      "-mass-form" $form
//
//      "-iter" $iter $tol, 
//      
//      "-integration" $Integration
//        - first try parsing $Integration as integer ($itag, form (iii))
//          if successfull, populate section_tags and continue
//        - next try parsing $Integration as basic quadrature ($ityp, form (i))
//        - finally, try parsing as full pattern spec
//
//      "-section" $Tag
//
//      "-transform" $Tag
//      "-vertical" {}
//
//      "-sections" $SectionTags
//        - if cannot split $SectionTags as list, then
//          mark "-sections" as positional and continue
//          with keyword loop
//        - Check if "-integration" was provided already; if so, it must have been in form (i);
//          otherwise throw an error.
//        - Parse $SectionTags
//          -sections {...} may occur anywhere
//          -sections ...   must occur after nIP is obtained
//
//
//  2) 
//     If pos[1] == "-sections" then command is Form (d):
//        nIP = pos[0]
//        trn = pos[2+nIP]
//
//     else
//
//     switch (pos.size())
//     case 1: // Form (e)
//        trn = pos[0]
//        if (section_tags.size() == 0)
//          ERROR
//
//     case 2: 
//        // Form (a) or (b)
//        trn = pos[0]
//        if GetInt(interp, pos[1]):
//           itag = pos[1]
//        else
//           ParseHingeScheme(pos[1])
//
//      case 3: 
//         // Form (c)
//         nip = int(pos[0])
//         sec = int(pos[1])
//         trn = int(pos[2])
//
int
TclBasicBuilder_addForceBeamColumn(ClientData clientData, Tcl_Interp *interp,
                                   int argc, TCL_Char **const argv)
{
  assert(clientData != nullptr);
  ModelRegistry *builder = static_cast<ModelRegistry*>(clientData);
  Domain *domain = builder->getDomain();
  assert(domain != nullptr);

  int status = TCL_OK;

  enum class FrameSyntax {
    None,
    A, // _GaussInline,
    B, // _GaussLookup,
    C, // _GaussPrismRule,
    D, // _GaussHingeRule,
    X
  } syntax = FrameSyntax::None;

  struct {
    bool transform   = false;
    bool integration = false;
    bool sections    = false;
  } received;

  // collect positional arguments
  std::vector<int> positions;


  //
  // Preliminary checks
  //
  int ndm = builder->getNDM();
  int ndf = builder->getNDF();

  { // Check dimension and DOFs of problem
    int ok = 0;
    if ((ndm == 2 && ndf == 3) || (ndm == 2 && ndf == 4))
      ok = 1;
    if (ndm == 3 && ndf >= 6)
      ok = 1;

    if (ok == 0) {
      opserr << OpenSees::PromptValueError << "ndm = " << ndm << " and ndf = " << ndf
             << " not compatible with Frame element" 
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }
  }

  if (argc < 6) {
    opserr << OpenSees::PromptValueError << "insufficient arguments\n";
    return TCL_ERROR;
  }

  //
  // Required positional arguments
  //
  FrameClass beam_type = FrameClass::Unknown;
  if (strcasecmp(argv[1], "dispBeamColumn") == 0) {
    if (ndm == 2)
      beam_type = FrameClass::DispBeamColumn2d;
    else
      beam_type = FrameClass::DispBeamColumn3d;
  }
  else if ((strcasecmp(argv[1], "forceBeamColumn") == 0) ||
           (strcasecmp(argv[1], "nonlinearBeamColumn") == 0)) {
    if (ndm == 2)
      beam_type = FrameClass::ForceBeamColumn2d;
    else
      beam_type = FrameClass::ForceBeamColumn3d;
  }
  else if (strcasecmp(argv[1], "timoshenkoBeamColumn") == 0) {
    if (ndm == 2)
      beam_type = FrameClass::TimoshenkoBeamColumn2d;
    else {
      opserr << OpenSees::PromptValueError 
             << "timoshenkoBeamColumn not available in 3d\n"
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }
  }
  else if (strcasecmp(argv[1], "ForceFrame") == 0) {
    if (ndm == 2)
      beam_type = FrameClass::ForceBeamColumn2d;
    else
      beam_type = FrameClass::ForceFrame3d;
  }
  // else {
  //   opserr << OpenSees::PromptValueError 
  //          << "unknown beam element type " << argv[1] 
  //          << OpenSees::SignalMessageEnd;
  //   return TCL_ERROR;
  // }

  int tag;
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "invalid " << " tag " << tag 
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  int argi = 5;
  std::vector<int> multi_nodes;
  {
    int list_argc;
    TCL_Char **list_argv;
    if (Tcl_SplitList(interp, argv[3], &list_argc, &list_argv) == TCL_OK && list_argc >= 2) {
      argi -= 1;
      
      for (int i = 0; i < list_argc; ++i) {
        int node;
        if (Tcl_GetInt(interp, list_argv[i], &node) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid node\n";
          return TCL_ERROR;
        }
        multi_nodes.push_back(node); 
      }
      Tcl_Free((char *)list_argv);
    }
    else {
      int iNode, jNode;
      if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid iNode\n";
        return TCL_ERROR;
      }
  
      if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid jNode\n";
        return TCL_ERROR;
      }
      multi_nodes.push_back(iNode);
      multi_nodes.push_back(jNode);
    }
  }

  int max_iter = 10;
  double tol  = 1.0e-12;
  double mass = 0.0;
  int transfTag;
  std::array<double,2> shear_center = {0.0, 0.0};

  //
  // Defaults
  //
  struct Options options;
  options.mass_flag  =  0;
  options.shear_flag = -1;
  options.geom_flag  =  0;
  options.use_mass   =  0;
  switch (beam_type) {
    case FrameClass::DispBeamColumn2d:
    case FrameClass::DispBeamColumn3d:
    case FrameClass::ForceBeamColumn2d:
    case FrameClass::ForceBeamColumn3d:
      options.shear_flag = 0;
      break;
    case FrameClass::TimoshenkoBeamColumn2d:
      options.shear_flag = 1;
      break;
    default:
      break;
  }

  //
  // Parse keywords
  //
  {
    while (argi < argc) {
      // Shear
      if (strcmp(argv[argi], "-shear") == 0) {
        if (argc < argi + 2) {
          opserr << OpenSees::PromptValueError << "not enough arguments, expected -shear $flag\n";
          return TCL_ERROR;
        }
        if (Tcl_GetInt(interp, argv[argi + 1], &options.shear_flag) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid shear_flag, expected integer\n";
          return TCL_ERROR;
        }
        argi += 2;
      }

      // Geometry
      else if (strcmp(argv[argi], "-order") == 0) {
        syntax = FrameSyntax::X;
        if (argc < argi + 2) {
          opserr << OpenSees::PromptValueError << "not enough arguments, expected -order $flag\n";
          return TCL_ERROR;
        }
        if (Tcl_GetInt(interp, argv[argi + 1], &options.geom_flag) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid geom_flag, expected integer\n";
          return TCL_ERROR;
        }
        argi += 2;
      }

      // -iter $max_iter $tol 
      else if (strcmp(argv[argi], "-iter") == 0) {
        if (argc < argi + 3) {
          opserr << OpenSees::PromptValueError << "not enough -iter args need -iter max_iter? tol?\n";
          return TCL_ERROR;
        }
        if (Tcl_GetInt(interp, argv[argi + 1], &max_iter) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid max_iter\n";
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[argi + 2], &tol) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid tol\n";
          return TCL_ERROR;
        }
        argi += 3;
      }
      // mass
      else if (strcmp(argv[argi], "-mass") == 0) {
        if (argc < argi + 2) {
          opserr << OpenSees::PromptValueError << "not enough arguments, expected -mass $mass\n";
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[argi + 1], &mass) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid mass\n";
          return TCL_ERROR;
        }
        argi += 2;
        options.use_mass = 1;

      // mass type
      } else if ((strcmp(argv[argi], "-lMass") == 0) ||
                 (strcmp(argv[argi], "lMass") == 0)) {
        options.mass_flag = 0;
        argi++;
      }
      else if ((strcmp(argv[argi], "-cMass") == 0) ||
                 (strcmp(argv[argi], "cMass") == 0)) {
        options.mass_flag = 1;
        argi++;
      }

      // Transform
      else if (strcmp(argv[argi], "-transform") == 0) {
        syntax = FrameSyntax::X;
        if (argc < argi + 2) {
          opserr << OpenSees::PromptValueError 
                 << "not enough arguments, expected -transform $transform\n";
          return TCL_ERROR;
        }

        argi++;
        if (Tcl_GetInt(interp, argv[argi], &transfTag) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid transform\n";
          return TCL_ERROR;
        }
        argi++;
        received.transform = true;
      }

      else if (strcmp(argv[argi], "-shearCenter") == 0) {
        if (argc < argi + 3) {
          opserr << OpenSees::PromptValueError
                 << "not enough arguments, expected -shearCenter $y $z"
                 << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[argi + 1], &shear_center[0]) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid shear_center y\n";
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[argi + 2], &shear_center[1]) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid shear_center z\n";
          return TCL_ERROR;
        }
        argi += 3;
      }

      // Quadrature
      else if (strcmp(argv[argi], "-integration") == 0) {
        // Either syntax c, d, or e.
        if (syntax == FrameSyntax::X) {
          opserr << OpenSees::PromptValueError
                 << "-integration cannot be used with Xara options\n";
          return TCL_ERROR;
        }
        syntax = FrameSyntax::C;
        for (int ii=0; ii<argc; ii++) {
          if (strcmp(argv[ii], "-sections") == 0) {
            syntax = FrameSyntax::D;
            break;
          }
        }
        argi += 2;
      }
      else if ((strcmp(argv[argi], "-gauss_points") == 0) ||
               (strcmp(argv[argi], "-n") == 0) || 
               (strcmp(argv[argi], "-gauss_type") == 0)) {
        syntax = FrameSyntax::X;
        argi += 2;
      }
      else if ((strcmp(argv[argi], "-gauss") == 0)) {
        if (argc < argi + 2) {
          opserr << OpenSees::PromptValueError 
                 << "not enough arguments, expected -gauss $type\n";
          return TCL_ERROR;
        }
        int itag;
        if (Tcl_GetInt(interp, argv[argi + 1], &itag) != TCL_OK)
          syntax = FrameSyntax::X;
        else
          syntax = FrameSyntax::A;
        argi += 2;
      }
      // Section
      else if (strcmp(argv[argi], "-section") == 0) {
        syntax = FrameSyntax::X;
        argi += 2;
      }

      // Positional argument
      else {
        positions.push_back(argi);
        argi++;
      }
    }
  }

  //
  // II Parse Positional Arguments
  //
  // If we get a BeamIntegration from a BeamIntegrationRule
  // then we dont own it and can't delete it
  bool deleteBeamIntegr = true;
  BeamIntegration   *beamIntegr   = nullptr;
  std::vector<int>   section_tags;
  if (syntax == FrameSyntax::None) {
    // a or b
    if (positions.size() == 2 || positions.size() > 3) {
      int itg_tag;
      if (Tcl_GetInt(interp, argv[positions[1]], &itg_tag) == TCL_OK)
        syntax = FrameSyntax::B;
      else
        syntax = FrameSyntax::A;
    }
    // c
    else if (positions.size() == 3)
      syntax = FrameSyntax::C;
  }

  // Version a or b
  if (syntax == FrameSyntax::A ||
      syntax == FrameSyntax::B) {
    // Here we create a BeamIntegrationRule (theRule) which is a pair of
    // section tags and a BeamIntegration. In this case we do not
    // delete the BeamIntegration because it is owned by theRule.
    deleteBeamIntegr = false;


    // Geometric transformation
    if (Tcl_GetInt(interp, argv[positions[0]], &transfTag) != TCL_OK) {
      opserr << OpenSees::PromptValueError 
             << "invalid transform " << argv[positions[0]] 
             << OpenSees::SignalMessageEnd;
      status = TCL_ERROR;
      goto clean_up;
    }
    received.transform = true;

    int itg_tag;
    // Version b)
    if (syntax == FrameSyntax::B) {
      if (Tcl_GetInt(interp, argv[positions[1]], &itg_tag) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "invalid integration tag " << argv[positions[1]] 
               << OpenSees::SignalMessageEnd;
        status = TCL_ERROR;
        goto clean_up;
      }
    }
    // Version a)
    else {
#if !defined(OPS_API)
      // If we fail to parse an integer tag for the integration,
      // then we assume that the integration is specified as an
      // inline BeamIntegration command
      builder->findFreeTag<BeamIntegrationRule>(itg_tag);
      std::string integrCommand{argv[positions[1]]};
      if (integrCommand.find(" ") == std::string::npos) {
        for (std::vector<int>::size_type i =2; i< positions.size(); i++) {
          integrCommand += " " + std::string(argv[positions[i]]);
        }
      }
      integrCommand.insert(integrCommand.find(" "), " "+std::to_string(itg_tag)+" ");
      integrCommand.insert(0, "beamIntegration ");
      if (Tcl_Eval(interp, integrCommand.c_str()) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to parse integration\n";
        return TCL_ERROR;
      }
#else
      return TCL_ERROR;
#endif
    }

    BeamIntegrationRule* theRule =
        builder->getTypedObject<BeamIntegrationRule>(itg_tag);
    if (theRule == nullptr)
      return TCL_ERROR;

    beamIntegr = theRule->getBeamIntegration()->getCopy();
    const ID& secTags = theRule->getSectionTags();
    received.integration = true;

    for (int i=0; i < secTags.Size(); i++)
      section_tags.push_back(secTags(i));
    received.sections = true;

    if (syntax == FrameSyntax::A) {
      builder->removeObject<BeamIntegrationRule>(itg_tag);
      delete theRule;
    }
  }

  // Version c)
  //
  // .. $nip $section $transf
  else if (syntax == FrameSyntax::C) {
    deleteBeamIntegr = true;

    int nIP;
    if (Tcl_GetInt(interp, argv[positions[0]], &nIP) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid nIP\n";
      status = TCL_ERROR;
      goto clean_up;
    }
    if (nIP <= 0) {
      opserr << OpenSees::PromptValueError << "invalid nIP, must be positive.\n";
      status = TCL_ERROR;
      goto clean_up;
    }

    //
    int secTag;
    if (Tcl_GetInt(interp, argv[positions[1]], &secTag) != TCL_OK) {
      opserr << OpenSees::PromptValueError 
             << "invalid secTag\n";
      status = TCL_ERROR;
      goto clean_up;
    }

    for (int i=0; i < nIP; i++)
      section_tags.push_back(secTag);
    received.sections = true;

    if ((beamIntegr = GetBeamIntegration("Lobatto", nIP)) == nullptr) {
      opserr << OpenSees::PromptValueError 
             << "invalid integration type or size\n";
      status = TCL_ERROR;
      goto clean_up;
    }
    received.integration = true;

    // Transform
    if (Tcl_GetInt(interp, argv[positions[2]], &transfTag) != TCL_OK) {
      opserr << OpenSees::PromptValueError 
             << "invalid transform"
             << OpenSees::SignalMessageEnd;
      status = TCL_ERROR;
      goto clean_up;
    }
    received.transform = true;
  }
  // Version d)
  else if (syntax == FrameSyntax::D) {
    // positional arguments are:
    //   0: nIP
    //   1: -sections
    //   2: secTag1
    //   3: secTag2...
    deleteBeamIntegr = true;
  
    int nIP;
    if (Tcl_GetInt(interp, argv[positions[0]], &nIP) != TCL_OK) {
      opserr << OpenSees::PromptValueError 
             << "invalid nIP\n";
      return TCL_ERROR;
    }
    if (nIP <= 0) {
      opserr << OpenSees::PromptValueError 
             << "invalid nIP, must be positive."
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }
    // TODO: Make sure 2+nIP < positions.size()

    // Get section tags
    for (int i = 0; i < nIP; i++) {
      int secTag;
      if (Tcl_GetInt(interp, argv[positions[2+i]], &secTag) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "invalid section"
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      section_tags.push_back(secTag);
    }
    if (section_tags.size() != static_cast<size_t>(nIP)) {
      opserr << OpenSees::PromptValueError 
             << "number of sections does not match nIP"
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }
    received.sections = true;

    if (Tcl_GetInt(interp, argv[positions[2+nIP]], &transfTag) != TCL_OK) {
      opserr << OpenSees::PromptValueError 
             << "invalid transform"
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }
    received.transform = true;

    for (int ii=0; ii<argc; ii++) {
      if (strcmp(argv[ii], "-integration") == 0) {
        if (argc < ii + 2) {
          opserr << OpenSees::PromptValueError 
                 << "not enough arguments, expected -integration <type>"
                 << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        }
        beamIntegr = GetBeamIntegration(argv[ii+1], nIP);
        if (beamIntegr == nullptr) {
          opserr << OpenSees::PromptValueError 
                 << "invalid integration type"
                 << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        }
        received.integration = true;
        break;
      }
    }
  }
  else if (syntax == FrameSyntax::X) {
    deleteBeamIntegr = true;
    int n = 5;
    int section_tag = -1;
    const char* gauss_type = "Legendre";

    for (int ii=0; ii<argc; ii++) {
      if ((strcmp(argv[ii], "-gauss") == 0) || 
          (strcmp(argv[ii], "-gauss_type") == 0)) {
        if (argc < ii + 2) {
          opserr << OpenSees::PromptValueError 
                 << "not enough arguments, expected -gauss <type>"
                 << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        }
        gauss_type = argv[ii + 1];
      }
      else if ((strcmp(argv[ii], "-n") == 0) ||
               (strcmp(argv[ii], "-gauss_points") == 0)) {
        if (argc < ii + 2) {
          opserr << OpenSees::PromptValueError 
                 << "not enough arguments, expected -n <n>"
                 << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        }
        if (Tcl_GetInt(interp, argv[ii + 1], &n) != TCL_OK) {
          opserr << OpenSees::PromptValueError 
                 << "invalid n: " << argv[ii + 1]
                 << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        }
      }
      else if (strcmp(argv[ii], "-section") == 0) {
        if (argc < ii + 2) {
          opserr << OpenSees::PromptValueError 
                 << "not enough arguments, expected -section <section>"
                 << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        }

        if (Tcl_GetInt(interp, argv[ii + 1], &section_tag) != TCL_OK) {
          opserr << OpenSees::PromptValueError 
                 << "invalid section tag " << argv[ii + 1]
                 << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        }
      }
    }
    if (section_tag == -1) {
      opserr << OpenSees::PromptValueError 
             << "missing required argument: -section <section>"
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }
    for (int ii=0; ii<n; ii++)
      section_tags.push_back(section_tag);
    received.sections = true;

    if (gauss_type == nullptr) {
      if (strcasecmp(argv[1], "ForceFrame") == 0) {
        // integration_type = "Legendre";
        gauss_type = "Lobatto";
      } else if (strstr(argv[1], "ispBeam") == 0) {
        // forceBeamColumn
        gauss_type = "Lobatto";
      } else {
        gauss_type = "Legendre";
      }
    }
    if ((beamIntegr = GetBeamIntegration(gauss_type, n)) == nullptr) {
      opserr << OpenSees::PromptValueError 
             << "invalid integration type or size\n";
      return TCL_ERROR;
    }
    received.integration = true;
  }
  else if (syntax == FrameSyntax::None) {
    opserr << OpenSees::PromptValueError 
           << "could not determine syntax for frame element\n";
    return TCL_ERROR;
  }

  //
  // Finalize options
  //
  if (!received.transform) {
    opserr << OpenSees::PromptValueError 
           << "missing required argument: transform\n";
    status = TCL_ERROR;
    goto clean_up;
  }
  if (!received.integration) {
    opserr << OpenSees::PromptValueError 
           << "missing required argument: integration\n";
    status = TCL_ERROR;
    goto clean_up;
  }
  if (!received.sections) {
    opserr << OpenSees::PromptValueError 
           << "missing required argument: sections\n";
    status = TCL_ERROR;
    goto clean_up;
  }
  assert(beamIntegr != nullptr);

  //
  //
  {
    Element *theElement = ndm == 2 
                        ? CreateFrame<2, CrdTransf, FrameSection>(*builder, argv[1], tag, multi_nodes, transfTag, 
                                                                  section_tags, *beamIntegr, mass, max_iter, tol, 
                                                                  shear_center, options)
                        : CreateFrame<3, CrdTransf, FrameSection>(*builder, argv[1], tag, multi_nodes, transfTag, 
                                                                  section_tags, *beamIntegr, mass, max_iter, tol, 
                                                                  shear_center, options);

                                                                        
    if (theElement == nullptr) {
      status = TCL_ERROR;
      goto clean_up;
    }
    if (domain->addElement(theElement) == false) {
      opserr << OpenSees::PromptValueError 
             << "could not add element to the domain"
             << OpenSees::SignalMessageEnd;
      delete theElement;
      status = TCL_ERROR;
      goto clean_up;
    }
  }


clean_up:
  //
  // Clean up
  //

  if (deleteBeamIntegr && beamIntegr != nullptr)
    delete beamIntegr;

  return status;
}



//
//  BeamWithHinges
//
//     element beamWithHinges tag? ndI? ndJ? secTagI? lenI? secTagJ? lenJ? 
//        E? A? I? transfTag? <-shear shearLength?> <-mass massDens?> 
//        <-iter maxIters tolerance>
//
int
TclBasicBuilder_addBeamWithHinges(ClientData clientData, Tcl_Interp *interp,
                                  int argc, TCL_Char ** const argv)
{
  ModelRegistry *builder = (ModelRegistry*)clientData;

  int NDM = builder->getNDM();
  int NDF = builder->getNDF();

  // Plane frame element
  if (NDM == 2 && NDF == 3) {
    if (argc < 13) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: element beamWithHinges tag? ndI? ndJ? secTagI? lenI? "
                "secTagJ? lenJ? ";
      opserr << "E? A? I? transfTag? <-shear shearLength?> <-mass massDens?> "
                "<-iter maxIters tolerance>"
             << "\n";
      return TCL_ERROR;
    }

    double massDens = 0.0;
    int max_iters = 10;
    double tol = 1.0e-10;
    int tag, ndI, ndJ, secTagI, secTagJ, transfTag;
    double lenI, lenJ, E, A, I;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid beamWithHinges tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &ndI) != TCL_OK) {
      opserr << "WARNING invalid ndI\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[4], &ndJ) != TCL_OK) {
      opserr << "WARNING invalid ndJ\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[5], &secTagI) != TCL_OK) {
      opserr << "WARNING invalid secTagI\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &lenI) != TCL_OK) {
      opserr << "WARNING invalid lenI\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[7], &secTagJ) != TCL_OK) {
      opserr << "WARNING invalid ndJ\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[8], &lenJ) != TCL_OK) {
      opserr << "WARNING invalid lenJ\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[9], &E) != TCL_OK) {
      opserr << "WARNING invalid E\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[10], &A) != TCL_OK) {
      opserr << "WARNING invalid A\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[11], &I) != TCL_OK) {
      opserr << "WARNING invalid I\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[12], &transfTag) != TCL_OK) {
      opserr << "WARNING invalid transfTag\n";
      return TCL_ERROR;
    }

    bool isShear = false;
    int shearTag = 0;

    if (argc > 13) {
      for (int i = 13; i < argc; ++i) {
        if (strcmp(argv[i], "-mass") == 0 && ++i < argc) {
          if (Tcl_GetDouble(interp, argv[i], &massDens) != TCL_OK) {
            opserr << "WARNING invalid massDens\n";
            opserr << "BeamWithHinges: " << tag << "\n";
            return TCL_ERROR;
          }
        }

        if (strcmp(argv[i], "-constHinge") == 0 && ++i < argc) {
          if (Tcl_GetInt(interp, argv[i], &shearTag) != TCL_OK) {
            opserr << "WARNING invalid constHinge tag\n";
            opserr << "BeamWithHinges: " << tag << "\n";
            return TCL_ERROR;
          }
          isShear = true;
        }

        if (strcmp(argv[i], "-iter") == 0 && i + 2 < argc) {
          if (Tcl_GetInt(interp, argv[++i], &max_iters) != TCL_OK) {
            opserr << "WARNING invalid maxIters\n";
            opserr << "BeamWithHinges: " << tag << "\n";
            return TCL_ERROR;
          }
          if (Tcl_GetDouble(interp, argv[++i], &tol) != TCL_OK) {
            opserr << "WARNING invalid tolerance\n";
            opserr << "BeamWithHinges: " << tag << "\n";
            return TCL_ERROR;
          }
        }
      }
    }

    // Retrieve section I from the model builder
    FrameSection *sectionI = builder->getTypedObject<FrameSection>(secTagI);
    if (sectionI == nullptr)
      return TCL_ERROR;

    // Retrieve section J from the model builder
    FrameSection *sectionJ = builder->getTypedObject<FrameSection>(secTagJ);
    if (sectionJ == nullptr)
      return TCL_ERROR;


    CrdTransf *theTransf = builder->getTypedObject<CrdTransf>(transfTag);
    if (theTransf == nullptr)
      return TCL_ERROR;

    Element *theElement = nullptr;
    int numSections = 0;
    SectionForceDeformation *sections[10];
    BeamIntegration *theBeamIntegr = nullptr;

    ElasticSection2d elastic(8, E, A, I);

    if (strcmp(argv[1], "beamWithHinges1") == 0) {
      theBeamIntegr = new HingeMidpointBeamIntegration(lenI, lenJ);

      numSections = 4;

      sections[0] = sectionI;
      sections[1] = &elastic;
      sections[2] = &elastic;
      sections[3] = sectionJ;

    } else if (strcmp(argv[1], "beamWithHinges2") == 0) {
      theBeamIntegr = new HingeRadauTwoBeamIntegration(lenI, lenJ);

      numSections = 6;
      sections[0] = sectionI;
      sections[1] = sectionI;
      sections[2] = &elastic;
      sections[3] = &elastic;
      sections[4] = sectionJ;
      sections[5] = sectionJ;

    } else if (strcmp(argv[1], "beamWithHinges3") == 0 ||
               strcmp(argv[1], "beamWithHinges") == 0) {
      theBeamIntegr = new HingeRadauBeamIntegration(lenI, lenJ);

      numSections = 6;
      sections[0] = sectionI;
      sections[1] = &elastic;
      sections[2] = &elastic;
      sections[3] = &elastic;
      sections[4] = &elastic;
      sections[5] = sectionJ;

    } else if (strcmp(argv[1], "beamWithHinges4") == 0) {
      theBeamIntegr = new HingeEndpointBeamIntegration(lenI, lenJ);

      numSections = 4;
      sections[0] = sectionI;
      sections[1] = &elastic;
      sections[2] = &elastic;
      sections[3] = sectionJ;
    }

    if (theBeamIntegr == nullptr) {
      opserr << "Unknown element type: " << argv[1] << "\n";
      return TCL_ERROR;
    }

    if (isShear) {
      FrameSection *sectionL = builder->getTypedObject<FrameSection>(shearTag);
      if (sectionL == nullptr)
        return TCL_ERROR;

      sections[numSections++] = sectionL;
    }

    theElement = new ForceBeamColumn2d(tag, ndI, ndJ, numSections, 
                                       sections,
                                       *theBeamIntegr, *theTransf, massDens,
                                       max_iters, tol);

    delete theBeamIntegr;

    if (builder->getDomain()->addElement(theElement) == false) {
      opserr << "WARNING could not add element to domain.\n";
      return TCL_ERROR;
    }
  }

  else if (NDM == 3 && NDF == 6) {
    if (argc < 16) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: element beamWithHinges tag? ndI? ndJ? secTagI? lenI? "
                "secTagJ? lenJ? ";
      opserr << "E? A? Iz? Iy? G? J? transfTag? <-shear shearLength?> <-mass "
                "massDens?> <-iter maxIters tolerance>"
             << "\n";
      return TCL_ERROR;
    }

    int tag, ndI, ndJ, secTagI, secTagJ, transfTag;
    double lenI, lenJ, E, A, Iz, Iy, G, J;
    double massDens = 0.0;
    int max_iters = 10;
    double tol = 1.0e-10;
    double shearLength = 1.0;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid beamWithHinges tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &ndI) != TCL_OK) {
      opserr << "WARNING invalid ndI\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[4], &ndJ) != TCL_OK) {
      opserr << "WARNING invalid ndJ\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[5], &secTagI) != TCL_OK) {
      opserr << "WARNING invalid secTagI\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &lenI) != TCL_OK) {
      opserr << "WARNING invalid lenI\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[7], &secTagJ) != TCL_OK) {
      opserr << "WARNING invalid ndJ\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[8], &lenJ) != TCL_OK) {
      opserr << "WARNING invalid lenJ\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[9], &E) != TCL_OK) {
      opserr << "WARNING invalid E\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[10], &A) != TCL_OK) {
      opserr << "WARNING invalid A\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[11], &Iz) != TCL_OK) {
      opserr << "WARNING invalid Iz\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[12], &Iy) != TCL_OK) {
      opserr << "WARNING invalid Iy\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[13], &G) != TCL_OK) {
      opserr << "WARNING invalid G\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[14], &J) != TCL_OK) {
      opserr << "WARNING invalid J\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[15], &transfTag) != TCL_OK) {
      opserr << "WARNING invalid transfTag\n";
      return TCL_ERROR;
    }


    if (argc > 16) {
      for (int i = 16; i < argc; ++i) {
        if (strcmp(argv[i], "-mass") == 0 && ++i < argc) {
          if (Tcl_GetDouble(interp, argv[i], &massDens) != TCL_OK) {
            opserr << "WARNING invalid massDens\n";
            opserr << "BeamWithHinges: " << tag << "\n";
            return TCL_ERROR;
          }
        }

        if (strcmp(argv[i], "-shear") == 0 && ++i < argc) {
          if (Tcl_GetDouble(interp, argv[i], &shearLength) != TCL_OK) {
            opserr << "WARNING invalid shear\n";
            return TCL_ERROR;
          }
        }

        if (strcmp(argv[i], "-iter") == 0 && i + 2 < argc) {
          if (Tcl_GetInt(interp, argv[++i], &max_iters) != TCL_OK) {
            opserr << "WARNING invalid maxIters\n";
            return TCL_ERROR;
          }
          if (Tcl_GetDouble(interp, argv[++i], &tol) != TCL_OK) {
            opserr << "WARNING invalid tolerance\n";
            return TCL_ERROR;
          }
        }
      }
    }

    // Retrieve section I from the model builder
    SectionForceDeformation *sectionI = builder->getTypedObject<FrameSection>(secTagI);
    if (sectionI == nullptr)
      return TCL_ERROR;

    // Retrieve section J from the model builder
    SectionForceDeformation *sectionJ = builder->getTypedObject<FrameSection>(secTagJ);
    if (sectionJ == nullptr)
      return TCL_ERROR;


    CrdTransf *theTransf = builder->getTypedObject<CrdTransf>(transfTag);
    if (theTransf == nullptr)
      return TCL_ERROR;


    Element *theElement = nullptr;
    int numSections = 0;
    SectionForceDeformation *sections[10];
    BeamIntegration *theBeamIntegr = nullptr;

    ElasticSection3d elastic(0, E, A, Iz, Iy, G, J);

    if (strcmp(argv[1], "beamWithHinges1") == 0) {
      theBeamIntegr = new HingeMidpointBeamIntegration(lenI, lenJ);

      numSections = 4;
      sections[0] = sectionI;
      sections[1] = &elastic;
      sections[2] = &elastic;
      sections[3] = sectionJ;

    } else if (strcmp(argv[1], "beamWithHinges2") == 0) {
      theBeamIntegr = new HingeRadauTwoBeamIntegration(lenI, lenJ);

      numSections = 6;
      sections[0] = sectionI;
      sections[1] = sectionI;
      sections[2] = &elastic;
      sections[3] = &elastic;
      sections[4] = sectionJ;
      sections[5] = sectionJ;

    } else if (strcmp(argv[1], "beamWithHinges3") == 0 ||
               strcmp(argv[1], "beamWithHinges") == 0) {
      theBeamIntegr = new HingeRadauBeamIntegration(lenI, lenJ);

      numSections = 6;
      sections[0] = sectionI;
      sections[1] = &elastic;
      sections[2] = &elastic;
      sections[3] = &elastic;
      sections[4] = &elastic;
      sections[5] = sectionJ;

    } else if (strcmp(argv[1], "beamWithHinges4") == 0) {
      theBeamIntegr = new HingeEndpointBeamIntegration(lenI, lenJ);

      numSections = 4;
      sections[0] = sectionI;
      sections[1] = &elastic;
      sections[2] = &elastic;
      sections[3] = sectionJ;
    }

    if (theBeamIntegr == nullptr) {
      opserr << "Unknown element type: " << argv[1] << "\n";
      return TCL_ERROR;
    }

    // TODO fix shear for beamWithHinges
    /*
    if (isShear) {
      SectionForceDeformation *sectionL = builder->getTypedObject<SectionForceDeformation>(shearTag);

      if (sectionL == 0) {
        opserr << "WARNING section L does not exist\n";
        opserr << "section: " << shearTag;
        opserr << "\nBeamWithHinges: " << tag << "\n";
        return TCL_ERROR;
      }
      sections[numSections++] = sectionL;
    }
    */

    theElement = new ForceBeamColumn3d(tag, ndI, ndJ, numSections, sections,
                                       *theBeamIntegr, *theTransf, massDens,
                                       max_iters, tol);

    delete theBeamIntegr;

    // Add to the domain
    if (builder->getDomain()->addElement(theElement) == false) {
      opserr << "WARNING could not add "
                "element to domain ";
      opserr << tag << "\n";
      return TCL_ERROR;
    }
  }

  else {
    opserr << "ERROR -- model dimension: " << NDM
           << " and nodal degrees of freedom: " << NDF
           << " are incompatible for BeamWithHinges element" << "\n";
    return TCL_ERROR;
  }

  return TCL_OK;
}
