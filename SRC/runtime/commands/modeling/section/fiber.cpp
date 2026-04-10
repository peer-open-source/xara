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
#include <ModelRegistry.h>
#include <Parsing.h>
#include <Logging.h>
#include <ArgumentTracker.h>

#include <QuadFiberPatch.h>
#include <CircPatch.h>
#include <StraightFiberLayer.h>
#include <CircFiberLayer.h>

#include <FiberSection2d.h>
#include <NDFiberSection2d.h>
#include <NDFiberSection3d.h>
#include <NDFiberSectionWarping2d.h>
#include <FiberSection2dInt.h>
#include <FiberSection3d.h>
#include <MixedFrameSection.h>
#include <FrameFiberSection3d.h>
#include <FrameSolidSection3d.h>
#include <FrameTraceSection3d.h>
#include <FiberSectionAsym3d.h>
#include <FiberSection2dThermal.h>
#include <FiberSection3dThermal.h>

#include <NDMaterial.h>
#include <UniaxialMaterial.h>
#include <ElasticMaterial.h>
#include <FiberSectionBuilder.h>

struct FiberSectionConfig {
  bool isND            = false;
  bool isAsym          = false;
  bool isWarping       = false;
  bool isThermal       = false;
  bool isMixed         = false;
  bool isNew           = false; // use new FrameFiberSection class
  bool computeCentroid = true;
  bool use_twist = false;
  bool use_density = false;
  int  reserve = 30;
  bool wagner = false;
};

SectionBuilder* 
findSectionBuilder(ModelRegistry* builder, Tcl_Interp *interp, int argc, const char** const argv)
{
  int tag;
  bool section_passed = false;
  for (int i = 0; i<argc; ++i) {
    if (strcmp(argv[i], "-section") == 0) {
      if (Tcl_GetInt(interp, argv[i+1], &tag) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to parse section tag \"" << argv[i+1] << "\"\n";
        return nullptr;
      } else {
        section_passed = true;
        break;
      }
    }
  }

  if (!section_passed)
   if (builder->getCurrentSectionBuilder(tag) != 0) {
     return nullptr;
   }

  if (tag == -1)
    return nullptr;

  return builder->getTypedObject<SectionBuilder>(tag);
}




// build the section
// This function assumes torsion is not NULL when num==3
static int
initSectionCommands(ClientData clientData,
                    Tcl_Interp *interp,
                    int secTag, 
                    UniaxialMaterial *theTorsion, 
                    const Frame::Shape& shape_data,
                    const FiberSectionConfig& options)
{
  assert(clientData != nullptr);
  ModelRegistry *builder = static_cast<ModelRegistry*>(clientData);

  // Dimension of the structure
  int ndm = builder->getNDM();

  SectionBuilder  *sbuilder = nullptr;
  FrameSection    *section  = nullptr;
  // Create 2d section
  if (ndm == 2) {
    if (options.isND) {
      if (options.isNew) {
        auto sec = new FrameSolidSection3d(secTag, options.reserve);
        sbuilder = new FiberSectionBuilder<2, NDMaterial, FrameSolidSection3d>(*builder, *sec);
        section = sec;
      }
      else if (options.isWarping) {
        auto sec = new NDFiberSectionWarping2d(secTag, options.reserve, shape_data.alpha);
        sbuilder = new FiberSectionBuilder<2, NDMaterial, NDFiberSectionWarping2d>(*builder, *sec);
        section = sec;
      }
      else {
        auto sec = new NDFiberSection2d(secTag, options.reserve, shape_data.alpha, options.computeCentroid);
        sbuilder = new FiberSectionBuilder<2, NDMaterial, NDFiberSection2d>(*builder, *sec);
        section = sec;
      }
    } else {
      if (options.isThermal) {
        auto sec = new FiberSection2dThermal(secTag, options.computeCentroid);
        sbuilder = new FiberSectionBuilder<2, UniaxialMaterial, FiberSection2dThermal>(*builder, *sec);
        section = sec;
      } else {
        auto sec = new FiberSection2d(secTag, options.reserve, options.computeCentroid);
        sbuilder = new FiberSectionBuilder<2, UniaxialMaterial, FiberSection2d>(*builder, *sec);
        section = sec;
      }
    }
  }

  else if (ndm == 3) {

    if (options.isND) {
      if (options.isMixed) {
        auto sec = new MixedFrameSection(secTag, options.reserve, shape_data.mixed_form, options.wagner);
        sbuilder = new FiberSectionBuilder<3, NDMaterial, MixedFrameSection>(*builder, *sec);
        section = sec;
      }
      else if (options.isNew) {
        if (!getenv("XARA_OLD_WARP")) {
          auto sec = new MixedFrameSection(secTag, options.reserve, shape_data.mixed_form, options.wagner);
          sbuilder = new FiberSectionBuilder<3, NDMaterial, MixedFrameSection>(*builder, *sec);
          section = sec;
        } else {
          auto sec = new FrameTraceSection3d(secTag, options.reserve, options.wagner);
          sbuilder = new FiberSectionBuilder<3, NDMaterial, FrameTraceSection3d>(*builder, *sec);
          section = sec;
          // auto sec = new FrameSolidSection3d(secTag, options.reserve);
          // sbuilder = new FiberSectionBuilder<3, NDMaterial, FrameSolidSection3d>(*builder, *sec);
          // section = sec;
        }
      }
      else {
        auto sec = new NDFiberSection3d(secTag,
                                        options.reserve,
                                        shape_data.alpha,
                                        options.computeCentroid);
        sbuilder = new FiberSectionBuilder<3, NDMaterial, NDFiberSection3d>(*builder, *sec);
        section = sec;
      }
    } 
    else {

      if (options.isThermal) {
        if (theTorsion == nullptr) {
          opserr << OpenSees::PromptValueError 
                 << "FiberThermal section requires torsion"
                 << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        }
        auto sec = new FiberSection3dThermal(secTag, 30, *theTorsion,
                                             options.computeCentroid);
        sbuilder = new FiberSectionBuilder<3, UniaxialMaterial, FiberSection3dThermal>(*builder, *sec);
        section = sec;

      }
      else if (options.isAsym) {
        auto sec = new FiberSectionAsym3d(secTag, 30, theTorsion, 
                                          shape_data.mixed.shear_origin[0], 
                                          shape_data.mixed.shear_origin[1]);
        sbuilder = new FiberSectionBuilder<3, UniaxialMaterial, FiberSectionAsym3d>(*builder, *sec);
        section = sec;
      }
      else {
        if (theTorsion == nullptr) {
          opserr << OpenSees::PromptValueError 
                << "Fiber section requires torsion in 3D\n";
          return TCL_ERROR;
        }
        if (options.isNew) {
          auto sec = new FrameFiberSection3d(secTag, 
                                             30,
                                             shape_data,
                                             *shape_data.density, 
                                             options.use_density);
          sbuilder = new FiberSectionBuilder<3, UniaxialMaterial, FrameFiberSection3d>(*builder, *sec);
          section = sec;
        }
        else {
          auto sec = new FiberSection3d(secTag, 30, *theTorsion, options.computeCentroid);
          sbuilder = new FiberSectionBuilder<3, UniaxialMaterial, FiberSection3d>(*builder, *sec);
          section = sec;
        }
      }
    }

  } else {
    opserr << OpenSees::PromptValueError << "Model dimension (ndm = " << ndm
           << ") is incompatible with available frame elements\n";
    return TCL_ERROR;
  }

  //
  // Set extra data
  //
  if (shape_data.mixed.shear_align) {
    for (int i=0; i<2; ++i){
      for (int j=0; j<2; ++j) {
        Parameter p(-1, nullptr, nullptr, 0);
        const char* argv_ij[] {
          "shift_shear",
          i? "2" : "1",
          j? "2" : "1",
        };
        if (section->setParameter(argv_ij, 3, p) < 0) {
          opserr << OpenSees::PromptValueError 
                 << "Failed to set shear alignment parameter "
                 << argv_ij[1] << argv_ij[2] << "\n";
          return TCL_ERROR;
        }
        p.update((*shape_data.mixed.shear_align)(i,j));
      }
    }
  }
  if (shape_data.mixed.shift_twist) {
    for (int i=0; i<2; ++i){
      Parameter p(-1, nullptr, nullptr, 0);
      const char* argv_i[] {
        "shift_twist",
        i? "2" : "1",
      };
      if (section->setParameter(argv_i, 2, p) < 0) {
        opserr << OpenSees::PromptValueError 
               << "Failed to set twist shift parameter "
               << argv_i[1] << "\n";
        return TCL_ERROR;
      }
      p.update((*shape_data.mixed.shift_twist)(i));
    }
  }
  if (shape_data.mixed.shift_axial) {
    for (int i=0; i<2; ++i){
      Parameter p(-1, nullptr, nullptr, 0);
      const char* argv_i[] {
        "shift_axial",
        i? "2" : "1",
      };
      if (section->setParameter(argv_i, 2, p) < 0) {
        opserr << OpenSees::PromptValueError 
               << "Failed to set axial shift parameter "
               << argv_i[1] << "\n";
        return TCL_ERROR;
      }
      p.update((*shape_data.mixed.shift_axial)(i));
    }
  }


  // In 2D truss elements still look for FrameSections
  if (builder->addTaggedObject<FrameSection>(*section) < 0) {
    return TCL_ERROR;
  }
  // if (ndm == 2) {
  //   if (builder->addTaggedObject<SectionForceDeformation>(*section->getCopy()) < 0) {
  //     return TCL_ERROR;
  //   }
  // }

  if (builder->addTypedObject<SectionBuilder>(secTag, sbuilder) < 0) {
    opserr << OpenSees::PromptValueError << "Faled to add section\n";
    return TCL_ERROR;
  }

  return TCL_OK;
}

int
TclCommand_addFiberSection(ClientData clientData, Tcl_Interp *interp, int argc,
                           TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  ModelRegistry* builder = static_cast<ModelRegistry*>(clientData);

  // Check if we are being invoked from Python or Tcl
  bool openseespy = false;
  if (Tcl_GetVar(interp, "opensees::pragma::openseespy", 0) != nullptr) { 
    openseespy = true;
  }

  int ndm = builder->getNDM();

  // cmp - Check argument counts; In Tcl we require the brace argument
  //       in Python it can be omitted, but only for legacy reasons
  if (argc < 4 && !openseespy) {
    opserr << "Insufficient arguments, expected at least 4\n";
    return TCL_ERROR;
  }
  else if (argc < 3 && openseespy) {
    opserr << "Insufficient arguments, expected at least 3\n";
    return TCL_ERROR;
  }


  int secTag;
  if (Tcl_GetInt(interp, argv[2], &secTag) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "failed to parse section tag \"" << argv[2] << "\""
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  builder->setCurrentSectionBuilder(secTag);

  FiberSectionConfig options;
  if (strcmp(argv[1], "NDFiber") == 0 ||
      strcmp(argv[1], "ShearFiber") == 0 ||
      strcmp(argv[1], "MixedFiber") == 0)
    options.isND = true;
  
  if (strcmp(argv[1], "MixedFiber") == 0)
    options.isMixed = true;

  if (strcmp(argv[1], "NDFiberWarping") == 0) {
    options.isND = true;
    options.isWarping = true;
  }
  else if (strcmp(argv[1], "FrameFiber") == 0 ||
           strcmp(argv[1], "FiberFrame") == 0 ||
           strcmp(argv[1], "AxialFiber") == 0 ||
           strcmp(argv[1], "ShearFiber") == 0 ||
           strcmp(argv[1], "MixedFiber") == 0
    )
    options.isNew = true;

  else if (strcmp(argv[1], "FiberThermal") == 0 ||
            strcmp(argv[1], "fiberSecThermal") == 0)
    options.isThermal = true;

  else if (strstr(argv[1], "Asym") != nullptr)
    options.isAsym    = true;


  int iarg  = 3;
  // FiberSectionData data;
  Frame::Shape shape_data(builder->getNDM(), builder->getNDF());

  // bool shape_done = false;

  if (builder->getNDF() <= 6)
    shape_data.mixed_form = MixedFrameSection::MixedType::UT;//Energetic;
  else 
    shape_data.mixed_form = MixedFrameSection::MixedType::None;

  UniaxialMaterial *torsion = nullptr;
  bool deleteTorsion = false;
  bool shearParsed = false;

//// Interaction parameters
// int NStrip1, NStrip2, NStrip3;
// double t1, t2, t3;

  while (iarg < argc) {

    if (strcmp(argv[iarg], "-noCentroid") == 0) {
      options.computeCentroid = false;
      iarg += 1;
    }

    else if (strcmp(argv[iarg], "-reserve") == 0 && iarg + 1 < argc) {
      int reserve;
      if (Tcl_GetInt(interp, argv[iarg + 1], &reserve) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid reserve" << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      options.reserve = reserve;
      iarg += 2;
    }

    else if (strcmp(argv[iarg], "-wagner") == 0) {
      options.wagner = true;
      iarg += 1;
    }

    else if (strcmp(argv[iarg], "-mass") == 0 && iarg + 1 < argc) {
      if (argc < iarg + 2) {
        opserr << OpenSees::PromptValueError << "not enough -mass args need -mass mass?\n";
        return TCL_ERROR;
      }
      double density;
      if (Tcl_GetDouble(interp, argv[iarg + 1], &density) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid density";
        return TCL_ERROR;
      }
      shape_data.density = density;
      options.use_density = true;

      iarg  += 2;
    }

    else if (strcmp(argv[iarg], "-GJ") == 0 && iarg + 1 < argc) {
      double G;
      if (Tcl_GetDouble(interp, argv[iarg + 1], &G) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid GJ";
        return TCL_ERROR;
      }
      deleteTorsion = true;
      torsion = new ElasticMaterial(0, G);

      iarg  += 2;
    }

    else if (strcmp(argv[iarg], "-alpha") == 0) {
      // -alpha <value>
      if (iarg + 1 >= argc) {
        opserr << OpenSees::PromptValueError 
               << "alpha value missing after -alpha"
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      double alpha;
      if (Tcl_GetDouble(interp, argv[iarg + 1], &alpha) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "invalid alpha value after -alpha"
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      shape_data.alpha = alpha;
      iarg += 2;
    }

    else if (strcmp(argv[iarg], "-mixed_type") == 0) {
      // -mixed_type {None|constant|energetic}
      if (iarg + 1 >= argc) {
        opserr << OpenSees::PromptValueError 
               << "mixed type value missing after -mixed_type"
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      if ((strcmp(argv[iarg + 1], "None") == 0) || 
          (strcmp(argv[iarg + 1], "NV") == 0) ||
          (strcmp(argv[iarg + 1], "NT") == 0)) {
        shape_data.mixed_form = MixedFrameSection::MixedType::None;
      } else if ((strcmp(argv[iarg + 1], "constant") == 0) || 
                 (strcmp(argv[iarg + 1], "geometric") == 0) ||
                 (strcmp(argv[iarg + 1], "UG") == 0)) {
        shape_data.mixed_form = MixedFrameSection::MixedType::Constant;
      } else if (strcmp(argv[iarg + 1],  "UT") == 0) {
        shape_data.mixed_form = MixedFrameSection::MixedType::UT;
      } else if (strcmp(argv[iarg + 1],  "U02") == 0) {
        shape_data.mixed_form = MixedFrameSection::MixedType::U02;
      } else if ((strcmp(argv[iarg + 1], "energetic") == 0) || 
                 (strcmp(argv[iarg + 1], "UE") == 0) ) {
        shape_data.mixed_form = MixedFrameSection::MixedType::Energetic;
      } else if (strcmp(argv[iarg + 1],  "NR") == 0) {
        shape_data.mixed_form = MixedFrameSection::MixedType::Equilibrium;
      } else {
        opserr << OpenSees::PromptValueError 
               << "invalid mixed type value: " << argv[iarg + 1]
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      iarg += 2;
    }
    else if (strcmp(argv[iarg], "-align") == 0) {
      // -align {yy yz zy zz}
      if (iarg + 1 >= argc) {
        opserr << OpenSees::PromptValueError 
               << "shearAlign matrix missing after -align"
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      shape_data.mixed.shear_align = MatrixND<2,2>();
      int argc_sa;
      TCL_Char** argv_sa;
      if (Tcl_SplitList(interp, argv[iarg + 1], &argc_sa,  &argv_sa) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "failed to parse shearAlign matrix after -align"
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      if (argc_sa != 4) {
        opserr << OpenSees::PromptValueError 
               << "shearAlign matrix after -align needs 4 entries"
               << OpenSees::SignalMessageEnd;
        Tcl_Free((char*)argv_sa);
        return TCL_ERROR;
      }
      for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
          double val;
          if (Tcl_GetDouble(interp, argv_sa[i * 2 + j], &val) != TCL_OK) {
            opserr << OpenSees::PromptValueError
                   << "invalid shearAlign matrix entry"
                   << OpenSees::SignalMessageEnd;
            Tcl_Free((char*)argv_sa);
            return TCL_ERROR;
          }
          (*shape_data.mixed.shear_align)(i, j) = val;
        }
      }
      Tcl_Free((char*)argv_sa);
      iarg += 2;
    }
    else if (strcmp(argv[iarg], "-shiftTwist") == 0) {
      // -shiftTwist {tz tx}
      if (iarg + 1 >= argc) {
        opserr << OpenSees::PromptValueError 
               << "shiftTwist vector missing after -shiftTwist"
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      shape_data.mixed.shift_twist = VectorND<2>();
      int argc_st;
      TCL_Char** argv_st;
      if (Tcl_SplitList(interp, argv[iarg + 1], &argc_st,  &argv_st) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "failed to parse shiftTwist vector after -shiftTwist"
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      if (argc_st != 2) {
        opserr << OpenSees::PromptValueError 
               << "shiftTwist vector after -shiftTwist needs 2 entries"
               << OpenSees::SignalMessageEnd;
        Tcl_Free((char*)argv_st);
        return TCL_ERROR;
      }
      for (int i = 0; i < 2; ++i) {
        double val;
        if (Tcl_GetDouble(interp, argv_st[i], &val) != TCL_OK) {
          opserr << OpenSees::PromptValueError
                 << "invalid shiftTwist vector entry"
                 << OpenSees::SignalMessageEnd;
          Tcl_Free((char*)argv_st);
          return TCL_ERROR;
        }
        (*shape_data.mixed.shift_twist)(i) = val;
      }
      Tcl_Free((char*)argv_st);
      iarg += 2;
    }
    else if (strcmp(argv[iarg], "-shiftAxial") == 0) {
      // -shiftAxial {az ax}
      if (iarg + 1 >= argc) {
        opserr << OpenSees::PromptValueError 
               << "shiftAxial vector missing after -shiftAxial"
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      shape_data.mixed.shift_axial = VectorND<2>();
      int argc_sa;
      TCL_Char** argv_sa;
      if (Tcl_SplitList(interp, argv[iarg + 1], &argc_sa,  &argv_sa) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "failed to parse shiftAxial vector after -shiftAxial"
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      if (argc_sa != 2) {
        opserr << OpenSees::PromptValueError 
               << "shiftAxial vector after -shiftAxial needs 2 entries"
               << OpenSees::SignalMessageEnd;
        Tcl_Free((char*)argv_sa);
        return TCL_ERROR;
      }
      for (int i = 0; i < 2; ++i) {
        double val;
        if (Tcl_GetDouble(interp, argv_sa[i], &val) != TCL_OK) {
          opserr << OpenSees::PromptValueError
                 << "invalid shiftAxial vector entry"
                 << OpenSees::SignalMessageEnd;
          Tcl_Free((char*)argv_sa);
          return TCL_ERROR;
        }
        (*shape_data.mixed.shift_axial)(i) = val;
      }
      Tcl_Free((char*)argv_sa);
      iarg += 2;
    }

    else if (strcmp(argv[iarg], "-torsion") == 0 && iarg + 1 < argc) {
      int torsionTag = 0;
      if (Tcl_GetInt(interp, argv[iarg + 1], &torsionTag) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid torsionTag";
        return TCL_ERROR;
      }

      torsion = builder->getTypedObject<UniaxialMaterial>(torsionTag);
      if (torsion == nullptr) {
        opserr << OpenSees::PromptValueError << "uniaxial material does not exist\n";
        opserr << "uniaxial material: " << torsionTag;
        opserr << "\nFiberSection3d: " << secTag << "\n";
        return TCL_ERROR;
      }

      iarg += 2;
    }

    else if (strstr(argv[1], "Asym") != nullptr && !shearParsed) {
      if (iarg + 1 >= argc) {
        opserr << OpenSees::PromptValueError << "Asym sections require shear center before fiber block.\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[iarg], &shape_data.mixed.shear_origin[0]) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid Ys";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[iarg+1], &shape_data.mixed.shear_origin[1]) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid Zs";
        return TCL_ERROR;
      }
      shearParsed = true;
      iarg  += 2;
    }

    else {
      // braces; skip and handle later
      break;
    }
  }

  if (torsion == nullptr && ndm == 3 && !options.isNew && !options.isND) {
    opserr << OpenSees::PromptValueError
           << "missing required torsion for 3D fiber section, use -GJ or "
              "-torsion\n";
    return TCL_ERROR;
  }
  
  // try preallocating fiber space.
  if (iarg < argc) {
    std::string_view sv{argv[iarg]};
    options.reserve = std::count(sv.begin(), sv.end(), '\n');
  }

  // initialize  the fiber section (for building)
  if (initSectionCommands(clientData, interp, secTag, torsion, shape_data, options) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "error constructing the section\n";
    return TCL_ERROR;
  }

  //
  // Execute the commands inside the braces (fibers, patches, and reinforcing layers)
  //
#if !defined(OPS_API)
  if (iarg < argc && Tcl_Eval(interp, argv[iarg]) != TCL_OK) {
    // Assume the subcommands have printed a message regarding the error
    return TCL_ERROR;
  }
#endif
  if (deleteTorsion)
    delete torsion;

  return TCL_OK;
}


//
// add patch to fiber section
//
int
TclCommand_addPatch(ClientData clientData, 
                    Tcl_Interp *interp, 
                    int argc,
                    TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  ModelRegistry* builder = static_cast<ModelRegistry*>(clientData);

  SectionBuilder* fiberSectionRepr = findSectionBuilder(builder, interp, argc, argv);
  if (fiberSectionRepr == nullptr) {
    opserr << OpenSees::PromptValueError << "cannot retrieve section\n";
    return TCL_ERROR;
  }


  // make sure at least one other argument to contain patch type
  if (argc < 2) {
    opserr << OpenSees::PromptValueError << "need to specify a patch type \n";
    return TCL_ERROR;
  }

  // check argv[1] for type of patch  and create the object
  if (strcmp(argv[1], "quad") == 0 || strcmp(argv[1], "quadr") == 0) {
    int numSubdivIJ, numSubdivJK, matTag;
    MatrixND<4,2> vertexCoords{};

    if (argc < 13) {
      opserr << OpenSees::PromptValueError << "invalid number of parameters: patch quad matTag "
                "numSubdivIJ numSubdivJK yVertI zVertI yVertJ zVertJ yVertK "
                "zVertK yVertL zVertL\n";
      return TCL_ERROR;
    }

    int argi = 2;

    if (Tcl_GetInt(interp, argv[argi++], &matTag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid material tag\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &numSubdivIJ) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid numSubdivIJ\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &numSubdivJK) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid numSubdivJK\n";
      return TCL_ERROR;
    }

    for (int j = 0; j < 4; j++) {
      double vertexCoordY, vertexCoordZ;
      if (Tcl_GetDouble(interp, argv[argi++], &vertexCoordY) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid Coordinate y\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[argi++], &vertexCoordZ) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid Coordinate z\n";
        return TCL_ERROR;
      }

      vertexCoords(j, 0) = vertexCoordY;
      vertexCoords(j, 1) = vertexCoordZ;
    }

    // Done parsing
    QuadFiberPatch patch(matTag, numSubdivIJ, numSubdivJK, vertexCoords);
    int error = fiberSectionRepr->addPatch(patch);
    if (error != 0)
      return TCL_ERROR;

  }

  // check argv[1] for type of patch  and create the object
  else if (strcmp(argv[1], "rect") == 0 ||
           strcmp(argv[1], "rectangular") == 0) {

    int numSubdivIJ, numSubdivJK, matTag;
    MatrixND<4,2> vertexCoords;

    if (argc < 9) {
      opserr << OpenSees::PromptValueError 
             << "invalid number of parameters: patch quad matTag "
                "numSubdivIJ numSubdivJK yVertI zVertI yVertK zVertK\n";
      return TCL_ERROR;
    }

    int argi = 2;
    if (Tcl_GetInt(interp, argv[argi++], &matTag) != TCL_OK) {
      opserr << OpenSees::PromptValueError 
             << "invalid matTag:  " << argv[argi-1]
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &numSubdivIJ) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid numSubdivIJ: patch quad matTag numSubdivIJ "
                "numSubdivJK yVertI zVertI yVertJ zVertJ yVertK zVertK yVertL "
                "zVertL\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &numSubdivJK) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid numSubdivJK\n";
      return TCL_ERROR;
    }

    for (int j = 0; j < 2; j++) {
      double vertexCoordY, vertexCoordZ;
      if (Tcl_GetDouble(interp, argv[argi++], &vertexCoordY) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "invalid Coordinate y\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[argi++], &vertexCoordZ) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "invalid Coordinate z\n";
        return TCL_ERROR;
      }

      vertexCoords(j * 2, 0) = vertexCoordY;
      vertexCoords(j * 2, 1) = vertexCoordZ;
    }

    vertexCoords(1, 0) = vertexCoords(2, 0);
    vertexCoords(1, 1) = vertexCoords(0, 1);
    vertexCoords(3, 0) = vertexCoords(0, 0);
    vertexCoords(3, 1) = vertexCoords(2, 1);

    // create and add patch
    QuadFiberPatch patch(matTag, numSubdivIJ, numSubdivJK, vertexCoords);
    int error = fiberSectionRepr->addPatch(patch);
    if (error) {
      return TCL_ERROR;
    }
  }

  else if (strcmp(argv[1], "circ") == 0) {
    int numSubdivRad, numSubdivCirc, matTag;
    double yCenter, zCenter;
    Vector centerPosition(2);
    double intRad, extRad;
    double startAng, endAng;

    int argi = 2;
    if (argc < 11) {
      opserr << OpenSees::PromptValueError << "invalid number of parameters: patch circ matTag "
                "numSubdivCirc numSubdivRad yCenter zCenter intRad extRad "
                "startAng endAng\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &matTag) != TCL_OK) {
      opserr << OpenSees::PromptValueError 
             << "invalid matTag: patch circ matTag numSubdivCirc "
                "numSubdivRad yCenter zCenter intRad extRad startAng endAng\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &numSubdivCirc) != TCL_OK) {
      opserr
          << OpenSees::PromptValueError << "invalid numSubdivCirc\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &numSubdivRad) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid numSubdivRad\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &yCenter) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid yCenter\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &zCenter) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid zCenter\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &intRad) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid intRad\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &extRad) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid extRad\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &startAng) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid startAng\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &endAng) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid endAng\n";
      return TCL_ERROR;
    }

    centerPosition(0) = yCenter;
    centerPosition(1) = zCenter;

    // create patch
    CircPatch patch(matTag, numSubdivCirc, numSubdivRad, centerPosition,
                    intRad, extRad, startAng, endAng);

    // add patch to section
    int error = fiberSectionRepr->addPatch(patch);
    if (error) {
      return TCL_ERROR;
    }
  }

  else {
    opserr << OpenSees::PromptValueError << "patch type is not available\n";
    return TCL_ERROR;
  }

  return TCL_OK;
}

// Add a fiber to a fiber section
int
TclCommand_addFiber(ClientData clientData, Tcl_Interp *interp, int argc,
                    TCL_Char ** const argv)
{
  enum class Position : int {
    Y, Z, Area, Material, End
  };
  ArgumentTracker<Position> tracker;
  std::set<int> positional;

  assert(clientData != nullptr);
  ModelRegistry* builder = static_cast<ModelRegistry*>(clientData);

  SectionBuilder* fiberSectionRepr = findSectionBuilder(builder, interp, argc, argv);
  if (fiberSectionRepr == nullptr) {
    opserr << OpenSees::PromptValueError << "cannot retrieve a section builder\n";
    return TCL_ERROR;
  }

  double yLoc, zLoc=0, area;
  int matTag;
  bool warn_2d_z = false;
  static constexpr int WarpModeCount = 3;
  Vector3D warp[WarpModeCount]{};
  int warp_arg = -1;
  for (int i=1; i<argc; i++) {
    if (strcmp(argv[i], "-section") == 0) {
      ++i;
    }
    else if (strcmp(argv[i], "-warp") == 0) {
      if (i + 1 >= argc) {
        opserr << OpenSees::PromptValueError << "missing warp argument\n";
        return TCL_ERROR;
      }
      warp_arg = i+1;
      i++;
    }
    else if (strcmp(argv[i], "-warn-2d-z") == 0) {
        warn_2d_z = true;
        i++;
    }
    else if (strcmp(argv[i], "-material") == 0) {
      if (argc == ++i || Tcl_GetInt(interp, argv[i], &matTag) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid material tag\n";
        return TCL_ERROR;
      }
      tracker.consume(Position::Material);
    }
    else if ((strcmp(argv[i], "-area") == 0) || (strcmp(argv[i],"-A") == 0)) {
      if (argc == ++i || Tcl_GetDouble(interp, argv[i], &area) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "invalid area: " << argv[i]
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      tracker.consume(Position::Area);
    }
    else if (strcmp(argv[i], "-y") == 0) {
      if (argc == ++i || Tcl_GetDouble(interp, argv[i], &yLoc) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid y coordinate\n";
        return TCL_ERROR;
      }
      tracker.consume(Position::Y);
    }
    else if (strcmp(argv[i], "-z") == 0) {
      if (argc == ++i || Tcl_GetDouble(interp, argv[i], &zLoc) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid z coordinate\n";
        return TCL_ERROR;
      }
      tracker.consume(Position::Z);
    }
    else {
      positional.insert(i);
    }
  }

  for (int i: positional) {
    switch (tracker.current()) {
      case Position::Y:
        if (Tcl_GetDouble(interp, argv[i], &yLoc) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid y coordinate\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::Y);
        break;
      case Position::Z:
        if (Tcl_GetDouble(interp, argv[i], &zLoc) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid z coordinate\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::Z);
        break;
      case Position::Area:
        if (Tcl_GetDouble(interp, argv[i], &area) != TCL_OK) {
          opserr << OpenSees::PromptValueError 
                 << "invalid area: " << argv[i]
                 << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        }
        tracker.consume(Position::Area);
        break;
      case Position::Material:
        if (Tcl_GetInt(interp, argv[i], &matTag) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid material tag\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::Material);
        break;
      default:
        opserr << OpenSees::PromptValueError << "unexpected argument at position " << i << "\n";
        return TCL_ERROR;
    }
  }

  if (builder->getNDM() == 2) {
    if (warn_2d_z && zLoc != 0.0) {
      opswrn << OpenSees::SignalWarning << "z coordinate ignored in 2D\n";
    }
  }
  if (tracker.current() != Position::End) {
    opserr << OpenSees::PromptValueError << "missing required arguments: ";
    while (tracker.current() != Position::End) {
      switch (tracker.current()) {
        case Position::Y:
          opserr << "y ";
          break;
        case Position::Z:
          opserr << "z ";
          break;
        case Position::Area:
          opserr << "area ";
          break;
        case Position::Material:
          opserr << "material ";
          break;
        case Position::End:
          break;
      }
      if (tracker.current() == Position::End)
        break;
      tracker.consume(tracker.current());
    }
    opserr << "\n";
    return TCL_ERROR;
  }

  //
  // process warping
  //
  int i_warp = 0;
  int          split_1_argc;
  const char **split_1_argv = nullptr;
  if (warp_arg >= 0 && Tcl_SplitList(interp, argv[warp_arg], &split_1_argc, &split_1_argv) == TCL_OK) {
    int argi = 0;
    for (; i_warp<WarpModeCount; i_warp++) {
      if (argi >= split_1_argc)
        break;

      //
      int          split_argc;
      const char **split_argv;
      if (Tcl_SplitList(interp, split_1_argv[argi], &split_argc, &split_argv) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid warp\n";
        return TCL_ERROR;
      }

      if (split_argc != 3) {
        opserr << "WARNING warp parameter expected list of 3 floats\n";
          Tcl_Free((char *) split_argv);
          return TCL_ERROR;
      }

      for (int j = 0; j < 3; j++) {
        double val;
        if (Tcl_GetDouble(interp, split_argv[j], &val) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid warp\n";
          Tcl_Free((char *) split_argv);
          return TCL_ERROR;
        }
        warp[i_warp][j] = val;
      }

      // Free memory allocated by Tcl_SplitList.
      Tcl_Free((char *) split_argv);
      argi++;
    }
  }

  //
  // Add fiber to section builder
  //
  int ndm = builder->getNDM();
  int id = -1;
  if (ndm == 2) {
    Vector pos(2);
    pos(0) = yLoc;
    pos(1) = zLoc;
    id = fiberSectionRepr->addFiber(0, matTag, area, pos);
  } else if (ndm == 3) {
    Vector pos(2);
    pos(0) = yLoc;
    pos(1) = zLoc;
    id = fiberSectionRepr->addFiber(0, matTag, area, pos);
  }
  if (id < 0) {
    opserr << OpenSees::PromptValueError << "Failed to add fiber to section\n";
    return TCL_ERROR;
  }

  // set warping
  while (i_warp > 0) {
    if (0 > fiberSectionRepr->setWarping(id, i_warp-1, warp[i_warp-1])) {
      opserr << OpenSees::PromptValueError << "failed to set warping for fiber\n";
      return TCL_ERROR;
    }
    i_warp--;
  }

  Tcl_SetObjResult(interp, Tcl_NewIntObj(id));

  return TCL_OK;
}



// add layers of reinforcing bars to fiber section

int
TclCommand_addFiberLayer(ClientData clientData, Tcl_Interp *interp, int argc,
                         TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  ModelRegistry *builder = static_cast<ModelRegistry*>(clientData);

  SectionBuilder* fiberSectionRepr = findSectionBuilder(builder, interp, argc, argv);
  if (fiberSectionRepr == nullptr) {
    opserr << OpenSees::PromptValueError << "cannot retrieve section\n";
    return TCL_ERROR;
  }

  // make sure at least one other argument to contain layer type
  if (argc < 2) {
    opserr << OpenSees::PromptValueError << "need to specify a layer type \n";
    return TCL_ERROR;
  }

  // check argv[1] for type of layer and create the object
  if (strcmp(argv[1], "straight") == 0 ||
      strcmp(argv[1], "line")     == 0) {
    if (argc < 9) {
      opserr << OpenSees::PromptValueError << "invalid number of parameters: layer straight matTag "
                "numReinfBars reinfBarArea yStartPt zStartPt yEndPt zEndPt\n";
      return TCL_ERROR;
    }

    int matTag, numReinfBars;
    double area;
    VectorND<2> startPt{};
    VectorND<2> endPt{};
  
    int argi = 2;

    if (Tcl_GetInt(interp, argv[argi++], &matTag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid matTag\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &numReinfBars) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid numReinfBars\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &area) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid area\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &startPt[0]) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid yStartPt\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &startPt[1]) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid zStartPt\n"
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &endPt[0]) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid yEndPt\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &endPt[1]) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid zEndPt\n";
      return TCL_ERROR;
    }

    //

    StraightFiberLayer reinfLayer(matTag, numReinfBars, area, startPt, endPt);

    // add reinfLayer to section
    int error = fiberSectionRepr->addLayer(reinfLayer);
    if (error)
      return TCL_ERROR;

  } else if (strcmp(argv[1], "circ") == 0) {
    if (argc < 8) {
      opserr << OpenSees::PromptValueError << "invalid number of parameters: layer circ matTag "
                "numReinfBars reinfBarArea yCenter zCenter arcRadius <startAng "
                "endAng>\n";
      return TCL_ERROR;
    }

    int  matTag, numReinfBars;
    double area;
    double radius, startAng, endAng;

    VectorND<2> center{};
    int argi = 2;

    if (Tcl_GetInt(interp, argv[argi++], &matTag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid matTag\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[argi++], &numReinfBars) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid numReinfBars\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &area) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid reinfBarArea\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &center[0]) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid yCenter\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &center[1]) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid zCenter\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[argi++], &radius) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid radius\n";
      return TCL_ERROR;
    }

    bool anglesSpecified = false;

    if (argc > 9) {
      if (Tcl_GetDouble(interp, argv[argi++], &startAng) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid startAng\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[argi++], &endAng) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid endAng\n";
        return TCL_ERROR;
      }

      anglesSpecified = true;
    }

    // create the reinforcing layer

    int error = -1;
    // construct and add to section
    if (anglesSpecified) {
      // Construct arc
      CircFiberLayer reinfLayer(matTag, numReinfBars, area,
                                center, radius, startAng, endAng);
      error = fiberSectionRepr->addLayer(reinfLayer);
    } else {
      // Construct circle
      CircFiberLayer reinfLayer(matTag, numReinfBars, area,
                                center, radius);
      error = fiberSectionRepr->addLayer(reinfLayer);
    }

    if (error)
      return TCL_ERROR;

  } else {
    opserr << OpenSees::PromptValueError
           << "reinforcing layer type is not available\n";
    return TCL_ERROR;
  }

  return TCL_OK;
}



// add Hfiber to fiber section
int
TclCommand_addHFiber(ClientData clientData, Tcl_Interp *interp, int argc,
                     TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  ModelRegistry* builder = static_cast<ModelRegistry*>(clientData);


  SectionBuilder* fiberSectionRepr = findSectionBuilder(builder, interp, argc, argv);
  if (fiberSectionRepr == nullptr) {
    opserr << OpenSees::PromptValueError << "cannot retrieve section\n";
    return TCL_ERROR;
  }

  // make sure at least one other argument to contain patch type
  if (argc < 5) {
    opserr << OpenSees::PromptValueError << "invalid num args: Hfiber yLoc zLoc area matTag\n";
    return TCL_ERROR;
  }


  int matHTag;
  double yHLoc, zHLoc, Harea;

  if (Tcl_GetDouble(interp, argv[1], &yHLoc) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid yLoc: Hfiber yLoc zLoc area matTag\n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[2], &zHLoc) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid zLoc.\n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &Harea) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid area.\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[4], &matHTag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid matTag.\n";
    return TCL_ERROR;
  }

  Vector fiberHPosition(2);
  fiberHPosition(0) = yHLoc;
  fiberHPosition(1) = zHLoc;

  // add patch to section builder
  int error = fiberSectionRepr->addHFiber(0, matHTag, Harea, fiberHPosition);

  if (error)
    return TCL_ERROR;

  return TCL_OK;
}
