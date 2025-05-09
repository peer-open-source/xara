//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the function to parse the Tcl input
//              for prismatic frame elements
//
// Written: fmk,cmp
// Created: 07/99
//
#if 1
#include <array>
#include <vector>
#include <stdlib.h>
#include <string.h>

#include <tcl.h>
#include <Logging.h>
#include <ArgumentTracker.h>

#include <Domain.h>
#include <BasicFrameTransf.h>
#include <LinearFrameTransf.hpp>

#include "ElasticBeam2d.h"
#include "ElasticBeam3d.h"
#include "PrismFrame2d.h"
#include "PrismFrame3d.h"

#include <FrameSection.h>

#include <BasicModelBuilder.h>
#endif
//
// element NAME 
//    $tag $iNode $jNode
//    $E $G $A $J $Iy $Iz $Avy $Avz
//    -section $tag
//    [$tag | -transform $tag | -vert $x $y $z]
//    < -delta $delta >
//    < -shear $Avy $Avz [$G] >
//    < -mass $m> < -cMass >
//    < -warp $warp>
//
// For a two-dimensional problem:
//   element elasticBeamColumn $eleTag $iNode $jNode $A $E $Iz $transfTag <-mass $massDens> <-cMass>
//
// For a three-dimensional problem:
//   element elasticBeamColumn $eleTag $iNode $jNode $A $E $G $J $Iy $Iz $transfTag <-mass $massDens> <-cMass>
//
// For two- and three-dimensional problems:
//   element elasticBeamColumn $eleTag $iNode $jNode $secTag $transfTag <-mass $massDens> <-cMass>
//


template <typename Position, typename Section>
int
Parse_ElasticBeam(ClientData clientData, Tcl_Interp *interp, int argc,
                               TCL_Char ** const argv)
{

  // The program is fundamentally incorrect if this function
  // is invoked with NULL clientData
  assert(clientData != nullptr);
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;

  //
  // Declare the variables we're parsing for
  //
  struct FrameData {
    double E, G;                         // Material data
    double A, Iz, Iy, Ixy, J, Cw;        // Section data
    double thermal_coeff = 0.0,          // Thermal
           thermal_depth = 0.0;
  } beam_data;
  double mass = 0.0;
  bool use_mass = false;


  struct FrameFlags {
    int                          shear_flag = 0;
    int                          geom_flag  = 0;
    int                          relz_flag  = 0;
    int                          rely_flag  = 0;
    enum  {NoWarp=0, HaveWarp}   warp_flag  = NoWarp;
    enum  {
      RigidMass=0,
      PrismMass
    } mass_type  = RigidMass;
  } options;

  ArgumentTracker<Position> tracker;
  std::vector<int> positions;



  FrameTransform2d *theTrans2d = nullptr;
  FrameTransform3d *theTrans3d = nullptr;
  Section* theSection     = nullptr;



  int ndm = builder->getNDM();
  int ndf = builder->getNDF();

  //
  // Parse common parameters
  //
  if (argc < 5) {
    opserr << OpenSees::PromptValueError << "insufficient number of arguments.\n";
    return TCL_ERROR;
  }

  // Parse tag, iNode, jNode
  int tag, iNode, jNode, transTag;
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid tag " << argv[2] << "\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid iNode " << argv[3] << "\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid jNode " << argv[4] << "\n";
    return TCL_ERROR;
  }

  //
  // Parse keyword arguments
  //
  for (int argi = 5; argi < argc; argi++) {

    // Section
    if (strcmp(argv[argi], "-section") == 0) {
      if (argc < argi + 2) {
        opserr << OpenSees::PromptValueError
               << "not enough arguments, expected -section $section?\n";
        return TCL_ERROR;
      }

      int section;
      if (Tcl_GetInt(interp, argv[argi+1], &section) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "invalid section tag.\n";
        return TCL_ERROR;
      }

      theSection = builder->getTypedObject<Section>(section);
      if (theSection == nullptr)
        return TCL_ERROR;

      argi += 1;
      // Fast-forward through section parameters
      while (tracker.current() != Position::Transform && 
             tracker.current() != Position::End)
        tracker.increment();
    }

    // Release
    else if ((strcmp(argv[argi], "-release") == 0) ||
             (strcmp(argv[argi], "-releasez") == 0)) {
      if (argc < argi + 2) {
        opserr << OpenSees::PromptValueError
               << "not enough arguments, expected -release $release?\n";
        return TCL_ERROR;
      }

      if (Tcl_GetInt(interp, argv[argi+1], &options.relz_flag) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "invalid release flag.\n";
        return TCL_ERROR;
      }
      argi += 1;
    }
    else if ((strcmp(argv[argi], "-releasey") == 0) ||
             (strcmp(argv[argi], "-releasey") == 0)) {
      if (argc < argi + 2) {
        opserr << OpenSees::PromptValueError
               << "not enough arguments, expected -releasey $release?\n";
        return TCL_ERROR;
      }

      if (Tcl_GetInt(interp, argv[argi+1], &options.rely_flag) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "invalid release flag.\n";
        return TCL_ERROR;
      }
      argi += 1;
    }

    // Transform
    else if (strcmp(argv[argi], "-transform") == 0) {
      if (argc < argi + 2) {
        opserr << OpenSees::PromptValueError
               << "not enough arguments, expected -transform $transform?\n";
        return TCL_ERROR;
      }

      int transform;
      if (Tcl_GetInt(interp, argv[argi+1], &transform) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "invalid transform tag.\n";
        return TCL_ERROR;
      }
      if (ndm == 2) {
        theTrans2d = builder->getTypedObject<FrameTransform2d>(transform);
        if (theTrans2d == nullptr)
          return TCL_ERROR;
      }
      else if (ndm == 3) {
        theTrans3d = builder->getTypedObject<FrameTransform3d>(transform);
        if (theTrans3d == nullptr)
          return TCL_ERROR;
      }
      argi += 1;
      tracker.consume(Position::Transform);
    }
    else if (strcmp(argv[argi], "-vertical") == 0) {
      //
      // -vertical {$x $y $z}
      //
      if (argc < argi + 2) {
        opserr << OpenSees::PromptValueError
               << "not enough arguments, expected -vertical $vertical?\n";
        return TCL_ERROR;
      }
      int yargc;
      TCL_Char ** yargv;
      if (Tcl_SplitList(interp, argv[argi+1], &yargc, &yargv) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "could not split list\n";
        return TCL_ERROR;
      }
      Vector3D vertical;
      for (int j=0; j<3; j++)
        if (Tcl_GetDouble(interp, yargv[j], &vertical[j]) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "could not read vector\n";
          return TCL_ERROR;
        }
      theTrans3d = new BasicFrameTransf3d(new LinearFrameTransf<2,6>(0, vertical));
      builder->addTaggedObject<FrameTransform3d>(*theTrans3d);
      Tcl_Free((char *)yargv);
      argi += 1;
      tracker.consume(Position::Transform);
    }

    // Mass density
    else if (strcmp(argv[argi], "-mass") == 0) {
      if (argc < argi + 2) {
        opserr << OpenSees::PromptValueError << "not enough arguments, expected -mass $mass?\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[argi + 1], &mass) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "invalid mass\n";
        return TCL_ERROR;
      }
      // If mass is specified directly in the element command, tell the element
      // not to ask the cross section for it's mass.
      use_mass = true;
      argi += 1;
    }

    // Shear flag
    else if (strcmp(argv[argi], "-shear") == 0) {
      if (argc < argi + 2) {
        opserr << OpenSees::PromptValueError
               << "not enough arguments, expected -shear $flag\n";
        return TCL_ERROR;
      }
      if (Tcl_GetInt(interp, argv[argi + 1], &options.shear_flag) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "invalid shear flag\n";
        return TCL_ERROR;
      }
      argi += 1;
    }

    // Geometry flag
    else if (strcmp(argv[argi], "-order") == 0) {
      if (argc < argi + 2) {
        opserr << OpenSees::PromptValueError
               << "not enough arguments, expected -delta $delta\n";
        return TCL_ERROR;
      }
      if (Tcl_GetInt(interp, argv[argi + 1], &options.geom_flag) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "invalid delta\n";
        return TCL_ERROR;
      }
      argi += 1;
    }

    // Mass flags
    else if ((strcmp(argv[argi], "-lMass") == 0) ||
               (strcmp(argv[argi], "lMass") == 0)) {
      options.mass_type = FrameFlags::RigidMass;
    }
    else if ((strcmp(argv[argi], "-cMass") == 0) ||
               (strcmp(argv[argi], "cMass") == 0)) {
      options.mass_type = FrameFlags::PrismMass;
    }

    // Thermal
    else if (strcmp(argv[argi], "-alpha") == 0) {
      if (argc == ++argi ||
          Tcl_GetDouble(interp, argv[argi], &beam_data.thermal_coeff) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "invalid thermal alpha" << tag
               << " iNode jNode A E I alpha depth \n";
        return TCL_ERROR;
      }
    }

    else if (strcmp(argv[argi], "-depth") == 0) {
      if (argc == ++argi ||
          Tcl_GetDouble(interp, argv[argi], &beam_data.thermal_depth) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "invalid thermal depth\n";
        return TCL_ERROR;
      }
    } 

    else
      positions.push_back(argi);
  }

  //
  // Parse positional arguments
  //

  // Case 1
  if (positions.size() <= 2) {
    std::size_t i=0;
    // Section
    if (theSection == nullptr) {
      if (i >= positions.size()) {
        opserr << OpenSees::PromptValueError
               << "Missing required positional argument $section\n";
        return TCL_ERROR;
      }
      int section;
      if (Tcl_GetInt(interp, argv[positions[i]], &section) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid secTag\n";
        return TCL_ERROR;
      }
      theSection = builder->getTypedObject<FrameSection>(section);
      if (theSection == nullptr)
        return TCL_ERROR;

      // Fast-forward through section parameters
      while (tracker.current() != Position::Transform && 
             tracker.current() != Position::End)
        tracker.increment();
      i++;
    }

    // Coordinate transformation
    if (tracker.contains(Position::Transform)) {
      if (i >= positions.size()) {
        opserr << OpenSees::PromptValueError
               << "Missing required positional argument $transform\n";
        return TCL_ERROR;
      }
      if (Tcl_GetInt(interp, argv[positions[i]], &transTag) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid transTag " << argv[positions[i]] << "\n";
        return TCL_ERROR;
      }

      // Check that the builder has a transform with tag; error will be
      // printed from builder
      if (ndm == 2) {
        theTrans2d = builder->getTypedObject<FrameTransform2d>(transTag);
        if (theTrans2d == nullptr)
          return TCL_ERROR;
      }
      else if (ndm == 3) {
        theTrans3d = builder->getTypedObject<FrameTransform3d>(transTag);
        if (theTrans3d == nullptr)
          return TCL_ERROR;
      }
      tracker.consume(Position::Transform);
    }
  }

  else {
    //
    // Parse remaining section properties as positional arguments
    //
    for (int i : positions) {

      switch (tracker.current()) {
        case Position::E:
          if (Tcl_GetDouble(interp, argv[i], &beam_data.E) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid E\n";
              return TCL_ERROR;
          } else {
            tracker.increment();
            break;
          }

        case Position::G:
          if (Tcl_GetDouble(interp, argv[i], &beam_data.G) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid G\n";
              return TCL_ERROR;
          } else {
            tracker.increment();
            break;
          }

        case Position::A:
          if (Tcl_GetDouble(interp, argv[i], &beam_data.A) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid A\n";
              return TCL_ERROR;
          } else {
            tracker.increment();
            break;
          }

        case Position::Iz:
          if (Tcl_GetDouble(interp, argv[i], &beam_data.Iz) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid Iz\n";
              return TCL_ERROR;
          }
          tracker.increment();
          break;

        case Position::Iy:
          if (Tcl_GetDouble (interp, argv[i], &beam_data.Iy) != TCL_OK || beam_data.Iy < 0) {
              opserr << OpenSees::PromptParseError << "invalid Iy\n";
              return TCL_ERROR;
          }
          tracker.increment();
          break;


        case Position::J:
          if (Tcl_GetDouble (interp, argv[i], &beam_data.J) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid J\n";
              return TCL_ERROR;
          }
          tracker.increment();
          break;

        case Position::Transform:
          if (Tcl_GetInt(interp, argv[i], &transTag) != TCL_OK) {
            opserr << OpenSees::PromptValueError << "invalid transTag\n";
            return TCL_ERROR;
          }

          if (ndm == 2) {
            theTrans2d = builder->getTypedObject<FrameTransform2d>(transTag);
            if (theTrans2d == nullptr)
              return TCL_ERROR;
          }
          else if (ndm == 3) {
            theTrans3d = builder->getTypedObject<FrameTransform3d>(transTag);
            if (theTrans3d == nullptr)
              return TCL_ERROR;
          }

          tracker.increment();
          break;

        case Position::End:
          opserr << OpenSees::PromptParseError
                 << "unexpected argument " << argv[i] << "\n";
          return TCL_ERROR;
      }
    }
  }

  // Ensure that all positional arguments have been consumed
  if (tracker.current() != Position::End) {
      opserr << OpenSees::PromptParseError
             << "Missing required positional arguments\n";

      while (tracker.current() != Position::End) {
        switch (tracker.current()) {
          case Position::E:
            opserr << "E ";
            break;
          case Position::A:
            opserr << "A ";
            break;
          case Position::Iz:
            opserr << "Iz ";
            break;
          case Position::Iy:
            opserr << "Iy ";
            break;
          case Position::G:
            opserr << "G ";
            break;
          case Position::J:
            opserr << "J ";
            break;
          case Position::Transform:
            opserr << "transform ";
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


  Element *theBeam = nullptr;

  if (ndm == 2) {
    // check plane frame problem has 3 dof per node
    if (ndf != 3) {
      opserr << OpenSees::PromptValueError << "invalid ndf: " << ndf;
      opserr << ", for plane problem need 3 \n";
      return TCL_ERROR;
    }

    if (theSection != nullptr) {
      // now create the beam and add it to the Domain

      theBeam = new PrismFrame2d(tag, iNode, jNode, 
                                 *theSection, *theTrans2d, 
                                 beam_data.thermal_coeff, beam_data.thermal_depth, 
                                 mass,
                                 options.mass_type,
                                 use_mass,
                                 options.relz_flag,
                                 options.geom_flag);
    }

    else if (strcmp(argv[1], "PrismFrame") == 0) {

      theBeam = new PrismFrame2d(tag, beam_data.A, beam_data.E, beam_data.Iz, 
                                 iNode, jNode, *theTrans2d,
                                 beam_data.thermal_coeff, beam_data.thermal_depth, 
                                 mass, options.mass_type, 
                                 options.relz_flag, 
                                 options.geom_flag);
    } else {

      theBeam = new ElasticBeam2d(tag, beam_data.A, beam_data.E, beam_data.Iz, 
                                  iNode, jNode, *theTrans2d,
                                  beam_data.thermal_coeff, beam_data.thermal_depth, 
                                  mass, options.mass_type,
                                  options.relz_flag);
    }

  } else {
    // ndm == 3

    //
    // Some final validation
    //
    // check space frame problem has 6 dof per node
    if ((ndf != 6 && !options.warp_flag) || (ndf != 7 && options.warp_flag)) {
      opserr << OpenSees::PromptValueError << "invalid ndof: " << ndf;
      opserr << ", for 3d problem  need 6\n";
      return TCL_ERROR;
    }


    if (theSection != nullptr) {
      // now create the beam and add it to the Domain

      std::array<int, 2> nodes {iNode, jNode};
      theBeam = new PrismFrame3d(tag, nodes, 
                                 *theSection, *theTrans3d, 
                                 mass,
                                 options.mass_type,
                                 use_mass,
                                 options.relz_flag, 
                                 options.rely_flag,
                                 options.geom_flag,
                                 options.shear_flag);

    } else {
      if (strcmp(argv[1], "PrismFrame") == 0) {
        std::array<int, 2> nodes {iNode, jNode};
        theBeam = new PrismFrame3d(tag, 
                                   nodes,
                                   beam_data.A, beam_data.E, beam_data.G, 
                                   beam_data.J, beam_data.Iy, beam_data.Iz,
                                   *theTrans3d, mass,
                                   options.mass_type,
                                   options.relz_flag, 
                                   options.rely_flag,
                                   options.geom_flag);

      } else {
        theBeam = new ElasticBeam3d(tag, 
                                    beam_data.A, beam_data.E, beam_data.G, 
                                    beam_data.J, beam_data.Iy, beam_data.Iz,
                                    iNode, jNode,
                                    *theTrans3d, 
                                    mass, options.mass_type);
      }
    }
  }

  // now add the beam to the domain
  if (builder->getDomain()->addElement(theBeam) == false) {
    opserr << OpenSees::PromptValueError 
           << "could not add beam to domain.\n";
    delete theBeam;
    return TCL_ERROR;
  }

  return TCL_OK;
}


int
TclBasicBuilder_addElasticBeam(ClientData clientData, Tcl_Interp *interp, int argc,
                               TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;

  const int ndm = builder->getNDM();

  enum class Args2D : int {
    A, E, Iz, Transform, End, 
    // values coming after End wont be handled
    // by ArgumentTracker; these need to be here for
    // the template to compile.
    G, J, Iy
  };

  enum class Args3D : int {
    A, E, G, J, Iy, Iz, Transform, End
  };

  if (ndm == 2)
    return Parse_ElasticBeam<Args2D, FrameSection>(clientData, interp, argc, argv);

  else
    return Parse_ElasticBeam<Args3D, FrameSection>(clientData, interp, argc, argv);
}

