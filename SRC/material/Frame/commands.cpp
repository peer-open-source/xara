//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Written: cmp
//
// TODO:
// section ElasticFrame 1 -linearize 2 
// 
#include <tcl.h>
#include <string.h>
#include <ArgumentTracker.h>
#include <Parsing.h>
#include <Logging.h>
#include <BasicModelBuilder.h>
//
#include <ElasticLinearFrameSection3d.h>
#include <ElasticSection3d.h>
#include <SectionAggregator.h>


static int
validate(FrameSectionConstants& section, int ndm, int shear)
{
  int status = 0;
  if (section.A <= 0.0) {
    status -= 1;
    //opserr << "A <= 0.0\n";
  }
  
  if (section.Iz <= 0.0) {
    status -= 1;
    //opserr << "Iz <= 0.0\n";
  }
  
  
//if (shear)
//  if (section.G <= 0.0) {
//    status -= 1;
//    //opserr << "G <= 0.0\n";
//  }
  
  if (ndm == 3) {
//  if (section.J <= 0.0) {
//    //opserr << "J <= 0.0\n";
//    status -= 1;
//  }
    if (section.Iy <= 0.0)  {
      //opserr << "Iy <= 0.0\n";
      status -= 1;
    }
  }
  return status;
}


// section ElasticFrame tag E? A? Iz? <Iy? G? J?>
//                          E  A  I
//                          E  A  I  G  J
//    -E     $E
//    -G     $G
//
//    -A     $A
//          {$A $Ay}      if ndm == 2
//          {$A $Ay $Az}
//    -Ay    $Ay
//    -Az    $Az
//    -I/B   $Iz
//          {$Iy $Iz}
//          {$Iy $Iz $Iyz}
//    -J     $J
//    -Cw    $Cw
//          {$Cw $Ca}
//    -Q    {$Qy $Qz <$Qyx $Qyz>}
//    -R    {$Qy $Qz}
//
int
TclCommand_newElasticSection(ClientData clientData, Tcl_Interp *interp,
                                  int argc, TCL_Char ** const argv)
{
    enum class Position : int {
      Tag, E, A, Iz, Iy, G, J, End
    };

    assert(clientData != nullptr);
    BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);
    FrameSectionConstants consts {};

    ArgumentTracker<Position> tracker;
    std::set<int> positional;

    FrameSection* theSection = nullptr;
    bool construct_full = false;

    if (argc < 5) {
        opserr << OpenSees::PromptParseError << "insufficient arguments\n";
        opserr << "Want: section Elastic tag? E? A? Iz? <Iy? G? J?>.\n";
        return TCL_ERROR;
    }

    int tag;
    double E, G, J;

    bool use_mass = false;
    double mass=0.0;

    if (builder->getNDM() == 2) {
        tracker.consume(Position::G);
        tracker.consume(Position::J);
        tracker.consume(Position::Iy);
    }

    int i;
    for (i=2; i<argc; i++) {

      if (strcmp(argv[i], "-mass") == 0) {
        if (argc == ++i || Tcl_GetDouble (interp, argv[i], &mass) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid mass.\n";
          return TCL_ERROR;
        }
        use_mass = true;
      }
      else if ((strcmp(argv[i], "-youngs-modulus") == 0) ||
               (strcmp(argv[i], "-E") == 0)) {
        if (argc == ++i || Tcl_GetDouble (interp, argv[i], &E) != TCL_OK || E <= 0.0) {
          opserr << OpenSees::PromptParseError << "invalid Young's modulus.\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::E);
      }
      else if ((strcmp(argv[i], "-shear-modulus") == 0) ||
               (strcmp(argv[i], "-G") == 0)) {
        if (argc == ++i || Tcl_GetDouble (interp, argv[i], &G) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid shear modulus.\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::G);
      }
      //
      // Section constants
      //
      else if ((strcmp(argv[i], "-area") == 0) ||
               (strcmp(argv[i], "-A") == 0)) {
        if (argc == ++i || Tcl_GetDouble (interp, argv[i], &consts.A) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid area.\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::A);
      }

      else if ((strcmp(argv[i], "-shear-y") == 0) ||
               (strcmp(argv[i], "-Ay") == 0)) {
        if (argc == ++i || Tcl_GetDouble (interp, argv[i], &consts.Ay) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid shear area.\n";
          return TCL_ERROR;
        }
        construct_full = true;
      }

      else if ((strcmp(argv[i], "-shear-z") == 0) ||
               (strcmp(argv[i], "-Az") == 0)) {
        if (argc == ++i || Tcl_GetDouble (interp, argv[i], &consts.Az) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid shear area.\n";
          return TCL_ERROR;
        }
        construct_full = true;
      }

      else if ((strcmp(argv[i], "-inertia") == 0) ||
               (strcmp(argv[i], "-I") == 0) ||
               (strcmp(argv[i], "-Iz") == 0)) {
        if (argc == ++i || Tcl_GetDouble (interp, argv[i], &consts.Iz) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid inertia Iz\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::Iz);
      }

      else if ((strcmp(argv[i], "-inertia-y") == 0) ||
               (strcmp(argv[i], "-Iy") == 0)) {
        if (argc == ++i || Tcl_GetDouble (interp, argv[i], &consts.Iy) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid inertia Iy\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::Iy);
      }

      else if ((strcmp(argv[i], "-inertia-yz") == 0) ||
               (strcmp(argv[i], "-Iyz") == 0)) {
        if (argc == ++i || Tcl_GetDouble (interp, argv[i], &consts.Iyz) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid product of inertia\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::Iy);
      }

      else if ((strcmp(argv[i], "-venant") == 0) ||
               (strcmp(argv[i], "-J") == 0)) {
        if (argc == ++i || Tcl_GetDouble (interp, argv[i], &J) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid St. Venant constant.\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::J);
      }

      else if (strcmp(argv[i], "-Qy")==0) {
        if (argc == ++i || Tcl_GetDouble (interp, argv[i], &consts.Qy) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid Qy.\n";
          return TCL_ERROR;
        }
        construct_full = true;
      }
      else if (strcmp(argv[i], "-Qz")==0) {
        if (argc == ++i || Tcl_GetDouble (interp, argv[i], &consts.Qz) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid Qz.\n";
          return TCL_ERROR;
        }
        construct_full = true;
      }
      else if (strcmp(argv[i], "-Sa")==0) {
        if (argc == ++i || Tcl_GetDouble (interp, argv[i], &consts.Sa) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid Sa.\n";
          return TCL_ERROR;
        }
        construct_full = true;
      }
      else if (strcmp(argv[i], "-Sy")==0) {
        if (argc == ++i || Tcl_GetDouble (interp, argv[i], &consts.Sy) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid Sy.\n";
          return TCL_ERROR;
        }
        construct_full = true;
      }
      else if (strcmp(argv[i], "-Sz")==0) {
        if (argc == ++i || Tcl_GetDouble (interp, argv[i], &consts.Sz) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid Sz.\n";
          return TCL_ERROR;
        }
        construct_full = true;
      }
      else if (strcmp(argv[i], "-Rw")==0) {
        if (argc == ++i || Tcl_GetDouble (interp, argv[i], &consts.Rw) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid Rw.\n";
          return TCL_ERROR;
        }
        construct_full = true;
      }
      else if (strcmp(argv[i], "-Ry")==0) {
        if (argc == ++i || Tcl_GetDouble (interp, argv[i], &consts.Ry) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid Ry.\n";
          return TCL_ERROR;
        }
        construct_full = true;
      }
      else if (strcmp(argv[i], "-Rz")==0) {
        if (argc == ++i || Tcl_GetDouble (interp, argv[i], &consts.Rz) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid Rz.\n";
          return TCL_ERROR;
        }
        construct_full = true;
      }
      else if (strcmp(argv[i], "-Cw")==0) {
        if (argc == ++i || Tcl_GetDouble (interp, argv[i], &consts.Cw) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid Cw.\n";
          return TCL_ERROR;
        }
        construct_full = true;
      }
      else
        positional.insert(i);

    }

    //
    // Positional arguments
    //
    for (int i : positional) {
      switch (tracker.current()) {
        case Position::Tag :
          if (Tcl_GetInt(interp, argv[i], &tag) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid section Elastic tag.\n";
              return TCL_ERROR;           
          } else {
            tracker.increment();
            break;
          }

        case Position::E:
          if (Tcl_GetDouble (interp, argv[i], &E) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid E.\n";
              return TCL_ERROR;
          } else {
            tracker.increment();
            break;
          }

        case Position::A:
          if (Tcl_GetDouble (interp, argv[i], &consts.A) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid A.\n";
              return TCL_ERROR;
          } else {
            tracker.increment();
            break;
          }

        case Position::Iz:
          if (Tcl_GetDouble (interp, argv[i], &consts.Iz) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid Iz.\n";
              return TCL_ERROR;
          } else {
            tracker.increment();
            break;
          }

        case Position::Iy:
          if (Tcl_GetDouble (interp, argv[i], &consts.Iy) != TCL_OK || consts.Iy < 0) {
              opserr << OpenSees::PromptParseError << "invalid Iy.\n";
              return TCL_ERROR;
          } else {
            tracker.increment();
            break;
          }

        case Position::G:
          if (Tcl_GetDouble (interp, argv[i], &G) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid G.\n";
              return TCL_ERROR;
          } else {
            tracker.increment();
            break;
          }

        case Position::J:
          if (Tcl_GetDouble (interp, argv[i], &J) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid J.\n";
              return TCL_ERROR;
          } else {
            tracker.increment();
            break;
          }

        case Position::End:
          opserr << OpenSees::PromptParseError << "unexpected argument " << argv[i] << ".\n";
          return TCL_ERROR;
      }
    }

    if (tracker.current() != Position::End) {
      opserr << OpenSees::PromptParseError
             << "missing required arguments: ";
      while (tracker.current() != Position::End) {
        switch (tracker.current()) {
          case Position::Tag :
            opserr << "tag ";
            break;
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



    if (construct_full) {
      consts.Ca =   consts.Iy + consts.Iz - J;
      consts.Sa = -(consts.Iy + consts.Iz - J);
      theSection = new ElasticLinearFrameSection3d(tag,
          E, G,
          consts.A,                                     // n-n
          consts.Iy, consts.Iz, consts.Iyz,             // m-m
          consts.Cw, consts.Ca,                         // w-w
          consts.Qy, consts.Qz,                         // n-m
          consts.Rw, consts.Ry, consts.Rz,              // n-w|v
          consts.Sa, consts.Sy, consts.Sz,              // m-v|w
          mass, use_mass                                // mass
      );
    }
    else if (strcmp(argv[1], "Elastic") == 0)
      theSection = new ElasticSection3d(tag, E, consts.A, consts.Iz, consts.Iy, G, J);

    else if (builder->getNDM() == 3) { //argc > 8) {
      theSection = new ElasticLinearFrameSection3d(tag, E,  consts.A,
                                                   consts.Iz, consts.Iy, 
                                                   G, J, mass, use_mass);
    }
    else
      theSection = new ElasticLinearFrameSection3d(tag, E, consts.A, 
                                                   consts.Iz, mass, use_mass);


    if (theSection == nullptr || builder->addTaggedObject<FrameSection>(*theSection) < 0) {
      if (theSection != nullptr)
        delete theSection;
      return TCL_ERROR;
    } else
      return TCL_OK;
}


int
TclCommand_addSectionAggregator(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char**const argv)
{
    assert(clientData != nullptr);
    BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

    if (argc < 5) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      opserr << "Want: section Aggregator tag? uniTag1? code1? ... <-section "
                "secTag?>"
             << endln;
      return TCL_ERROR;
    }

    int tag;
    int secTag;
    FrameSection *theSec = nullptr;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid Aggregator tag" << endln;
      return TCL_ERROR;
    }

    int nArgs = argc - 3;

    for (int ii = 5; ii < argc; ii++) {
      if (strcmp(argv[ii], "-section") == 0 && ++ii < argc) {
        if (Tcl_GetInt(interp, argv[ii], &secTag) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid Aggregator tag" << endln;
          return TCL_ERROR;
        }
        
        theSec = builder->getTypedObject<FrameSection>(secTag);
        if (theSec == 0)
          return TCL_ERROR;
        
        nArgs -= 2;
      }
    }

    int nMats = nArgs / 2;

    if (nArgs % 2 != 0) {
      opserr << G3_ERROR_PROMPT << "improper number of arguments for Aggregator" << endln;
      return TCL_ERROR;
    }

    ID codes(nMats);
    UniaxialMaterial **theMats = new UniaxialMaterial *[nMats];

    int i, j;
    for (i = 3, j = 0; j < nMats; i++, j++) {
      int tagI;
      if (Tcl_GetInt(interp, argv[i], &tagI) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid Aggregator matTag" << endln;
        return TCL_ERROR;
      }

      theMats[j] = builder->getTypedObject<UniaxialMaterial>(tagI);
      if (theMats[j] == 0)
        return TCL_ERROR;

      i++;

      if (strcmp(argv[i], "Mz") == 0)
        codes(j) = SECTION_RESPONSE_MZ;
      else if (strcmp(argv[i], "P") == 0)
        codes(j) = SECTION_RESPONSE_P;
      else if (strcmp(argv[i], "Vy") == 0)
        codes(j) = SECTION_RESPONSE_VY;
      else if (strcmp(argv[i], "My") == 0)
        codes(j) = SECTION_RESPONSE_MY;
      else if (strcmp(argv[i], "Vz") == 0)
        codes(j) = SECTION_RESPONSE_VZ;
      else if (strcmp(argv[i], "T") == 0)
        codes(j) = SECTION_RESPONSE_T;
      else {
        opserr << G3_ERROR_PROMPT << "invalid code" << endln;
        opserr << "\nsection Aggregator: " << tag << endln;
        return TCL_ERROR;
      }
    }

    FrameSection* theSection = nullptr;
    if (theSec)
      theSection = new SectionAggregator(tag, *theSec, nMats, theMats, codes);
    else
      theSection = new SectionAggregator(tag, nMats, theMats, codes);

    // Now add the material to the modelBuilder
    if (builder->addTaggedObject<FrameSection>(*theSection) < 0) {
      delete theSection;
      return TCL_ERROR;
    } else
      return TCL_OK;

    delete[] theMats;
}
