
#include <tcl.h>
#include <assert.h>
#include <string.h>
#include <Logging.h>
#include <Parsing.h>
#include <BoucWenMG.h>
#include <BasicModelBuilder.h>
#include <ArgumentTracker.h>

int
TclCommand_newBoucWenMG(ClientData clientData, Tcl_Interp *interp,
                                  int argc, TCL_Char ** const argv)
{
    assert(clientData != nullptr);
    BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

    // Declare order for positional arguments
    enum class Position : int {
        Tag,
        eta,
        k0,      // True elastic stiffness
        sy0,     // True yield stress
        sig,
        lam,
        mup,
        sigp,
        rsmax,
        n,
        alpha,
        alpha1,
        alpha2,
        betam1,
        End
    };

    // Create tracker for positional arguments
    ArgumentTracker<Position> tracker;

    // Declare parameters
    BoucWenMG::Params params {};

    if (argc < 16) {
        opserr << OpenSees::PromptParseError << "insufficient arguments\n";
        return TCL_ERROR;
    }

    int tag;
    std::set<int> positional;

    int i;
    for (i=2; i<argc; i++) {


      if ((strcmp(argv[i], "-eta") == 0)) {
        if (argc == ++i || Tcl_GetDouble (interp, argv[i], &params.eta) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid eta\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::eta);
      }
      else if (strcmp(argv[i], "-k0") == 0) {
        if (argc == ++i || Tcl_GetDouble(interp, argv[i], &params.k0) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid k0\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::k0);
      }
      else if (strcmp(argv[i], "-sy0") == 0) {
        if (argc == ++i || Tcl_GetDouble(interp, argv[i], &params.sy0) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid sy0\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::sy0);
      }
      else if (strcmp(argv[i], "-sig") == 0) {
        if (argc == ++i || Tcl_GetDouble(interp, argv[i], &params.sig) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid sig\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::sig);
      }
      else if (strcmp(argv[i], "-lam") == 0) {
        if (argc == ++i || Tcl_GetDouble(interp, argv[i], &params.lam) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid lam\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::lam);
      }
      else if (strcmp(argv[i], "-mup") == 0) {
        if (argc == ++i || Tcl_GetDouble(interp, argv[i], &params.mup) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid mup\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::mup);
      }
      else if (strcmp(argv[i], "-sigp") == 0) {
        if (argc == ++i || Tcl_GetDouble(interp, argv[i], &params.sigp) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid sigp\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::sigp);
      }
      else if (strcmp(argv[i], "-rsmax") == 0) {
        if (argc == ++i || Tcl_GetDouble(interp, argv[i], &params.rsmax) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid rsmax\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::rsmax);
      }
      else if (strcmp(argv[i], "-n") == 0) {
        if (argc == ++i || Tcl_GetDouble(interp, argv[i], &params.n) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid n\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::n);
      }
      else if (strcmp(argv[i], "-alpha") == 0) {
        if (argc == ++i || Tcl_GetDouble(interp, argv[i], &params.alpha) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid alpha\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::alpha);
      }
      else if (strcmp(argv[i], "-alpha1") == 0) {
        if (argc == ++i || Tcl_GetDouble(interp, argv[i], &params.alpha1) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid alpha1\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::alpha1);
      }
      else if (strcmp(argv[i], "-alpha2") == 0) {
        if (argc == ++i || Tcl_GetDouble(interp, argv[i], &params.alpha2) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid alpha2\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::alpha2);
      }
      else if (strcmp(argv[i], "-betam1") == 0) {
        if (argc == ++i || Tcl_GetDouble(interp, argv[i], &params.betam1) != TCL_OK) {
          opserr << OpenSees::PromptParseError << "invalid betam1\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::betam1);
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
              opserr << OpenSees::PromptParseError << "invalid tag.\n";
              return TCL_ERROR;           
          } else {
            tracker.increment();
            break;
          }

        case Position::eta:
          if (Tcl_GetDouble (interp, argv[i], &params.eta) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid eta.\n";
              return TCL_ERROR;
          } else {
            tracker.increment();
            break;
          }
        case Position::k0:
          if (Tcl_GetDouble(interp, argv[i], &params.k0) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid k0.\n";
              return TCL_ERROR;
          } else {
            tracker.increment();
            break;
          }
        case Position::sy0:
          if (Tcl_GetDouble(interp, argv[i], &params.sy0) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid sy0.\n";
              return TCL_ERROR;
          } else {
            tracker.increment();
            break;
          }
        
        case Position::sig:
          if (Tcl_GetDouble(interp, argv[i], &params.sig) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid sig.\n";
              return TCL_ERROR;
          } else {
            tracker.increment();
            break;
          }
        
        case Position::lam:
          if (Tcl_GetDouble(interp, argv[i], &params.lam) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid lam.\n";
              return TCL_ERROR;
          } else {
            tracker.increment();
            break;
          }
        
        case Position::mup:
          if (Tcl_GetDouble(interp, argv[i], &params.mup) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid mup.\n";
              return TCL_ERROR;
          } else {
            tracker.increment();
            break;
          }
        
        case Position::sigp:
          if (Tcl_GetDouble(interp, argv[i], &params.sigp) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid sigp.\n";
              return TCL_ERROR;
          } else {
            tracker.increment();
            break;
          }
        
        case Position::rsmax:
          if (Tcl_GetDouble(interp, argv[i], &params.rsmax) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid rsmax.\n";
              return TCL_ERROR;
          } else {
            tracker.increment();
            break;
          }
        
        case Position::n:
          if (Tcl_GetDouble(interp, argv[i], &params.n) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid n.\n";
              return TCL_ERROR;
          } else {
            tracker.increment();
            break;
          }
        
        case Position::alpha:
          if (Tcl_GetDouble(interp, argv[i], &params.alpha) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid alpha.\n";
              return TCL_ERROR;
          } else {
            tracker.increment();
            break;
          }
        
        case Position::alpha1:
          if (Tcl_GetDouble(interp, argv[i], &params.alpha1) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid alpha1.\n";
              return TCL_ERROR;
          } else {
            tracker.increment();
            break;
          }

        case Position::alpha2:
          if (Tcl_GetDouble(interp, argv[i], &params.alpha2) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid alpha2.\n";
              return TCL_ERROR;
          } else {
            tracker.increment();
            break;
          }
        
        case Position::betam1:
          if (Tcl_GetDouble(interp, argv[i], &params.betam1) != TCL_OK) {
              opserr << OpenSees::PromptParseError << "invalid betam1.\n";
              return TCL_ERROR;
          } else {
            tracker.increment();
            break;
          }

        case Position::End:
          opserr << OpenSees::PromptParseError << "unexpected argument" << argv[i] << ".\n";
          return TCL_ERROR;
      }
    }


    //
    //
    //
    if (tracker.current() != Position::End) {
      opserr << OpenSees::PromptParseError
             << "missing required arguments: ";
      while (tracker.current() != Position::End) {
        switch (tracker.current()) {
          case Position::Tag :
            opserr << "tag ";
            break;
          case Position::eta:
            opserr << "eta ";
            break;
          case Position::k0:
            opserr << "k0 ";
            break;
          case Position::sy0:
            opserr << "sy0 ";
            break;
          case Position::sig:
            opserr << "sig ";
            break;
          case Position::lam:
            opserr << "lam ";
            break;
          case Position::mup:
            opserr << "mup ";
            break;
          case Position::sigp:
            opserr << "sigp ";
            break;
          case Position::rsmax:
            opserr << "rsmax ";
            break;
          case Position::n:
            opserr << "n ";
            break;
          case Position::alpha:
            opserr << "alpha ";
            break;
          case Position::alpha1:
            opserr << "alpha1 ";
            break;
          case Position::alpha2:
            opserr << "alpha2 ";
            break;
          case Position::betam1:
            opserr << "betam1 ";
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
    //
    //

    params.eta2 = 1.0 - params.eta;

    UniaxialMaterial* theMaterial = new BoucWenMG(tag, params);


    if (theMaterial == nullptr || builder->addTaggedObject<UniaxialMaterial>(*theMaterial) < 0) {
      if (theMaterial != nullptr)
        delete theMaterial;
      return TCL_ERROR;

    } else
      return TCL_OK;
}
