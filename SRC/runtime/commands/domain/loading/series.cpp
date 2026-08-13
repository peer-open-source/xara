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
// Description: This file contains the function invoked when the user invokes
// the pattern command in the interpreter. It is invoked by the
// TclBasicBuilder_addPattern function.
//
// Written: fmk
// Created: 11/00
//
#include <sys/types.h>
#include <sys/stat.h>

#include <set>
#include <vector>
#include <string.h>
#include <Logging.h>
#include <Parsing.h>
#include <ModelRegistry.h>
#include <ArgumentTracker.h>

#include <Domain.h>
#include <LinearSeries.h>
#include <ConstantSeries.h>
#include <PathTimeSeries.h>
#include <PathSeries.h>
#include <TrigSeries.h>
#include <RectangularSeries.h>
#include <PulseSeries.h>
#include <TriangleSeries.h>


extern "C" int OPS_ResetInputNoBuilder(ClientData clientData, 
                                       Tcl_Interp *interp,
                                       int cArg,
                                       int mArg, TCL_Char ** const argv,
                                       Domain *domain);


static void *
TclDispatch_newLinearSeries(ClientData clientData, 
                            Tcl_Interp* interp,
                            Tcl_Size argc, 
                            TCL_Char ** const argv)
{
  int numRemainingArgs = argc;

  int tag = 0;
  double cFactor = 1.0;

  if (numRemainingArgs != 0) {

    if (numRemainingArgs == 1 || numRemainingArgs == 3) {
      if (Tcl_GetInt(interp, argv[0], &tag) != 0) {
        opserr << OpenSees::PromptValueError
               << "invalid series tag in LinearSeries tag? <-factor "
                  "factor?>"
               << OpenSees::SignalMessageEnd;
        return nullptr;
      }
      numRemainingArgs--;
    }

    if (numRemainingArgs > 1) {
      const char *argvS = argv[1];
      if (argvS == 0) {
        opserr << OpenSees::PromptValueError << "string error in LinearSeries with tag: " << tag
               << "\n";
        return nullptr;
      }
      if (Tcl_GetDouble(interp, argv[2], &cFactor) != 0) {
        opserr << OpenSees::PromptValueError << "invalid factor in LinearSeries with tag: " << tag
               << "\n";
        return nullptr;
      }
    }
  }

  return new LinearSeries(tag, cFactor);
}


static TimeSeries *
TclDispatch_newTimeSeries(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{

  // note the 1 instead of usual 2
  OPS_ResetInputNoBuilder(clientData, interp, 1, argc, argv, nullptr);

  TimeSeries *theSeries = nullptr;

  if ((strcmp(argv[0], "Constant") == 0) ||
      (strcmp(argv[0], "ConstantSeries") == 0)) {
    // LoadPattern and ConstantSeries - read args & create LinearSeries object
    double cFactor = 1.0;
        
    int endMarker = 1;
    if ((endMarker != argc) && (strcmp(argv[endMarker], "-factor") == 0)) {
      endMarker++;
      if (endMarker == argc || 
          Tcl_GetDouble(interp, argv[endMarker], &cFactor) != TCL_OK) {

        opserr << OpenSees::PromptValueError
               << "Invalid cFactor " 
               << (endMarker < argc ? argv[endMarker] : "")
               << OpenSees::SignalMessageEnd;
        Tcl_Free((char*)argv);
        return 0;
      }
      endMarker++;
    }

    theSeries = new ConstantSeries(cFactor);
  }

  else if ((strcmp(argv[0],"Trig") == 0 )|| 
           (strcmp(argv[0],"TrigSeries") == 0) ||
           (strcmp(argv[0],"SineSeries") == 0) ||
           (strcmp(argv[0],"Sine") == 0)) {

    // Trig tStart tFinish period <-shift shift> <-factor cFactor>
    enum class Args {
      Tag,
      TStart,
      TFinish,
      Period,
      EndRequired,
      EndPositional,
      Shift,
      Factor,
      End
    };
    ArgumentTracker<Args> tracker;
    std::set<int> positions;
    int tag = 0;
    double cFactor = 1.0;
    double tStart, tFinish, period;
    double shift = 0.0;

    for (int i = 1; i < argc; i++) {
      if (strcmp(argv[i], "-shift") == 0) {
        if (i + 1 == argc) {
          opserr << OpenSees::PromptValueError 
                 << "Missing value for -shift" 
                 << OpenSees::SignalMessageEnd;
          return nullptr;
        }
        if (Tcl_GetDouble(interp, argv[i + 1], &shift) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "Invalid value for -shift: " << argv[i + 1] << "\n";
          return nullptr;
        }
        tracker.consume(Args::Shift);
        i++; // Skip the next argument since it's the value for -shift
      } 
      else if (strcmp(argv[i], "-factor") == 0) {
        if (i + 1 == argc) {
          opserr << OpenSees::PromptValueError 
                << "Missing value for -factor" 
                << OpenSees::SignalMessageEnd;
          return nullptr;
        }
        if (Tcl_GetDouble(interp, argv[i + 1], &cFactor) != TCL_OK) {
          opserr << OpenSees::PromptValueError 
                << "Invalid value for -factor: " << argv[i + 1] 
                << OpenSees::SignalMessageEnd;
          return nullptr;
        }
        tracker.consume(Args::Factor);
        i++; // Skip the next argument since it's the value for -factor
      }
      else if (strcmp(argv[i], "-tag") == 0) {
        if (i + 1 == argc) {
          opserr << OpenSees::PromptValueError 
                 << "Missing value for -tag" 
                 << OpenSees::SignalMessageEnd;
          return nullptr;
        }
        if (Tcl_GetInt(interp, argv[i + 1], &tag) != TCL_OK) {
          opserr << OpenSees::PromptValueError 
                 << "Invalid value for -tag: " << argv[i + 1] 
                 << OpenSees::SignalMessageEnd;
          return nullptr;
        }
        tracker.consume(Args::Tag);
        i++; // Skip the next argument since it's the value for -tag
      }
      else if (strcmp(argv[i], "-tStart") == 0) {
        if (i + 1 == argc) {
          opserr << OpenSees::PromptValueError 
                 << "Missing value for -tStart" 
                 << OpenSees::SignalMessageEnd;
          return nullptr;
        }
        if (Tcl_GetDouble(interp, argv[i + 1], &tStart) != TCL_OK) {
          opserr << OpenSees::PromptValueError 
                 << "Invalid value for -tStart: " << argv[i + 1] 
                 << OpenSees::SignalMessageEnd;
          return nullptr;
        }
        tracker.consume(Args::TStart);
        i++; // Skip the next argument since it's the value for -tStart
      }
      else if (strcmp(argv[i], "-tFinish") == 0) {
        if (i + 1 == argc) {
          opserr << OpenSees::PromptValueError 
                 << "Missing value for -tFinish" 
                 << OpenSees::SignalMessageEnd;
          return nullptr;
        }
        if (Tcl_GetDouble(interp, argv[i + 1], &tFinish) != TCL_OK) {
          opserr << OpenSees::PromptValueError 
                 << "Invalid value for -tFinish: " << argv[i + 1] 
                 << OpenSees::SignalMessageEnd;
          return nullptr;
        }
        tracker.consume(Args::TFinish);
        i++; // Skip the next argument since it's the value for -tFinish
      }
      else if (strcmp(argv[i], "-period") == 0) {
        if (i + 1 == argc) {
          opserr << OpenSees::PromptValueError 
                 << "Missing value for -period" 
                 << OpenSees::SignalMessageEnd;
          return nullptr;
        }
        if (Tcl_GetDouble(interp, argv[i + 1], &period) != TCL_OK) {
          opserr << OpenSees::PromptValueError 
                 << "Invalid value for -period: " << argv[i + 1] 
                 << OpenSees::SignalMessageEnd;
          return nullptr;
        }
        tracker.consume(Args::Period);
        i++; // Skip the next argument since it's the value for -period
      }
      else {
        positions.insert(i);
      }
    }
  
    int missing = 0;
    missing += tracker.contains(Args::TStart);
    missing += tracker.contains(Args::TFinish);
    missing += tracker.contains(Args::Period);

    bool hasPositionalTag = false;

    if (tracker.contains(Args::Tag)) {
      if ((int)positions.size() == missing + 1) {
        hasPositionalTag = true;
      } else if ((int)positions.size() == missing) {
        tracker.consume(Args::Tag); // no tag supplied; default tag = 0
      } else {
        opserr << OpenSees::PromptValueError
              << "wrong number of positional arguments for Trig series"
              << OpenSees::SignalMessageEnd;
        return nullptr;
      }
    } else {
      if ((int)positions.size() != missing) {
        opserr << OpenSees::PromptValueError
              << "unexpected positional argument after -tag"
              << OpenSees::SignalMessageEnd;
        return nullptr;
      }
    }

    // Parse positional arguments
    for (int i : positions) {
      switch (tracker.current()) {
        case Args::Tag:
          if (Tcl_GetInt(interp, argv[i], &tag) != TCL_OK) {
            opserr << OpenSees::PromptValueError 
                   << "Invalid series tag in Trig tag?" 
                   << OpenSees::SignalMessageEnd;
            return nullptr;
          }
          tracker.consume(Args::Tag);
          break;
        case Args::TStart:
          if (Tcl_GetDouble(interp, argv[i], &tStart) != TCL_OK) {
            opserr << OpenSees::PromptValueError 
                   << "Invalid tStart: " << argv[i] 
                   << OpenSees::SignalMessageEnd;
            return nullptr;
          }
          tracker.consume(Args::TStart);
          break;
        case Args::TFinish:
          if (Tcl_GetDouble(interp, argv[i], &tFinish) != TCL_OK) {
            opserr << OpenSees::PromptValueError 
                   << "Invalid tFinish: " << argv[i] 
                   << OpenSees::SignalMessageEnd;
            return nullptr;
          }
          tracker.consume(Args::TFinish);
          break;
        case Args::Period:
          if (Tcl_GetDouble(interp, argv[i], &period) != TCL_OK) {
            opserr << OpenSees::PromptValueError 
                   << "Invalid period: " << argv[i] 
                   << OpenSees::SignalMessageEnd;
            return nullptr;
          }
          tracker.consume(Args::Period);
          break;
        default:
          opserr << OpenSees::PromptValueError 
                 << "Unexpected argument: " << argv[i] 
                 << OpenSees::SignalMessageEnd;
          return nullptr;
      }
    }

    if ((tracker.current() < Args::EndRequired)) {
      opserr << OpenSees::PromptValueError 
             << "Missing required arguments for Trig series: " ;
      while (tracker.current() != Args::EndRequired) {
        switch (tracker.current()) {
          case Args::Tag:
            opserr << "tag ";
            break;
          case Args::TStart:
            opserr << "tStart ";
            break;
          case Args::TFinish:
            opserr << "tFinish ";
            break;
          case Args::Period:
            opserr << "period ";
            break;
          default:
            break;
        }
        tracker.increment();
      }
      opserr << OpenSees::SignalMessageEnd;
      return nullptr;
    }

    // Validation
    if (tFinish < tStart) {
      opserr << OpenSees::PromptValueError 
             << "tFinish must be greater than or equal to tStart" 
             << OpenSees::SignalMessageEnd;
      return nullptr;
    }
    if (period <= 0.0) {
      opserr << OpenSees::PromptValueError 
             << "period must be greater than 0" 
             << OpenSees::SignalMessageEnd;
      return nullptr;
    }
    theSeries = new TrigSeries(tag, tStart, tFinish, period, shift, cFactor);     
   }

  else if ((strcmp(argv[0], "Linear") == 0) ||
           (strcmp(argv[0], "LinearSeries") == 0)) {

    void *theResult = TclDispatch_newLinearSeries(clientData, interp, argc - 1, &argv[1]);

    if (theResult != nullptr)
      theSeries = (TimeSeries *)theResult;
    else
      opserr << "ERROR\n";
  }

  else if (strcmp(argv[0],"Pulse") == 0)  {
    // LoadPattern and PulseSeries - read args & create PulseSeries object
    // Pulse tStart tFinish period <-width pulseWidth> <-shift shift> <-factor cFactor>
    double cFactor = 1.0;
    double tStart, tFinish, period;
    double width = 0.5;
    double shift = 0.0;
      
    if (argc < 4) {
      opserr << "WARNING not enough args for Pulse"
             << OpenSees::SignalMessageEnd;
      Tcl_Free((char*)argv);
      return 0; 
    }
    if (Tcl_GetDouble(interp, argv[1], &tStart) != TCL_OK) {
      opserr << "WARNING invalid tStart " << argv[1]
             << OpenSees::SignalMessageEnd;
      Tcl_Free((char*)argv);
      return 0;                         
    }
    if (Tcl_GetDouble(interp, argv[2], &tFinish) != TCL_OK) {
      opserr << "WARNING invalid tFinish " << argv[2] << "\n";
      Tcl_Free((char*)argv);
      return 0; 
    }
    if (Tcl_GetDouble(interp, argv[3], &period) != TCL_OK) {
      opserr << "WARNING invalid period " << argv[3] << "\n";
      Tcl_Free((char*)argv);
      return 0; 
    }

    int endMarker = 4;

    while (endMarker < argc && endMarker < argc) {
      if (strcmp(argv[endMarker], "-factor") == 0) {
        // allow user to specify the factor
        endMarker++;
        if (endMarker == argc || 
            Tcl_GetDouble(interp, argv[endMarker], &cFactor) != TCL_OK) {
          opserr << "WARNING invalid cFactor " 
                 << (endMarker < argc ? argv[endMarker] : "")
                 << OpenSees::SignalMessageEnd;
          Tcl_Free((char*)argv);
          return 0;
        }
      }

      else if (strcmp(argv[endMarker],"-width") == 0) {
        // allow user to specify pulse width
        endMarker++;
        if (endMarker == argc || 
            Tcl_GetDouble(interp, argv[endMarker], &width) != TCL_OK) {
            
          opserr << "WARNING invalid pulse width " 
                 << (endMarker < argc ? argv[endMarker] : "")
                 << OpenSees::SignalMessageEnd;
          Tcl_Free((char*)argv);
          return 0;
        }
      }

      else if (strcmp(argv[endMarker],"-shift") == 0) {
        // allow user to specify phase shift
        endMarker++;
        if (endMarker == argc || 
            Tcl_GetDouble(interp, argv[endMarker], &shift) != TCL_OK) {
            
          opserr << "WARNING invalid phase shift " 
                 << (endMarker < argc ? argv[endMarker] : "")
                 << OpenSees::SignalMessageEnd;
          Tcl_Free((char*)argv);
          return 0;
        }
      }
      endMarker++;
    }

    theSeries = new PulseSeries(tStart, tFinish, period, width, shift, cFactor);
  }     

  else if (strcmp(argv[0],"Triangle") == 0)  {
    // LoadPattern and TriangleSeries - read args & create TriangleSeries object
    double cFactor = 1.0;
    double tStart, tFinish, period;
    double shift = 0.0;
      
    if (argc < 4) {
      opserr << "WARNING not enough TriangleSeries args"
             << OpenSees::SignalMessageEnd;
      Tcl_Free((char*)argv);
      return 0; 
    }   
    if (Tcl_GetDouble(interp, argv[1], &tStart) != TCL_OK) {
      opserr << "WARNING invalid tStart " 
             << argv[1] 
             << OpenSees::SignalMessageEnd;
      Tcl_Free((char*)argv);
      return 0;                         
    }
    if (Tcl_GetDouble(interp, argv[2], &tFinish) != TCL_OK) {
      opserr << "WARNING invalid tFinish " 
             << argv[2] 
             << OpenSees::SignalMessageEnd;
      Tcl_Free((char*)argv);
      return 0; 
    }
    if (Tcl_GetDouble(interp, argv[3], &period) != TCL_OK) {
      opserr << "WARNING invalid period " << argv[3] << OpenSees::SignalMessageEnd;
      Tcl_Free((char*)argv);
      return 0; 
    }     
    
    int endMarker = 4;
    
    while (endMarker < argc && endMarker < argc) {
      if (strcmp(argv[endMarker],"-factor") == 0) {
        // Scale factor
        endMarker++;
        if (endMarker == argc || 
            Tcl_GetDouble(interp, argv[endMarker], &cFactor) != TCL_OK) {
          
          opserr << "WARNING invalid scale factor " 
                 << (endMarker < argc ? argv[endMarker] : "")
                 << OpenSees::SignalMessageEnd;
          Tcl_Free((char*)argv);
          return 0;
        }
      }

      else if (strcmp(argv[endMarker],"-shift") == 0) {
        // Phase shift
        endMarker++;
        if (endMarker == argc || 
            Tcl_GetDouble(interp, argv[endMarker], &shift) != TCL_OK) {
            
          opserr << "WARNING invalid phase shift "
                 << (endMarker < argc ? argv[endMarker] : "")
                 << OpenSees::SignalMessageEnd;
          Tcl_Free((char*)argv);
          return 0;
        }
      }
      endMarker++;
    }

    theSeries = new TriangleSeries(tStart, tFinish, period, shift, cFactor);    
  }

  else if (strcmp(argv[0],"Rectangular") == 0) {
    // LoadPattern and RectangularSeries
    double tStart, tFinish;
    double cFactor = 1.0;
    if (argc < 3) {
      opserr << "WARNING not enogh args - ";
      opserr << " Rectangular tStart tFinish <-factor cFactor>\n";
      Tcl_Free((char*)argv);
      return 0; 
    }   
    if (Tcl_GetDouble(interp, argv[1], &tStart) != TCL_OK) {
      opserr << "WARNING invalid tStart " << argv[1] << OpenSees::SignalMessageEnd;
      Tcl_Free((char*)argv);
      return 0;
    }
    if (Tcl_GetDouble(interp, argv[2], &tFinish) != TCL_OK) {
      opserr << "WARNING invalid tStart " << argv[2] << OpenSees::SignalMessageEnd;
      Tcl_Free((char*)argv);
      return 0; 
    }     

    int endMarker =  3;
    
    if ((endMarker != argc) && (strcmp(argv[endMarker],"-factor") == 0)) {
      // allow user to specify the factor
      endMarker++;
      if (endMarker == argc || 
          Tcl_GetDouble(interp, argv[endMarker], &cFactor) != TCL_OK) {

        opserr << "WARNING invalid cFactor "
               << (endMarker < argc ? argv[endMarker] : "")
               << OpenSees::SignalMessageEnd;
        Tcl_Free((char*)argv);
        return 0;
      }
      endMarker++;
    }  

    theSeries = new RectangularSeries(tStart, tFinish, cFactor); 
  }

#if 0
  else if ((strcmp(argv[0], "Pulse") == 0) ||
           (strcmp(argv[0], "PulseSeries") == 0)) {

    void *theResult = OPS_PulseSeries(rt, argc, argv);
    if (theResult != 0)
      theSeries = (TimeSeries *)theResult;

  }
#endif

  else if ((strcmp(argv[0], "Series") == 0) || (strcmp(argv[0], "Path") == 0)) {

    double cFactor = 1.0;

    if (argc < 3) {
      opserr << OpenSees::PromptValueError << "not enough args - ";
      opserr << " Series -dt timeIncr -values {list of points }\n";
      return 0;
    }

    int tag = 0;
    double timeIncr = 0.0;
    int endMarker = 1;
    bool done = false;
    int fileName = 0;
    int fileTimeName = 0;
    int filePathName = 0;
    Vector *dataPath = nullptr;
    Vector *dataTime = nullptr;
    bool useLast = false;
    bool prependZero = false;
    double startTime = 0.0;

    struct stat fileInfo;

    if (Tcl_GetInt(interp, argv[endMarker], &tag) == TCL_OK) {
      endMarker++;
    }

    while (endMarker < argc && done == false) {

      if (strcmp(argv[endMarker], "-dt") == 0) {
        // allow user to specify the time increment
        endMarker++;
        if (endMarker == argc ||
            Tcl_GetDouble(interp, argv[endMarker], &timeIncr) != TCL_OK) {

          opserr << OpenSees::PromptValueError 
                 << "invalid dt " << argv[endMarker]
                 << OpenSees::SignalMessageEnd;
          return 0;
        }
      }

      else if (strcmp(argv[endMarker], "-tag") == 0) {
        // allow user to specify the tag
        endMarker++;
        if (endMarker == argc ||
            Tcl_GetInt(interp, argv[endMarker], &tag) != TCL_OK) {

          opserr << OpenSees::PromptValueError
                 << (endMarker < argc ? argv[endMarker] : "")
                 << OpenSees::SignalMessageEnd;
          return 0;
        }
      }

      else if (strcmp(argv[endMarker], "-factor") == 0) {
        // allow user to specify the factor
        endMarker++;
        if (endMarker == argc ||
            Tcl_GetDouble(interp, argv[endMarker], &cFactor) != TCL_OK) {

          opserr << OpenSees::PromptValueError
                 << "invalid scale factor "
                 << (endMarker < argc ? argv[endMarker] : "")
                 << OpenSees::SignalMessageEnd;
          return 0;
        }
      }

      else if (strcmp(argv[endMarker], "-file") == 0) {
        // File name containing time and data points
        endMarker++;
        if (endMarker != argc) {
          fileName = endMarker; // argv[endMarker];
          if (stat(argv[endMarker], &fileInfo ) != 0) {
            opserr << OpenSees::PromptValueError
                   << "Cannot open file "
                   << argv[endMarker]
                   << OpenSees::SignalMessageEnd;
            return nullptr;
          }
        }
        else {
          opserr << OpenSees::PromptValueError
                 << "Missing argument"
                 << OpenSees::SignalMessageEnd;
        }
      }

      else if (strcmp(argv[endMarker], "-filePath") == 0) {
        // File name containing the data points
        endMarker++;
        if (endMarker != argc) {
          filePathName = endMarker; // argv[endMarker];
          if (stat(argv[endMarker], &fileInfo ) != 0) {
            opserr << OpenSees::PromptValueError << "Cannot open file "
                   << argv[endMarker] << "\n";
            return nullptr;
          }
        }
      }

      else if (strcmp(argv[endMarker], "-fileTime") == 0) {
        // File name containing the data points
        endMarker++;
        if (endMarker != argc) {
          fileTimeName = endMarker; // argv[endMarker];
          if (stat(argv[endMarker], &fileInfo ) != 0) {
            opserr << OpenSees::PromptValueError << "Cannot open file "
                   << argv[endMarker] << "\n";
            return nullptr;
          }
        }
      }

      else if (strcmp(argv[endMarker], "-values") == 0) {
        // data points in tcl list
        endMarker++;
        if (endMarker == argc) {
          opserr << OpenSees::PromptValueError 
                 << "no data values after -values flag"
                 << OpenSees::SignalMessageEnd;
          return nullptr;
        }

        int pathSize;
        TCL_Char **pathStrings;
        if (Tcl_SplitList(interp, argv[endMarker], &pathSize, &pathStrings) != TCL_OK) {
          opserr << OpenSees::PromptValueError
                 << "problem splitting path list " << argv[endMarker]
                 << OpenSees::SignalMessageEnd;
          return nullptr;
        }

        double check_double;
        if ((pathSize == 1) && (endMarker +1 < argc) &&
            (Tcl_GetDouble(interp, argv[endMarker + 1], &check_double) == TCL_OK)) {
          // -values is not enclosed in {}
          Tcl_Free((char *)pathStrings);
          int count = 0;
          std::vector<double> values;
          for (int i = endMarker; i < argc; i++) {
            double value;
            if (Tcl_GetDouble(interp, argv[i], &value) != TCL_OK) {
              break;
            }
            values.push_back(value);
            count++;
          }
          dataPath = new Vector(count);
          for (int i = 0; i < count; i++)
            (*dataPath)(i) = values[i];

        }
        else {
          dataPath = new Vector(pathSize);
          for (int i = 0; i < pathSize; ++i) {
            double value;
            if (Tcl_GetDouble(interp, pathStrings[i], &value) != TCL_OK) {
              opserr << OpenSees::PromptValueError
                     << "Invalid path value "
                      << pathStrings[i]
                      << OpenSees::SignalMessageEnd;
              Tcl_Free((char *)pathStrings);
              return nullptr;
            }
            (*dataPath)(i) = value;
          }
          // free up the array of pathsStrings
          Tcl_Free((char *)pathStrings);
        }
      }

      else if (strcmp(argv[endMarker], "-time") == 0) {
        // time points in one of two forms:
        // a) -time 0.0 0.1 ...
        // b) -time {0.0 0.1 ...}
        endMarker++;
        if (endMarker >= argc) {
          opserr << OpenSees::PromptValueError
                 << "Missing required argument for 'time'"
                 << OpenSees::SignalMessageEnd;
          return nullptr;
        }
        
        // check for case a)
        double time_i;
        if (Tcl_GetDouble(interp, argv[endMarker], &time_i) == TCL_OK) {
          // it is case a)
          std::vector<double> time_vector{};
          do {
            time_vector.push_back(time_i);
            if ((endMarker+1 >= argc) || (Tcl_GetDouble(interp, argv[++endMarker], &time_i) != TCL_OK))
              break;
          } while (true);
          dataTime = new Vector(time_vector.size());
          for (std::size_t i=0; i<time_vector.size(); i++)
            (*dataTime)(i) = time_vector[i];
            
          endMarker--;
        }
        else
        {
          // case b), time in a Tcl list
          Tcl_Size pathSize;
          TCL_Char **pathStrings;

          if (Tcl_SplitList(interp, argv[endMarker], &pathSize, &pathStrings) !=
              TCL_OK) {

            opserr << OpenSees::PromptValueError 
                   << "problem spltting time path " << argv[endMarker]
                   << OpenSees::SignalMessageEnd;
            return nullptr;
          }

          dataTime = new Vector(pathSize);
          for (int i = 0; i < pathSize; ++i) {
            double value;
            if (Tcl_GetDouble(interp, pathStrings[i], &value) != TCL_OK) {
              opserr << OpenSees::PromptValueError 
                     << "Invalid time path value "
                     << pathStrings[i]
                     << OpenSees::SignalMessageEnd;
              Tcl_Free((char *)pathStrings);
              return nullptr;
            }
            (*dataTime)(i) = value;
          }
          // free up the array of pathsStrings .. see tcl man pages as to why
          Tcl_Free((char *)pathStrings);
        }
      }

      else if (strcmp(argv[endMarker], "-useLast") == 0) {
        useLast = true;
      }

      else if (strcmp(argv[endMarker], "-prependZero") == 0) {
        prependZero = true;
      }

      else if (strcmp(argv[endMarker], "-startTime") == 0 ||
               strcmp(argv[endMarker], "-tStart") == 0) {
        // start time
        endMarker++;
        if (endMarker == argc ||
            Tcl_GetDouble(interp, argv[endMarker], &startTime) != TCL_OK) {

          opserr << OpenSees::PromptValueError << "invalid tStart " 
                 << (endMarker < argc ? argv[endMarker] : "")
                 << OpenSees::SignalMessageEnd;
          return 0;
        }
      }

      endMarker++;
    }

    if (filePathName != 0 && fileTimeName == 0 && timeIncr != 0.0) {
      theSeries = new PathSeries(tag, argv[filePathName], timeIncr, cFactor,
                                 useLast, prependZero, startTime);
    }

    else if (fileName != 0) {
      theSeries = new PathTimeSeries(tag, argv[fileName], cFactor, useLast);

    } else if (filePathName != 0 && fileTimeName != 0) {
      theSeries = new PathTimeSeries(tag, argv[filePathName],
                                     argv[fileTimeName], cFactor, useLast);

    } else if (dataPath != 0 && dataTime == 0 && timeIncr != 0.0) {
      theSeries = new PathSeries(tag, *dataPath, timeIncr, cFactor, useLast,
                                 prependZero, startTime);
      delete dataPath;

    }
    else if ((dataPath != nullptr) && (dataTime != nullptr)) {
      if (dataTime->Size() != dataPath->Size()) {
        opserr << OpenSees::PromptValueError << "size of time vector (" << dataTime->Size()
               << ") must be equal to size of values (" << dataPath->Size() << ")"
               << OpenSees::SignalMessageEnd;
        delete dataPath;
        delete dataTime;
        return nullptr;
      }
      theSeries =  new PathTimeSeries(tag, *dataPath, *dataTime, cFactor, useLast);
      delete dataPath;
      delete dataTime;

    }
    else {
      opserr << OpenSees::PromptValueError 
             << "choice of options for Path Series invalid - valid "
                "options for ";
      opserr << " Path are\n";
      opserr << " \t -fileT fileTimeName -fileP filePathName \n";
      opserr << " \t -dt constTimeIncr -file filePathName\n";
      opserr << " \t -dt constTimeIncr -values {list of points on path}\n";
      opserr << " \t -time {list of time points} -values {list of points on path}\n";
      return 0;
    }

  }

#if 0
  else if ((strcmp(argv[0], "PeerDatabase") == 0) ||
           (strcmp(argv[0], "PeerMotion") == 0)) {

    void *theResult = OPS_PeerMotion(rt, argc, argv);
    if (theResult != 0)
      theSeries = (TimeSeries *)theResult;

    PeerMotion *thePeerMotion = (PeerMotion *)theSeries;

    if (argc > 4 && theSeries != 0) {
      int argCount = 4;

      while (argCount + 1 < argc) {
        if ((strcmp(argv[argCount], "-dT") == 0) ||
            (strcmp(argv[argCount], "-dt") == 0) ||
            (strcmp(argv[argCount], "-DT") == 0)) {
          const char *variableName = argv[argCount + 1];
          double dT = thePeerMotion->getDt();
          char string[30];
          sprintf(string, "set %s %.18e", variableName, dT);
          if (Tcl_Eval(interp, string) != TCL_OK) {
            opserr << G3_WARN_PROMPT << Tcl_GetStringResult(interp);
            Tcl_Exit(TCL_ERROR);
          }
          argCount += 2;
        } else if ((strcmp(argv[argCount], "-nPts") == 0) ||
                   (strcmp(argv[argCount], "-NPTS") == 0)) {
          const char *variableName = argv[argCount + 1];
          int nPts = thePeerMotion->getNPts();
          char string[30];
          sprintf(string, "set %s %d", variableName, nPts);
          if (Tcl_Eval(interp, string) != TCL_OK) {
            opserr << G3_WARN_PROMPT << Tcl_GetStringResult(interp);
            Tcl_Exit(TCL_ERROR);
          }
          argCount += 2;
        } else
          argCount++;
      }
    }
  }
#endif

  else {
    // type unknown
    opserr << "WARNING unknown Series type " << argv[0] << " - ";
    opserr << " valid types: Linear, Rectangular, Path, Constant, Trig, Sine\n";
    return nullptr;
  }

  return theSeries;
}


TimeSeries *
TclSeriesCommand(ClientData clientData, Tcl_Interp *interp, TCL_Char * const arg)
{
  int argc;
  TCL_Char ** argv;
  TimeSeries *series;
  int timeSeriesTag = 0;

  if (Tcl_GetInt(interp, arg, &timeSeriesTag) == TCL_OK) {
    if (clientData && (series = ((ModelRegistry*)clientData)->getTypedObject<TimeSeries>(timeSeriesTag)))
      return series->getCopy();
  }

  // split the list
  if (Tcl_SplitList(interp, arg, &argc, &argv) != TCL_OK) {
    opserr << "WARNING could not split series list " << arg << "\n";
    return nullptr;
  }

  TimeSeries *theSeries = TclDispatch_newTimeSeries(clientData, interp, argc, argv);

  // clean up after ourselves and return the series
  Tcl_Free((char *)argv);
  return theSeries;
}


int
TclCommand_addTimeSeries(ClientData clientData,
                         Tcl_Interp *interp, 
                         Tcl_Size argc,
                         TCL_Char ** const argv)
{
  TimeSeries *theSeries = TclDispatch_newTimeSeries(clientData, interp, argc - 1, &argv[1]);

  ModelRegistry *builder = static_cast<ModelRegistry*>(clientData);

  if (theSeries != nullptr) {
    int tag;
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "failed to read tag \"" << argv[2] << "\"\n";
      return TCL_ERROR;
    }
    if (builder->addTypedObject<TimeSeries>(tag, theSeries) < 0)
      return TCL_ERROR;
    else
      return TCL_OK;
  }
  return TCL_ERROR;
}

