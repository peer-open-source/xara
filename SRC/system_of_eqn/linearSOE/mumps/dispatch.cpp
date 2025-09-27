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
#include <Logging.h>
#include <Parsing.h>

#if defined(_PARALLEL_PROCESSING) || defined(_PARALLEL_INTERPRETERS)
#  include <MumpsParallelSOE.h>
#  include <MumpsParallelSolver.h>
#else
#  include <MumpsSOE.h>
#  include <MumpsSolver.h>
#endif


struct MumpsOptions {
  int icntl14, icntl7;
};


LinearSOE*
TclDispatch_newMumpsLinearSOE(ClientData clientData, 
                              Tcl_Interp* interp,
                              Tcl_Size argc, 
                              G3_Char ** const argv)
{

  int icntl14 = 20;
  int icntl7  = 7;
  int matType = 0; // 0: unsymmetric, 
                   // 1: symmetric positive definite, 
                   // 2: symmetric general

  int currentArg = 2;
  while (currentArg < argc) {
    if (argc > 2) {
      if (strcmp(argv[currentArg], "-ICNTL14") == 0) {
        if (Tcl_GetInt(interp, argv[currentArg + 1], &icntl14) != TCL_OK)
          ;
        currentArg += 2;
      }

      else if (strcmp(argv[currentArg], "-ICNTL7") == 0) {
        if (argc < currentArg + 1) {
          opserr << OpenSees::PromptValueError
                 << "ICNTL7 option requires an argument\n";
          return nullptr;
        }
        if (Tcl_GetInt(interp, argv[currentArg + 1], &icntl7) != TCL_OK)
        ;
        currentArg += 2;
      }

      else if (strcmp(argv[currentArg], "-matrixType") == 0) {
        if (currentArg + 1 >= argc) {
          opserr << OpenSees::PromptValueError
                 << "matrixType option requires an argument\n";
          return nullptr;
        }
        if (Tcl_GetInt(interp, argv[currentArg + 1], &matType) != TCL_OK)
          opserr << OpenSees::PromptValueError 
                 << "failed to get -matrixType. Unsymmetric "
                    "matrix assumed\n";
        if (matType < 0 || matType > 2) {
          opserr << OpenSees::PromptValueError 
                 << "invalid matrixType (" << matType
                 << "). Unsymmetric matrix assumed\n";
          matType = 0;
        }
        currentArg += 2;
      }
      else
        currentArg++;
    }
  }

  //
  // Construct the SOE
  //
  LinearSOE* theSOE = nullptr;
#if defined(_PARALLEL_PROCESSING)
    MumpsParallelSolver *theSolver = new MumpsParallelSolver(icntl7, icntl14);
    theSOE = new MumpsParallelSOE(*theSolver);
#elif defined(_PARALLEL_INTERPRETERS)
    MumpsParallelSolver *theSolver = new MumpsParallelSolver(icntl7, icntl14);
    MumpsParallelSOE *theParallelSOE =
        new MumpsParallelSOE(*theSolver, matType);
    theParallelSOE->setProcessID(OPS_rank);
    theParallelSOE->setChannels(numChannels, theChannels);
    theSOE = theParallelSOE;
#else
  theSOE = new MumpsSOE(*new MumpsSolver(icntl7, icntl14), matType);
#endif

  return theSOE;
}
