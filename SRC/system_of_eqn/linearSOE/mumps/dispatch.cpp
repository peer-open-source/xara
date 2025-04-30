#include <tcl.h>

#if defined(_PARALLEL_PROCESSING) || defined(_PARALLEL_INTERPRETERS)
#  include <MumpsParallelSOE.h>
#  include <MumpsParallelSolver.h>
#else
#  include <MumpsSOE.h>
#  include <MumpsSolver.h>
#endif

struct MumpsOptions {int icntl14, icntl7;};

//struct MumpsOptions
LinearSOE*
TclDispatch_newMumpsLinearSOE(ClientData clientData, Tcl_Interp* interp, int argc, G3_Char ** const argv)
{
  if (strcmp(argv[1], "Mumps") == 0) {

    int icntl14 = 20;
    int icntl7 = 7;
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
        } else if (strcmp(argv[currentArg], "-ICNTL7") == 0) {
          if (Tcl_GetInt(interp, argv[currentArg + 1], &icntl7) != TCL_OK)
          ;
          currentArg += 2;
        } else if (strcmp(argv[currentArg], "-matrixType") == 0) {
          if (Tcl_GetInt(interp, argv[currentArg + 1], &matType) != TCL_OK)
            opserr << "Mumps Warning: failed to get -matrixType. Unsymmetric "
                      "matrix assumed\n";
          if (matType < 0 || matType > 2) {
            opserr << "Mumps Warning: wrong -matrixType value (" << matType
                   << "). Unsymmetric matrix assumed\n";
            matType = 0;
          }
          currentArg += 2;
        } else
          currentArg++;
      }
    }

#  ifdef _PARALLEL_PROCESSING
      MumpsParallelSolver *theSolver = new MumpsParallelSolver(icntl7, icntl14);
      theSOE = new MumpsParallelSOE(*theSolver);
#  elif _PARALLEL_INTERPRETERS
    if (G3_
      MumpsParallelSolver *theSolver = new MumpsParallelSolver(icntl7, icntl14);
      MumpsParallelSOE *theParallelSOE =
          new MumpsParallelSOE(*theSolver, matType);
      theParallelSOE->setProcessID(OPS_rank);
      theParallelSOE->setChannels(numChannels, theChannels);
      theSOE = theParallelSOE;
#  else
    MumpsSolver *theSolver = new MumpsSolver(icntl7, icntl14);
    theSOE = new MumpsSOE(*theSolver, matType);
#  endif
  }
  return theSOE;
}
