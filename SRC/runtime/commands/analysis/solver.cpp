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
// Description: This file implements commands that configure the linear solver.
//
#include <string>
#include <algorithm>
#ifdef _MSC_VER 
#  include <string.h>
#  define strcasecmp _stricmp
#else
#  include <strings.h>
#endif

#include <Logging.h>
#include <Parsing.h>
#include "solver.hpp"
#include "BasicAnalysisBuilder.h"

// system of eqn and solvers
#include <SProfileSPDLinSolver.h>
#include <SProfileSPDLinSOE.h>
#include <ProfileSPDLinDirectThreadSolver.h>
#include <SparseGenColLinSOE.h>
#include <SparseGenRowLinSOE.h>
#include <SymSparseLinSOE.h>
#include <SymSparseLinSolver.h>


#if defined(_PARALLEL_PROCESSING)
//  parallel soe & solvers
#  include <DistributedBandSPDLinSOE.h>
#  include <DistributedSparseGenColLinSOE.h>
#  include <DistributedSparseGenRowLinSOE.h>
#  include <DistributedBandGenLinSOE.h>
#  include <DistributedDiagonalSOE.h>
#  include <DistributedDiagonalSolver.h>

#  include <MPIDiagonalSOE.h>
#  include <MPIDiagonalSolver.h>

#  define MPIPP_H
#  include <DistributedSuperLU.h>
#  include <DistributedProfileSPDLinSOE.h>
#endif

using namespace Xara;

static LinearSOE*
G3Parse_newLinearSOE(ClientData, Tcl_Interp* interp, ArgSize, G3_Char **const);

Tcl_CmdProc ParseAppleSparse;

int
XaraCmd_systemSize(ClientData clientData, Tcl_Interp *interp, ArgSize argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  LinearSOE *theSOE = ((BasicAnalysisBuilder *)clientData)->getLinearSOE();

  if (theSOE == nullptr) {
    opserr << "No system has been set";
    return TCL_OK;
  }

  Tcl_SetObjResult(interp, Tcl_NewIntObj(theSOE->getNumEqn()));

  return TCL_OK;
}


int
specifySysOfEqnTable(ClientData clientData, Tcl_Interp *interp, ArgSize argc, G3_Char ** const argv)
{
  // make sure at least one other argument to contain type of system
  if (argc < 2) {
    opserr << OpenSees::PromptValueError
           << "need to specify a system type" 
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

#if defined(XARA_HAVE_APPLE_SPARSE)
  if (strcasecmp(argv[1], "AppleSparse") == 0) {
    return ParseAppleSparse(clientData, interp, argc, argv);
  }
#endif

  LinearSOE* theSOE = G3Parse_newLinearSOE(clientData, interp, argc, argv);

  if (theSOE == nullptr)
    return TCL_ERROR;

  BasicAnalysisBuilder* builder = (BasicAnalysisBuilder*)clientData;

  builder->set(theSOE);
  return TCL_OK;

}


static LinearSOE*
G3Parse_newLinearSOE(ClientData clientData, Tcl_Interp* interp, 
                     Tcl_Size argc,
                     G3_Char ** const argv)
{
  G3_Runtime* rt = G3_getRuntime(interp); 

  // argc is checked by calling function;

  LinearSOE *theSOE = nullptr;
  std::string sys_name{argv[1]};
  transform(sys_name.begin(), sys_name.end(), sys_name.begin(), ::tolower);
  auto ctor = soe_table.find(sys_name);

  if (ctor != soe_table.end()) {
    return ctor->second.ss(rt, argc, argv);
  }

#if defined(XARA_USE_MUMPS)
  else if (strcasecmp(argv[1], "mumps") == 0) {
    return TclDispatch_newMumpsLinearSOE(clientData, interp, argc, argv);
  }
#endif

  else if (strcasecmp(argv[1], "Umfpack")==0) {
    // TODO: if "umfpack" is in solver.hpp, this wont be reached
    return TclDispatch_newUmfpackLinearSOE(clientData, interp, argc, argv);
  }

#if 0
  else if (strcmp(argv[2],"Thread") == 0) {  
    int blockSize = 4;
    int numThreads = 1;
    if (argc == 5) {
      if (Tcl_GetInt(interp, argv[3], &blockSize) != TCL_OK)
        return nullptr; //TCL_ERROR;
      if (Tcl_GetInt(interp, argv[4], &numThreads) != TCL_OK)
        return nullptr; //TCL_ERROR;
    }
    return new ProfileSPDLinSOE(
        *new ProfileSPDLinDirectThreadSolver(numThreads,blockSize,1.0e-12)
    );
  } 
#endif

#if defined(OPS_PETSC)
  else if (strcmp(argv[1], "petsc")==0 ||
           strcmp(argv[1], "Petsc")==0) {
    theSOE = TclDispatch_newPetscSOE(clientData, interp, argc, argv);
  }
#endif

  else {
    opserr << OpenSees::PromptValueError 
           << " system '" 
           << argv[1] << "' is unknown or not installed"
           << OpenSees::SignalMessageEnd;
    return nullptr;
  }

  return theSOE;
}


LinearSOE*
specify_SparseSPD(G3_Runtime *rt, int argc, G3_Char ** const argv)
{
  Tcl_Interp *interp = G3_getInterpreter(rt);

  // determine ordering scheme
  //   1 -- MMD
  //   2 -- ND
  //   3 -- RCM

  int lSparse = 1;
  if (argc == 3) {
    if (Tcl_GetInt(interp, argv[2], &lSparse) != TCL_OK)
      return nullptr;
  }

  return new SymSparseLinSOE(*new SymSparseLinSolver(), lSparse);
}


#ifdef _THREADS
#  include "contrib/sys_of_eqn/ThreadedSuperLU/ThreadedSuperLU.h"
#else
#  include <SuperLU.h>
#endif
// TODO(cmp): Threaded SuperLU?


LinearSOE*
specifySparseGen(G3_Runtime* rt, int argc, G3_Char ** const argv)
{
  Tcl_Interp *interp = G3_getInterpreter(rt);

  SparseGenColLinSolver *theSolver = nullptr;
  int count = 2;
  double thresh = 0.0;
  int npRow = 1;
  int npCol = 1;
  int np = 1;
  char symmetric = 'N';
  double drop_tol = 0.0;
  int permSpec = 1;
  int panelSize = 6;
  int relax = 6;
  LinearSolveSpec spec {
    .renumber = LinearSolveSpec::Renumber::Automatic,
    .symmetry = LinearSolveSpec::Symmetry::General,
    .definite = LinearSolveSpec::Definite::Indefinite
  };

  // defaults for threaded SuperLU
  while (count < argc) {
    if ((strcmp(argv[count], "p")    == 0) ||
        (strcmp(argv[count], "piv")  == 0) ||
        (strcmp(argv[count], "-piv") == 0)) {
      thresh = 1.0;
      count++;
    }
    else if ((strcmp(argv[count], "-np") == 0) ||
                (strcmp(argv[count], "np")  == 0)) {
      count++;
      if (count >= argc) {
        opserr << OpenSees::PromptValueError
               << "np option requires an argument"
               << OpenSees::SignalMessageEnd;
        return nullptr;
      }
      if (Tcl_GetInt(interp, argv[count], &np) != TCL_OK)
        return nullptr;
      count++;
    }

    else if ((strcmp(argv[count], "npRow")  == 0) ||
                (strcmp(argv[count], "-npRow") == 0)) {
      count++;
      if (count >= argc) {
        opserr << OpenSees::PromptValueError
               << "npRow option requires an argument"
               << OpenSees::SignalMessageEnd;
        return nullptr;
      }
      if (Tcl_GetInt(interp, argv[count], &npRow) != TCL_OK)
        return nullptr;
      count++;
    }
    else if ((strcmp(argv[count], "npCol")  == 0) ||
                (strcmp(argv[count], "-npCol") == 0)) {
      count++;
      if (count >= argc) {
        opserr << OpenSees::PromptValueError
               << "npCol option requires an argument"
               << OpenSees::SignalMessageEnd;
        return nullptr;
      }
      if (Tcl_GetInt(interp, argv[count], &npCol) != TCL_OK)
        return nullptr;
      count++;
    }
    else if (strcmp(argv[count], "-symmetric") == 0) {
      if (count +1 < argc) {
        if (strcmp(argv[count+1], "False") == 0) {
          symmetric = 'N';
          spec.symmetry = LinearSolveSpec::Symmetry::General;
          count++;
        }
        else {
          symmetric = 'Y';
          spec.symmetry = LinearSolveSpec::Symmetry::Symmetric;
          // Only increase count if next arg is "True", otherwise
          // argv[count+1] belongs to another option
          if (strcmp(argv[count+1], "True") == 0)
            count++;
        }
      }
      else {
        // -symmetric is passed without an argument, assume True
        symmetric = 'Y';
        spec.symmetry = LinearSolveSpec::Symmetry::Symmetric;
      }
      count++;
    }
  }


#ifdef _THREADS
  if (np != 0)
    theSolver = new ThreadedSuperLU(np, permSpec, panelSize, relax, thresh);
  else
    return nullptr;
// #endif
// 
// #ifdef _PARALLEL_PROCESSING
//     if (theSolver != 0)
//       delete theSolver;
//     theSolver = 0;
// 
//     if (npRow != 0 && npCol != 0) {
//       theSolver = new DistributedSuperLU(npRow, npCol);
//     }
#else
  // TODO(cmp) : SuperLU
  theSolver = new SuperLU(permSpec, drop_tol, panelSize, relax, symmetric, spec);
#endif

#ifdef _PARALLEL_PROCESSING
  return new DistributedSparseGenColLinSOE(*theSolver);
#else
  return new SparseGenColLinSOE(*theSolver);
#endif
}
