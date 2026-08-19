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
#include <string.h>
#include <Logging.h>
#include <Parsing.h>
#include <UmfpackGenLinSOE.h>
#include <UmfpackGenLinSolver.h>
#include <UmfpackSolver02.h>
#include <UmfpackLinSOE02.h>


LinearSOE*
TclDispatch_newUmfpackLinearSOE(ClientData clientData,
                                Tcl_Interp* interp,
                                Tcl_Size argc, const char** const argv)
{
  int factLVALUE = 10;
  int factorOnce = 0;
  bool doDet = false;

  int count = 2;
  while (count < argc) {
    if ((strcmp(argv[count], "-lValueFact") == 0) ||
        (strcmp(argv[count], "-lvalueFact") == 0) ||
        (strcmp(argv[count], "-LVALUE") == 0)) {
      if (count+1 < argc && Tcl_GetInt(interp, argv[count + 1], &factLVALUE) != TCL_OK)
        return nullptr;
      count++;
    }
    else if ((strcmp(argv[count], "-factorOnce") == 0) ||
               (strcmp(argv[count], "-FactorOnce") == 0)) {
      factorOnce = 1;
      count++;
    }
    else if ((strcmp(argv[count], "-printTime") == 0) ||
                (strcmp(argv[count], "-time") == 0)) {
      count++;
    }
    else if (strcmp(argv[count], "-det") == 0) {
      doDet = true;
      count++;
    }
  }

//  return new UmfpackGenLinSOE(*theSolver, factLVALUE, factorOnce, false);
  if (strstr(argv[1], "02") != nullptr) {
    return new UmfpackLinSOE02(*new UmfpackSolver02(doDet));
  }
  return new UmfpackGenLinSOE(*new UmfpackGenLinSolver(doDet));
}

