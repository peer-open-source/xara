//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, Claudio M. Perez
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
#include <string.h>
#include <cassert>
#include <Logging.h>
#include <Parsing.h>
#include <BasicAnalysisBuilder.h>
#include "CreateAppleSparse.h"


int 
ParseAppleSparse(ClientData clientData, 
                 Tcl_Interp *interp, 
                 Tcl_Size argc, 
                 G3_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicAnalysisBuilder *theBuilder = static_cast<BasicAnalysisBuilder *>(clientData);

  if (argc < 2) {
    opserr << OpenSees::PromptValueError
           << "insufficient number of arguments"
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  if (strcmp(argv[1], "AppleSparse") != 0) {
    return -1;
  }

  FactorKind kind = FactorKind::General_QR; //FactorKind::SPD_Cholesky;
  int currentArg = 2;
  while (currentArg < argc) {
    if (strcmp(argv[currentArg], "-symmetric") == 0) {
      kind = FactorKind::Symmetric_LDLT;
      currentArg++;
    } else if (strcmp(argv[currentArg], "-unsymmetric") == 0) {
      kind = FactorKind::General_QR;
      currentArg++;
    } else {
      opserr << "unknown option " << argv[currentArg] << "\n";
      return TCL_ERROR;
    }
  }

  LinearSOE *theSOE = CreateAppleSparse(kind);

  theBuilder->set(theSOE);

  return TCL_OK;
}