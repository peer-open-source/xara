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
// Description: This file implements the selection of a Numberer object,
// which is used to optimally number the degrees of freedom of a problem.
//
#include <Parsing.h>
#include <Logging.h>
#include <assert.h>
#include <BasicAnalysisBuilder.h>
#include <numberer/DOF_Numberer.h>
#include <numberer/PlainNumberer.h>
#include <RCM.h>
#include <AMDNumberer.h>

#if defined(XARA_HAVE_PARALLEL_NUMBERING) // defined(_PARALLEL_PROCESSING) || defined(_PARALLEL_INTERPRETERS)
#  include <numberer/ParallelNumberer.h>
#endif

//
// command that sets the Numberer object
//
int
XaraCmd_numberer(ClientData clientData, 
                 Tcl_Interp* interp, 
                 ArgSize argc, 
                 TCL_Char ** const argv)
{
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)clientData;
  assert(builder != nullptr);


  DOF_Numberer *theNumberer = nullptr;

  // make sure at least one other argument to contain numberer
  if (argc < 2) {
    opserr << "WARNING need to specify a Numberer type \n";
    return TCL_ERROR;
  }


  // check argv[1] for type of Numberer and create the object
  if (strcmp(argv[1], "Plain") == 0) {
    theNumberer = new PlainNumberer();

  } else if (strcmp(argv[1], "RCM") == 0) {
    theNumberer = new DOF_Numberer(*new RCM(false));

  } else if (strcmp(argv[1], "AMD") == 0) {
    theNumberer = new DOF_Numberer(*new AMD());
  }

#if defined(XARA_HAVE_PARALLEL_NUMBERING)
  else if ((strcmp(argv[1], "ParallelPlain") == 0) ||
           (strcmp(argv[1], "Parallel") == 0)) {
    theNumberer = new ParallelNumberer;

  } else if (strcmp(argv[1], "ParallelRCM") == 0) {
    theNumberer = new ParallelNumberer(*new RCM(false));
  }
#endif

  else {
    opserr << "WARNING No Numberer type exists (Plain, RCM only) \n";
    return TCL_ERROR;
  }

  if (theNumberer == nullptr)
    return TCL_ERROR;

  builder->set(theNumberer);
  return TCL_OK;
}

int
XaraCmd_number(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char ** const argv)
{
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder *)clientData;
  assert(builder != nullptr);

  builder->domainChanged();

  return TCL_OK;
}
