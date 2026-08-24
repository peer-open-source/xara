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
#include <Parsing.h>
#include <Logging.h>
#include <ModelRegistry.h>
#include <Other/RosenbrockElement.h>

using namespace OpenSees;
int
XaraElem_Rosenbrock(ClientData clientData, Tcl_Interp *interp, ArgSize argc, TCL_Char ** const argv)
{
  if (argc < 4) {
    opserr << OpenSees::PromptValueError << "Missing required arguments: tag nodeTag" << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  Domain *theDomain = static_cast<ModelRegistry*>(clientData)->getDomain();

  Xara::Tag tag;
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "Invalid element tag, expected integer" << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  //
  std::array<int,1> nodes{};
  if (Tcl_GetInt(interp, argv[3], &nodes[0]) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "Invalid node tag, expected integer" << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }
  if (theDomain->getNode(nodes[0]) == nullptr) {
    opserr << OpenSees::PromptValueError << "Node with tag " << nodes[0] << " does not exist in the domain" << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  auto *theElement = new RosenbrockElement(tag, nodes);

  if (theDomain->addElement(theElement) != true) {
    opserr << OpenSees::PromptValueError << "Failed to add RosenbrockElement to domain" << OpenSees::SignalMessageEnd;
    delete theElement;
    return TCL_ERROR;
  }

  return TCL_OK;
}