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
// Description: This file implements the setResponse and getResponse commands 
// for the interpreter. These commands allow a Response object to 
// be reused in the interpreter, as opposed to the eleResponse command 
// which allocates a new Response object each time it is called. 
// This can be much more efficient, but requires extreme care at the interpreter
// level which deals with raw pointers to data of the response's Information object.
//
// The lifetime of the allocated Response is managed by the ModelRegistry,
// which is deleted when ModelRegistry is deleted.
//
/// Created: 06/2026
//
#include <Element.h>
#include <Domain.h>
#include <Parsing.h>
#include <Logging.h>
#include <string>
#include <ModelRegistry.h>
#include <DummyStream.h>
#include <recorder/response/Response.h>
#include <runtime/runtime/InterpreterResponse.h>
#include <cstdint>   // uintptr_t


static Response*
setElementResponse(Domain& domain, int eleTag, const char **argv, Tcl_Size argc)
{
  Element *theEle = domain.getElement(eleTag);

  if (theEle == nullptr)
    return nullptr;

  DummyStream dummy;
  return theEle->setResponse(argv, argc, dummy);
}


int 
XaraCmd_setResponse(ClientData context, 
                    Tcl_Interp *interp, 
                    ArgSize argc,
                    TCL_Char** const argv)
{
  assert(context != nullptr);
  ModelRegistry *theRegistry = static_cast<ModelRegistry*>(context);
  Domain *theDomain = theRegistry->getDomain();

  if (argc < 3) {
    opserr << OpenSees::PromptValueError << "want - setResponse eleTag? args...\n";
    return TCL_ERROR;
  }

  int tag;
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "Invalid element tag: " << argv[1]
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  Response *theResponse = setElementResponse(*theDomain, tag, argv + 2, argc - 2);
  if (theResponse == nullptr) {
    opserr << OpenSees::PromptValueError 
           << "Failed to set response for element " << tag << "\n";
    return TCL_ERROR;
  }

  if (theResponse->getResponse() < 0) {
    opserr << OpenSees::PromptValueError 
           << "Failed to get response" 
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }


  Information &info = theResponse->getInformation();

  InterpreterResponse* iresponse = new InterpreterResponse(theResponse, 
                                                           InterpreterResponse::Type::Vector, 
                                                           info.theVector.Size(), 
                                                           &info.theVector(0));

  int id = theRegistry->addResponse(iresponse);


  uintptr_t addr = reinterpret_cast<uintptr_t>(&info.theVector(0));

  Tcl_AppendResult(interp,
                   std::to_string(id).c_str(),                       " ",
                   std::to_string(info.theVector.Size()).c_str(),    " ",
                   std::to_string(addr).c_str(),
                   (char*)nullptr);
  return TCL_OK;
}


int
XaraCmd_getResponse(ClientData context, 
                    Tcl_Interp *interp,
                    ArgSize argc,
                    TCL_Char** const argv)
{
  assert(context != nullptr);
  ModelRegistry *theRegistry = static_cast<ModelRegistry*>(context);

  if (argc < 2) {
    opserr << OpenSees::PromptValueError << "want - getResponse responseID? \n";
    return TCL_ERROR;
  }

  int id;
  if (Tcl_GetInt(interp, argv[1], &id) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "Invalid responseID: " << argv[1]
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  //
  // Parsing complete.
  // Now get the response from the registry, update the data, and 
  // check that the data references created by setResponse are still valid.
  //
  InterpreterResponse *theResponse = theRegistry->getResponse(id);
  if (theResponse == nullptr) {
    opserr << OpenSees::PromptValueError 
           << "could not find response with ID " << id 
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  if (theResponse->response->getResponse() < 0) {
    opserr << OpenSees::PromptValueError 
           << "Failed to get response for ID " << id 
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  // Check the integrity of the response data
  Information &info = theResponse->response->getInformation();
  if (info.theVector.Size() != theResponse->size) {
    opserr << OpenSees::PromptValueError 
           << "response size mismatch" 
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  if (&info.theVector(0) != theResponse->data) {
    opserr << OpenSees::PromptValueError 
           << "response data mismatch" 
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  return TCL_OK;
}
