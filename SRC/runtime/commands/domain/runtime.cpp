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
// These commands expect a Domain* as their clientData.
// 
// Written: cmp
//
#include <assert.h>
#include <Domain.h>
#include <Parsing.h>
#include <Logging.h>
#include <LoadPattern.h>

int
XaraCmd_setLoadConst(ClientData clientData, Tcl_Interp *interp, ArgSize argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain* domain = (Domain*)clientData;

  bool const_done = false;

  int argi = 1;
  while (argi < argc) {
    if (strcmp(argv[argi], "-time") == 0) {
      argi++;
      if (argi == argc) {
        opserr << OpenSees::PromptValueError << "no time value after -time flag\n";
        return TCL_ERROR;
      }
      double newTime;
      if (Tcl_GetDouble(interp, argv[argi], &newTime) != TCL_OK) {
        opserr << "WARNING readingvalue - loadConst -time value \n";
        return TCL_ERROR;
      }
      argi++;
      domain->setCurrentTime(newTime);
      domain->setCommittedTime(newTime);
    }
    else if (strcmp(argv[argi], "-pattern") == 0) {
      argi++;
      if (argi == argc) {
        opserr << OpenSees::PromptValueError << "no load pattern tag after -pattern flag\n";
        return TCL_ERROR;
      }
      int loadPatternTag;
      if (Tcl_GetInt(interp, argv[argi], &loadPatternTag) != TCL_OK) {
        opserr << "WARNING reading value - loadConst -pattern tag \n";
        return TCL_ERROR;
      }
      argi++;
      LoadPattern* thePattern = domain->getLoadPattern(loadPatternTag);
      if (thePattern == nullptr) {
        opserr << "WARNING no load pattern with tag " << loadPatternTag << "\n";
        return TCL_ERROR;
      }
      thePattern->setLoadConstant();
      const_done = true;
    }
    else {
      opserr << OpenSees::PromptValueError << "invalid argument: " << argv[argi] << "\n";
      return TCL_ERROR;
    }
  }

  if (const_done == false)
    domain->setLoadConstant();

  return TCL_OK;
}


int
XaraCmd_setTime(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain* domain = (Domain*)clientData;

  if (argc < 2) {
    opserr << OpenSees::PromptValueError 
           << "missing required argument: time\n";
    return TCL_ERROR;
  }
  double newTime;
  if (Tcl_GetDouble(interp, argv[1], &newTime) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "invalid time value: " << argv[1] << "\n";
    return TCL_ERROR;
  }

  domain->setCurrentTime(newTime);
  domain->setCommittedTime(newTime);
  return TCL_OK;
}


int
XaraCmd_getTime(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain* domain = static_cast<Domain*>(clientData);

  Tcl_SetObjResult(interp, Tcl_NewDoubleObj(domain->getCurrentTime()));
  return TCL_OK;
}



int
XaraCmd_setCreep(ClientData clientData, Tcl_Interp *interp, ArgSize argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain* the_domain = (Domain*)clientData;

  if (argc < 2) {
    opserr << OpenSees::PromptValueError 
           << "missing required argument: creep value\n";
    return TCL_ERROR;
  }
  int newFlag;
  if (Tcl_GetInt(interp, argv[1], &newFlag) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "invalid creep value: " << argv[1] << "\n";
    return TCL_ERROR;
  }
  the_domain->setCreep(newFlag);

  return TCL_OK;
}
