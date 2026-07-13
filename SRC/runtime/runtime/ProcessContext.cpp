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
#include "ProcessContext.h"
#include <TclPackageClassBroker.h>
#include <Channel.h>
#include <Logging.h>
#include <runtime/interpreter/Interpreter.h>

using namespace OpenSees;

static int getPID(ClientData,  Tcl_Interp *, ArgSize, TCL_Char ** const argv);
static int getNP(ClientData,   Tcl_Interp *, ArgSize, TCL_Char ** const argv);


ProcessContext::~ProcessContext() 
{
  delete m_obroker;
}


ProcessContext::ProcessContext()
: m_obroker(new TclPackageClassBroker())
, theMachine(m_obroker, 0, nullptr)
{
  int OPS_rank = theMachine.getPID();
  int OPS_np = theMachine.getNP();

  if (OPS_rank == 0) {
    opserr << "Rank = " << OPS_rank << " / " << OPS_np << "\n";
    theChannels = new Channel *[OPS_np-1];
    numChannels = OPS_np-1;


    for (int j=0; j<OPS_np-1; j++) {
      Channel *otherChannel = theMachine.getRemoteProcess();
      theChannels[j] = otherChannel;
    //   otherChannel->sendID(0,0,data);
    //   otherChannel->sendMsg(0,0,msgChar);
    }
  }
  else {
    theChannels = new Channel *[1];
    numChannels = 1;

    Channel *myChannel = theMachine.getMyChannel();
    theChannels[0] = myChannel;
  }
}

int
ProcessContext::setup(Interpreter& interp)
{
  Tcl_CreateCommand(&interp, "getNP",     &getNP,   (ClientData)this, nullptr);
  Tcl_CreateCommand(&interp, "getPID",    &getPID,  (ClientData)this, nullptr);
  return 0;
}


static int
getPID(ClientData clientData, Tcl_Interp *interp, ArgSize argc, TCL_Char ** const argv)
{
  int pid = 0;

  // MPI_Comm_rank(MPI_COMM_WORLD, &pid);

  ProcessContext* context = (ProcessContext*)clientData;

  if (context != nullptr)
    pid = context->getProcessID();

  // set the returned integer
  Tcl_SetObjResult(interp, Tcl_NewIntObj(pid));

  return TCL_OK;
}


static int
getNP(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  int np = 1;
  // MPI_Comm_size(MPI_COMM_WORLD, &np);

  ProcessContext* context = (ProcessContext*)clientData;

  if (context != nullptr)
    np = context->getProcessCount();

  // set the returned integer
  Tcl_SetObjResult(interp, Tcl_NewIntObj(np));

  return TCL_OK;
}
