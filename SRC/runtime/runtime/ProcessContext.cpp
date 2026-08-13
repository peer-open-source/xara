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
#include <XaraClassBroker.h>
#include <Channel.h>
#include <Logging.h>
#include <Parsing.h>
// #include <Interpreter.h>

#if defined(XARA_ENABLE_MPI)
# include <mpi.h>
#endif

using namespace Xara;

static Tcl_CmdProc getPID;
static Tcl_CmdProc getNP;
static Tcl_CmdProc opsBarrier;


#if !defined(XARA_ENABLE_MPI)
ProcessContext::ProcessContext()
: m_obroker(new XaraClassBroker())
{

}
#else
ProcessContext::ProcessContext()
: m_obroker(new XaraClassBroker())
, theMachine(m_obroker, 0, nullptr)
{
  int OPS_rank = theMachine.getPID();
  int OPS_np   = theMachine.getNP();

  if (OPS_rank == 0) {
    opsdbg << "Rank = " << OPS_rank << " / " << OPS_np << "\n";

    theChannels = new Channel *[OPS_np-1];
    numChannels = OPS_np-1;


    for (int j=0; j<OPS_np-1; j++) {
      Channel *otherChannel = theMachine.getRemoteProcess();
      theChannels[j] = otherChannel;
    }
  }
  else {
    theChannels = new Channel *[1];
    numChannels = 1;

    Channel *myChannel = theMachine.getMyChannel();
    theChannels[0] = myChannel;
  }
}
#endif

ProcessContext::~ProcessContext() 
{
  // theMachine.shutdown();
  // MPI_Finalize();

  delete m_obroker;
}


int
ProcessContext::setup(Tcl_Interp* interp)
{
  Tcl_CreateCommand(interp, "getNP",     &getNP,   (ClientData)this, nullptr);
  Tcl_CreateCommand(interp, "getPID",    &getPID,  (ClientData)this, nullptr);
  Tcl_CreateCommand(interp, "barrier",   &opsBarrier,  (ClientData)this, nullptr);
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


static int
opsBarrier(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
#if defined(XARA_ENABLE_MPI)
  if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS) {
    opserr << OpenSees::PromptValueError 
           << "MPI_Barrier failed"
            << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }
#endif
  return TCL_OK;
}
