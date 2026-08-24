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
// Description: This file implements commands for interacting with nodes
// in the domain.
//
// All commands assume a Domain* is passed as clientData.
//
#include <math.h>
#include <assert.h>
#include <set>
#include <vector>
#include <algorithm>
// Framework
#include <Logging.h>
#include <Parsing.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Versor.h>
#include <Domain.h>
#include <DOF_Group.h>
#include <Node.h>
#include <NodeIter.h>
#include <Pressure_Constraint.h>
#include <MP_Constraint.h>
#include <MP_ConstraintIter.h>


int
getNodeTags(ClientData clientData,
            Tcl_Interp *interp, 
            ArgSize argc,
            TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = static_cast<Domain*>(clientData);

  NodeIter &nodeIter = the_domain->getNodes();

  Tcl_Obj* result = Tcl_NewListObj(the_domain->getNumNodes(), nullptr);

  Node *node;
  while ((node = nodeIter()) != nullptr)
    Tcl_ListObjAppendElement(interp, result, Tcl_NewIntObj(node->getTag()));

  Tcl_SetObjResult(interp, result);

  return TCL_OK;
}


int
findID(ClientData clientData, Tcl_Interp *interp, ArgSize argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *theDomain = (Domain*)clientData;

  if (argc < 2) {
    opserr << "WARNING want - findNodesWithID ?id\n";
    return TCL_ERROR;
  }

  int tag;
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING findNodesWithID eleTag? dof? - could not read nodeTag? \n";
    return TCL_ERROR;
  }

  NodeIter &theNodes = theDomain->getNodes();
  Node *theNode;
  char buffer[20] = {0};

  while ((theNode = theNodes()) != nullptr) {
    DOF_Group *theGroup = theNode->getDOF_GroupPtr();
    if (theGroup != nullptr) {
      const ID &nodeID = theGroup->getID();
      for (int i = 0; i < nodeID.Size(); ++i) {
        if (nodeID(i) == tag) {
          sprintf(buffer, "%d ", theNode->getTag());
          Tcl_AppendResult(interp, buffer, NULL);
          break;
        }
      }
    }
  }

  return TCL_OK;
}

int
XaraCmd_setNodeCoord(ClientData clientData, 
             Tcl_Interp *interp, 
             ArgSize argc,
             TCL_Char ** const argv)
{
  //
  // setNodeCoord nodeTag? dim? value?
  //
  assert(clientData != nullptr);
  Domain *domain = (Domain*)clientData;
  if (argc < 4) {
    opserr << OpenSees::PromptValueError 
           << "expected setNodeCoord nodeTag? dim? value?"
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  Xara::Tag tag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError
           << "could not read nodeTag"
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  int dim;
  double value;

  if (Tcl_GetInt(interp, argv[2], &dim) != TCL_OK) {
    opserr << OpenSees::PromptValueError
           << "could not read dim"
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &value) != TCL_OK) {
    opserr << "WARNING setNodeCoord nodeTag? dim? value? - could not read "
              "value? \n";
    return TCL_ERROR;
  }

  Node *theNode = domain->getNode(tag);

  if (theNode == nullptr) {
    opserr << OpenSees::PromptValueError 
           << "Unable to find node with tag '" << tag << "'"
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  //
  // TODO: Check dimensions
  //

  Vector coords(theNode->getCrds());
  coords(dim - 1) = value;
  theNode->setCrds(coords);
  domain->domainChange();
  return TCL_OK;
}


template <NodeData Response>
int
nodeResponseTemplate(ClientData clientData, Tcl_Interp *interp, ArgSize argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *domain = static_cast<Domain*>(clientData);

  if (argc < 2) {
    opserr << OpenSees::PromptValueError
           << "Insufficient arguments"
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  int tag;
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError
           << "Failed to read nodeTag"
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  int dof = -1;
  if (argc > 2) {
    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
      opserr << OpenSees::PromptValueError
             << "Failed to read dof"
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }
  }

  dof--;

  const Vector *response = domain->getNodeResponse(tag, Response);

  if (response == nullptr) {
    opserr << OpenSees::PromptValueError
           << "Failed to find node with tag " << tag
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  Tcl_Size size = response->Size();

  if (dof >= 0) {
    if (dof >= size) {
      opserr << OpenSees::PromptValueError
             << "dofTag too large"
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }

    Tcl_SetObjResult(interp, Tcl_NewDoubleObj((*response)(dof)));
  }
  else {
    Tcl_Obj* list = Tcl_NewListObj(size, nullptr);
    for (int i = 0; i < size; ++i)
      Tcl_ListObjAppendElement(interp, list, Tcl_NewDoubleObj((*response)(i)));

    Tcl_SetObjResult(interp, list);
  }

  return TCL_OK;
}

int
XaraCmd_nodeDisp(ClientData clientData, Tcl_Interp *interp, ArgSize argc, TCL_Char ** const argv)
{
  return nodeResponseTemplate<NodeData::Disp>(clientData, interp, argc, argv);
}

int
XaraCmd_nodeVel(ClientData clientData, Tcl_Interp *interp, ArgSize argc, TCL_Char ** const argv)
{
  return nodeResponseTemplate<NodeData::Vel>(clientData, interp, argc, argv);
}

int
XaraCmd_nodeAccel(ClientData clientData, Tcl_Interp *interp, ArgSize argc, TCL_Char ** const argv)
{
  return nodeResponseTemplate<NodeData::Accel>(clientData, interp, argc, argv);
}

int
XaraCmd_nodeUnbalance(ClientData clientData, Tcl_Interp *interp, ArgSize argc, TCL_Char ** const argv)
{
  return nodeResponseTemplate<NodeData::UnbalancedLoad>(clientData, interp, argc, argv);
}

int
XaraCmd_nodeReaction(ClientData clientData, Tcl_Interp *interp, ArgSize argc, TCL_Char ** const argv)
{
  return nodeResponseTemplate<NodeData::Reaction>(clientData, interp, argc, argv);
}

int
XaraCmd_nodeMass(ClientData clientData, Tcl_Interp *interp, ArgSize argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = static_cast<Domain*>(clientData);

  if (argc < 3) {
    opserr << "WARNING want - nodeMass nodeTag? nodeDOF?\n";
    return TCL_ERROR;
  }

  int tag, dof;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "Invalid tag " << argv[1] 
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "Invalid dof " << argv[2] 
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  char buffer[40];

  Node *theNode = the_domain->getNode(tag);
  if (theNode == nullptr) {
    opserr << OpenSees::PromptValueError
           << "Node with tag " << tag << " not found"
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }
  int numDOF = theNode->getNumberDOF();
  if (dof < 1 || dof > numDOF) {
    opserr << OpenSees::PromptValueError
           << "dof " << dof << " not in range [1, " << numDOF << "]"
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }
  else {
    const Matrix &mass = theNode->getMass();
    sprintf(buffer, "%35.20f", mass(dof - 1, dof - 1));
    Tcl_AppendResult(interp, buffer, NULL);
  }

  return TCL_OK;
}

int
XaraCmd_nodePressure(ClientData clientData, Tcl_Interp *interp, ArgSize argc,
                     TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = static_cast<Domain*>(clientData);

  if (argc < 2) {
    opserr << "WARNING: want - nodePressure nodeTag?\n";
    return TCL_ERROR;
  }
  int tag;
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "Invalid tag " << argv[1] 
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }
  double pressure = 0.0;
  Pressure_Constraint *thePC = the_domain->getPressure_Constraint(tag);
  if (thePC != nullptr) {
    pressure = thePC->getPressure();
  }
  Tcl_SetObjResult(interp, Tcl_NewDoubleObj(pressure));
  return TCL_OK;
}


int
nodeBounds(ClientData clientData, Tcl_Interp *interp, ArgSize argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = static_cast<Domain*>(clientData);

  // TODO: Clean these up
  char *resDataPtr  = nullptr;
  int   resDataSize = 0;

  const int requiredDataSize = 20*6;
  if (requiredDataSize > resDataSize) {
    if (resDataPtr != 0) {
      delete[] resDataPtr;
    }
    resDataPtr = new char[requiredDataSize];
    resDataSize = requiredDataSize;
  }

  for (int i = 0; i < requiredDataSize; ++i)
    resDataPtr[i] = '\n';

  const Vector &bounds = the_domain->getPhysicalBounds();

  int cnt = 0;
  for (int j = 0; j < 6; j++) {
    cnt += sprintf(&resDataPtr[cnt], "%.6e  ", bounds(j));
  }

  Tcl_SetResult(interp, resDataPtr, TCL_STATIC);

  return TCL_OK;
}


int
XaraCmd_setNodeVel(ClientData clientData, Tcl_Interp *interp, ArgSize argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = static_cast<Domain*>(clientData);

  if (argc < 4) {
    opserr << "WARNING want - setNodeVel nodeTag? dof? value? <-commit>\n";
    return TCL_ERROR;
  }

  int tag;
  int dof = -1;
  double value = 0.0;
  bool commit = false;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING setNodeVel nodeTag? dof? value?- could not read "
              "nodeTag? \n";
    return TCL_ERROR;
  }

  Node *theNode = the_domain->getNode(tag);
  if (theNode == nullptr) {
    opserr << OpenSees::PromptValueError
           << "node with tag " << tag << " not found"
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "Invalid dof " << argv[2] 
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &value) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "Invalid value " << argv[3] 
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }
  if (argc > 4 && strcmp(argv[4], "-commit") == 0)
    commit = true;

  dof--;

  int numDOF = theNode->getNumberDOF();

  if (dof >= 0 && dof < numDOF) {
    Vector vel(numDOF);
    vel = theNode->getVel();
    vel(dof) = value;
    theNode->setTrialVel(vel);
  }
  if (commit)
    theNode->commitState();

  return TCL_OK;
}


#if 0
int
setNodeDisp(ClientData clientData, Tcl_Interp *interp, ArgSize argc,
            TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = static_cast<Domain*>(clientData);

  if (argc < 4) {
    opserr << "WARNING want - setNodeDisp nodeTag? dof? value? <-commit>\n";
    return TCL_ERROR;
  }

  int tag;
  int dof = -1;
  double value = 0.0;
  bool commit = false;
  bool increment = false;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING setNodeDisp nodeTag? dof? value?- could not read "
              "nodeTag? \n";
    return TCL_ERROR;
  }

  Node *theNode = the_domain->getNode(tag);
  if (theNode == nullptr) {
    opserr << OpenSees::PromptValueError
           << "node with tag " << tag << " not found"
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "Invalid dof " << argv[2] 
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &value) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "Invalid value " << argv[3] 
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }
  if (argc > 4 && strcmp(argv[4], "-commit") == 0)
    commit = true;
  if (argc > 4 && strcmp(argv[4], "-increment") == 0)
    increment = true;

  dof--;

  int numDOF = theNode->getNumberDOF();

  if (dof >= 0 && dof < numDOF) {
    Vector vel(numDOF);
    vel = theNode->getDisp();
    vel(dof) = value;
    if (increment)
      theNode->incrTrialDisp(vel);
    else
      theNode->setTrialDisp(vel);
  }
  if (commit)
    theNode->commitState();

  return TCL_OK;
}

#else 
int
XaraCmd_setNodeDisp(ClientData clientData, Tcl_Interp *interp, ArgSize argc,
            TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = static_cast<Domain*>(clientData);

  bool commit    = false;
  bool increment = false;

  // First pass: separate optional flags from positional arguments.
  // Flags are permitted anywhere after the command name.
  constexpr int MAX_POS = 8;
  Tcl_Size pos[MAX_POS];
  int npos = 0;
  for (Tcl_Size i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "-commit") == 0) {
      commit = true;
    } else if (strcmp(argv[i], "-increment") == 0) {
      increment = true;
    } else {
      if (npos >= MAX_POS) {
        opserr << "WARNING setNodeDisp - too many arguments\n";
        return TCL_ERROR;
      }
      pos[npos++] = i;
    }
  }

  if (npos != 2 && npos != 3) {
    opserr << "WARNING want - setNodeDisp nodeTag? nodeValues? <-commit> <-increment>\n"
           << "       or   - setNodeDisp nodeTag? dof? value? <-commit> <-increment>\n";
    return TCL_ERROR;
  }

  int tag;
  if (Tcl_GetInt(interp, argv[pos[0]], &tag) != TCL_OK) {
    opserr << "WARNING setNodeDisp - could not read nodeTag\n";
    return TCL_ERROR;
  }

  Node *theNode = the_domain->getNode(tag);
  if (theNode == nullptr) {
    opserr << OpenSees::PromptValueError
           << "node with tag " << tag << " not found"
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  int numDOF = theNode->getNumberDOF();

  if (npos == 2) {
    // List form:  setNodeDisp nodeTag {v1 v2 ... vN}
    Tcl_Size listLen;
    TCL_Char **listArgv;
    if (Tcl_SplitList(interp, argv[pos[1]], &listLen, &listArgv) != TCL_OK) {
      opserr << OpenSees::PromptValueError
             << "could not parse nodeValues list"
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }
    if ((int)listLen != numDOF) {
      opserr << OpenSees::PromptValueError
             << "nodeValues has " << (int)listLen
             << " entries, expected " << numDOF
             << OpenSees::SignalMessageEnd;
      Tcl_Free((char*)listArgv);
      return TCL_ERROR;
    }

    Vector vel(numDOF);
    for (Tcl_Size i = 0; i < listLen; ++i) {
      double v;
      if (Tcl_GetDouble(interp, listArgv[i], &v) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "invalid value '" << listArgv[i] << "' in nodeValues"
               << OpenSees::SignalMessageEnd;
        Tcl_Free((char*)listArgv);
        return TCL_ERROR;
      }
      vel((int)i) = v;
    }
    Tcl_Free((char*)listArgv);

    if (increment)
      theNode->incrTrialDisp(vel);
    else
      theNode->setTrialDisp(vel);

  } else {
    // Single-DOF form:  setNodeDisp nodeTag dof value
    int    dof;
    double value;
    if (Tcl_GetInt(interp, argv[pos[1]], &dof) != TCL_OK) {
      opserr << OpenSees::PromptValueError
             << "Invalid dof " << argv[pos[1]]
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[pos[2]], &value) != TCL_OK) {
      opserr << OpenSees::PromptValueError
             << "Invalid value " << argv[pos[2]]
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }
    dof--;

    if (dof >= 0 && dof < numDOF) {
      Vector vel(numDOF);
      vel = theNode->getDisp();
      vel(dof) = value;
      if (increment)
        theNode->incrTrialDisp(vel);
      else
        theNode->setTrialDisp(vel);
    }
  }

  if (commit)
    theNode->commitState();

  return TCL_OK;
}
#endif 


int
XaraCmd_setNodeAccel(ClientData clientData, Tcl_Interp *interp, ArgSize argc,
             TCL_Char ** const argv)
{
  assert(clientData != nullptr);

  Domain *the_domain = (Domain *)clientData;

  if (argc < 4) {
    opserr << "WARNING want - setNodeAccel nodeTag? dof? value? <-commit>\n";
    return TCL_ERROR;
  }

  int tag;
  int dof = -1;
  double value = 0.0;
  bool commit = false;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING setNodeAccel nodeTag? dof? value?- could not read "
              "nodeTag? \n";
    return TCL_ERROR;
  }

  Node *theNode = the_domain->getNode(tag);
  if (theNode == nullptr) {
    opserr << OpenSees::PromptValueError
           << "node with tag " << tag << " not found"
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "Invalid dof " << argv[2] 
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &value) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "Invalid value " << argv[3] 
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  if (argc > 4 && strcmp(argv[4], "-commit") == 0)
    commit = true;

  dof--;

  int numDOF = theNode->getNumberDOF();

  if (dof >= 0 && dof < numDOF) {
    Vector vel(numDOF);
    vel = theNode->getAccel();
    vel(dof) = value;
    theNode->setTrialAccel(vel);
  }
  if (commit)
    theNode->commitState();

  return TCL_OK;
}

int 
XaraCmd_setNodePressure(ClientData clientData, Tcl_Interp *interp, ArgSize argc,
              TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = static_cast<Domain*>(clientData);

  if (argc < 3) {
    opserr << "WARNING want - setNodePressure nodeTag? value?\n";
    return TCL_ERROR;
  }

  int tag;
  double value;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError
           << "Invalid tag " << argv[1] 
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[2], &value) != TCL_OK) {
    opserr << OpenSees::PromptValueError
           << "Invalid value " << argv[2] 
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  Pressure_Constraint *thePC = the_domain->getPressure_Constraint(tag);
  if (thePC == nullptr) {
    opserr << OpenSees::PromptValueError
           << "Pressure constraint with tag " << tag << " not found"
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  thePC->setPressure(value);

  return TCL_OK;
}


int
XaraCmd_nodeRotation(ClientData clientData, Tcl_Interp *interp, ArgSize argc,
             TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = static_cast<Domain*>(clientData);

  if (argc < 2) {
    opserr << OpenSees::PromptValueError << "Missing required argument: tag" << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  int tag;
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError
           << "Invalid tag " << argv[1] 
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  // TODO: This may need some adjusting for PartitionedDomains and parallel
  Node *theNode = the_domain->getNode(tag);
  if (theNode == nullptr)
    return TCL_ERROR;

  const Versor rotation = theNode->getCommitRotation();

  Tcl_Obj* list = Tcl_NewListObj(4, nullptr);
  for (int i = 0; i < 3; ++i)
    Tcl_ListObjAppendElement(interp, list, Tcl_NewDoubleObj(rotation.vector[i]));
  Tcl_ListObjAppendElement(interp, list, Tcl_NewDoubleObj(rotation.scalar));

  Tcl_SetObjResult(interp, list);

  return TCL_OK;
}

int
XaraCmd_nodeResponse(ClientData clientData, Tcl_Interp *interp, ArgSize argc,
             TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain *)clientData;

  if (argc < 4) {
    opserr << "WARNING want - nodeResponse nodeTag? dof? responseID?\n";
    return TCL_ERROR;
  }

  int tag, dof, responseID;
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError
           << "Invalid tag " << argv[1] 
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
    opserr << OpenSees::PromptValueError
           << "Invalid dof " << argv[2] 
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[3], &responseID) != TCL_OK) {
    if (strcmp(argv[3], "displacement") == 0)
      responseID = (int)NodeData::Disp;
    else if (strcmp(argv[3], "velocity") == 0)
      responseID = (int)NodeData::Vel;
    else if (strcmp(argv[3], "acceleration") == 0)
      responseID = (int)NodeData::Accel;
    else if (strcmp(argv[3], "residual") == 0)
      responseID = (int)NodeData::UnbalancedLoad;
    else if (strcmp(argv[3], "reactionForce") == 0)
      responseID = (int)NodeData::Reaction;
    else {
      opserr << "WARNING unknown response " << argv[3] << "\n";
      return TCL_ERROR;
    }
  }

  dof--;

  const Vector *nodalResponse = nullptr;
  if (false) {
    // This would need some adjusting for PartitionedDomain
    Node *theNode = the_domain->getNode(tag);
    if (theNode == nullptr)
      return TCL_ERROR;

    nodalResponse = theNode->getResponse((NodeData)responseID);

  } else
    nodalResponse = the_domain->getNodeResponse(tag, (NodeData)responseID);


  if (nodalResponse == nullptr || nodalResponse->Size() < dof || dof < 0)
    // TODO: add error message
    return TCL_ERROR;

  Tcl_SetObjResult(interp, Tcl_NewDoubleObj((*nodalResponse)(dof)));

  return TCL_OK;
}

int
XaraCmd_nodeEigenvector(ClientData clientData, Tcl_Interp *interp, ArgSize argc,
                TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *domain = (Domain*)clientData;

  if (argc < 3) {
    opserr << "WARNING want - nodeEigenVector nodeTag? eigenVector? <dof?>\n";
    return TCL_ERROR;
  }

  int tag;
  int eigenvector = 0;
  int dof = -1;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr
        << "WARNING nodeEigenvector nodeTag? dof? - could not read nodeTag? \n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2], &eigenvector) != TCL_OK) {
    opserr << "WARNING nodeEigenvector nodeTag? dof? - could not read dof? \n";
    return TCL_ERROR;
  }

  if (argc > 3) {
    if (Tcl_GetInt(interp, argv[3], &dof) != TCL_OK) {
      opserr
          << "WARNING nodeEigenvector nodeTag? dof? - could not read dof? \n";
      return TCL_ERROR;
    }
  }

  dof--;
  eigenvector--;
  Node *theNode = domain->getNode(tag);
  const Matrix &theEigenvectors = theNode->getEigenvectors();

  int size = theEigenvectors.noRows();
  int numEigen = theEigenvectors.noCols();

  if (eigenvector < 0 || eigenvector >= numEigen) {
    opserr << "WARNING nodeEigenvector nodeTag? dof? - eigenvecor too large\n";
    return TCL_ERROR;
  }

  if (dof >= 0) {
    if (dof >= size) {
      opserr << "WARNING nodeEigenvector nodeTag? dof? - dofTag? too large\n";
      return TCL_ERROR;
    }

    double value = theEigenvectors(dof, eigenvector);

    // copy the value to the Tcl string that is returned
    char buffer[40];
    sprintf(buffer, "%35.20f", value);
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  } else {

    char buffer[40];
    for (int i = 0; i < size; ++i) {
      double value = theEigenvectors(i, eigenvector);
      sprintf(buffer, "%35.20f", value);
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }

  return TCL_OK;
}


int
XaraCmd_reactions(ClientData clientData, 
                  Tcl_Interp *interp, ArgSize argc,
                  TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *domain = static_cast<Domain *>(clientData);

  int incInertia = 0;

  if (argc == 2) {
    if ((strcmp(argv[1], "-incInertia") == 0) ||
        (strcmp(argv[1], "-dynamical") == 0) ||
        (strcmp(argv[1], "-Dynamic") == 0) ||
        (strcmp(argv[1], "-dynamic") == 0))

      incInertia = 1;

    else if ((strcmp(argv[1], "-rayleigh") == 0))

      incInertia = 2;
  }

  domain->calculateNodalReactions(incInertia);

  return TCL_OK;
}


int
nodeCoord(ClientData clientData, Tcl_Interp *interp, ArgSize argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = static_cast<Domain*>(clientData);

  if (argc < 2) {
    opserr << OpenSees::PromptValueError << "Missing required argument: nodeTag" << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  Xara::Tag tag;
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "Invalid node tag, expected integer" << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  int dim = -1;

  if (argc > 2) {
    if (strcmp(argv[2], "X") == 0 || strcmp(argv[2], "x") == 0 ||
        strcmp(argv[2], "1") == 0)
      dim = 0;
    else if (strcmp(argv[2], "Y") == 0 || strcmp(argv[2], "y") == 0 ||
             strcmp(argv[2], "2") == 0)
      dim = 1;
    else if (strcmp(argv[2], "Z") == 0 || strcmp(argv[2], "z") == 0 ||
             strcmp(argv[2], "3") == 0)
      dim = 2;
    else {
      opserr << OpenSees::PromptValueError 
             << "could not read dim"
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }
  }

  Node *theNode = the_domain->getNode(tag);

  if (theNode == nullptr) {
    opserr << OpenSees::PromptValueError 
           << "Unable to find node with tag '" << tag << "'"
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  const Vector &coords = theNode->getCrds();

  char buffer[40];
  int size = coords.Size();
  if (dim == -1) {
    for (int i = 0; i < size; ++i) {
      sprintf(buffer, "%35.20f", coords(i));
      Tcl_AppendResult(interp, buffer, NULL);
    }
    return TCL_OK;

  } else if (dim < size) {
    double value = coords(dim);
    Tcl_SetObjResult(interp, Tcl_NewDoubleObj(value));
    return TCL_OK;
  }

  return TCL_ERROR;
}


int
retainedNodes(ClientData clientData, Tcl_Interp *interp, ArgSize argc,
              TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *domain = (Domain*)clientData;
  bool all = 1;
  int cNode;
  if (argc > 1) {
    if (Tcl_GetInt(interp, argv[1], &cNode) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "Invalid cNode \n";
      return TCL_ERROR;
    }
    all = 0;
  }

  MP_Constraint *theMP;
  MP_ConstraintIter &mpIter = domain->getMPs();

  // get unique constrained nodes with set
  std::set<int> tags;
  int tag;
  while ((theMP = mpIter()) != nullptr) {
    tag = theMP->getNodeRetained();
    if (all || cNode == theMP->getNodeConstrained()) {
      tags.insert(tag);
    }
  }

  // assign set to vector and sort
  std::vector<int> tagv;
  tagv.assign(tags.begin(), tags.end());
  std::sort(tagv.begin(), tagv.end());

  // loop through unique, sorted tags, adding to output
  char buffer[20];
  for (int tag : tagv) {
    sprintf(buffer, "%d ", tag);
    Tcl_AppendResult(interp, buffer, NULL);
  }

  return TCL_OK;
}



int
nodeDOFs(ClientData clientData, Tcl_Interp *interp, ArgSize argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = static_cast<Domain*>(clientData);

  if (argc < 2) {
    opserr << OpenSees::PromptValueError 
           << "Missing required arguments"
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  int tag;
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "Failed to read nodeTag " << argv[1]
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }


  Node *theNode = the_domain->getNode(tag);
  if (theNode == nullptr) {
    opserr << OpenSees::PromptValueError 
           << "node with tag " << tag << " not found" 
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  int numDOF = theNode->getNumberDOF();
  DOF_Group *theDOFgroup = theNode->getDOF_GroupPtr();
  if (theDOFgroup == nullptr) {
    opserr << OpenSees::PromptValueError
           << "nodeDOFs DOF group null" 
           << OpenSees::SignalMessageEnd;
    return -1;
  }

  char buffer[40];
  const ID &eqnNumbers = theDOFgroup->getID();
  for (int i = 0; i < numDOF; ++i) {
    sprintf(buffer, "%d ", eqnNumbers(i));
    Tcl_AppendResult(interp, buffer, NULL);
  }

  return TCL_OK;
}

