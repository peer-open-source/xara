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
//
// Description: This file implements commands that configure Node objects
// for an analysis.
//
// Author: cmp
//
#include <string>
#include <assert.h>
#include <string.h>
#include <tcl.h>
#include <Logging.h>
#include <Parsing.h>
#include <Node.h>
#include <NodeND.h>
#include <Matrix.h>
#include <Domain.h>
#include <Parameter.h>
#include <ModelRegistry.h>

#define HeapNode Node

#define G3_NUM_DOF_BUFFER 20

using namespace Xara;


int
TclCommand_wipeNodes(ClientData context, Tcl_Interp *interp,
                    Tcl_Size argc, TCL_Char ** const argv)
{
  // TODO: Check that all nodes are deleted from the domain
  // assert(context != nullptr);
  Node::resetGlobalMatrices();
  return TCL_OK;
}


int
XaraCmd_node(ClientData context, 
                   Tcl_Interp *interp,
                   ArgSize argc,
                   TCL_Char ** const argv)
{
  assert(context != nullptr);

  ModelRegistry *builder = static_cast<ModelRegistry*>(context);

  Domain *theTclDomain = builder->getDomain();

  int ndm = builder->getNDM();
  int ndf = builder->getNDF();

  // make sure corect number of arguments on command line
  if (argc < 2 + 1) { // ndm) {
    opserr << OpenSees::PromptValueError 
           << "insufficient arguments"
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  Node *theNode = nullptr;

  // read the node id
  Xara::Tag tag;
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "invalid nodeTag " << argv[1]
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }


  Parameter* coord_params[3] = {nullptr, nullptr, nullptr};
  using namespace OpenSees::Parsing;

  // read in the coordinates and create the node
  double xLoc=0, yLoc=0, zLoc=0;
  if (ndm >= 1 && argc >= 3) {
    // create a node in 1d space
    if (GetDoubleParam(interp, *theTclDomain, argv[2], &xLoc, coord_params[0]) != TCL_OK) {
      opserr << OpenSees::PromptValueError 
             << "invalid coordinate " << argv[2] 
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }
  }

  if (ndm >= 2 && argc >= 4) {
    if (GetDoubleParam(interp, *theTclDomain, argv[3], &yLoc, coord_params[1]) != TCL_OK) {
    // if (Tcl_GetDouble(interp, argv[3], &yLoc) != TCL_OK) {
      opserr << OpenSees::PromptValueError 
             << "invalid 2nd coordinate " << argv[3]
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }
  }

  if (ndm >= 3 && argc >= 5) {
    if (GetDoubleParam(interp, *theTclDomain, argv[4], &zLoc, coord_params[2]) != TCL_OK) {
      opserr << OpenSees::PromptValueError 
             << "invalid 3rd coordinate " << argv[4]
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }
  }

  // check for -ndf override option
  int currentArg = 2 + ndm;
  if (currentArg < argc && strcmp(argv[currentArg], "-ndf") == 0) {
    if (Tcl_GetInt(interp, argv[currentArg + 1], &ndf) != TCL_OK) {
      opserr << OpenSees::PromptValueError 
             << "invalid nodal ndf given for node " << tag 
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }
    currentArg += 2;
  }

  //
  // create the node
  //
  switch (ndm) {
  case 1:
    theNode = new HeapNode(tag, ndf, xLoc);
    break;
  case 2:
    theNode = new HeapNode(tag, ndf, xLoc, yLoc);
    break;
  case 3:
#if 0
    if (getenv("NODE")) {
      switch (ndf) {
        case 3:
          theNode = new NodeND<3, 3>(tag, xLoc, yLoc, zLoc);
          break;
        case 6:
          theNode = new NodeND<3, 6>(tag, xLoc, yLoc, zLoc);
          break;
        default:
          theNode = new HeapNode(tag, ndf, xLoc, yLoc, zLoc);
          break;
      }
    } else
#endif
      theNode = new HeapNode(tag, ndf, xLoc, yLoc, zLoc, builder->getRotationType());
    break;
  }

  while (currentArg < argc) {
    if (strcmp(argv[currentArg], "-mass") == 0) {
      currentArg++;
      if (argc < currentArg + ndf) {
        opserr << OpenSees::PromptValueError 
               << "incorrect number of nodal mass terms for "
               << "node: " << tag 
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }

      double theMass;
      Matrix mass(ndf, ndf);
      for (int i = 0; i < ndf; ++i) {
        if (Tcl_GetDouble(interp, argv[currentArg++], &theMass) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid nodal mass term";
          opserr << " at dof " << i + 1 << "\n";
          return TCL_ERROR;
        }
        mass(i, i) = theMass;
      }
      theNode->setMass(mass);
    } 
    else if (strcmp(argv[currentArg], "-dispLoc") == 0) {
      
      opswrn << G3_WARN_PROMPT
             << "-dispLoc option is no longer supported\n"
             << OpenSees::SignalMessageEnd;
      currentArg += 1+ndm;
    } 
    else if (strcmp(argv[currentArg], "-disp") == 0) {
      currentArg++;
      if (argc < currentArg + ndf) {
        opserr << OpenSees::PromptValueError << "incorrect number of nodal disp terms\n";
        opserr << "node: " << tag << "\n";
        return TCL_ERROR;
      }
      Vector disp(ndf);
      double theDisp;
      for (int i = 0; i < ndf; ++i) {
        if (Tcl_GetDouble(interp, argv[currentArg++], &theDisp) != TCL_OK) {
          opserr << OpenSees::PromptValueError 
                 << "invalid nodal disp term\n"
                 << "node: " << tag << ", dof: " << i + 1 
                 << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        }
        disp(i) = theDisp;
      }
      theNode->setTrialDisp(disp);
      theNode->commitState();
    } 
    else if (strcmp(argv[currentArg], "-vel") == 0) {
      currentArg++;
      if (argc < currentArg + ndf) {
        opserr << OpenSees::PromptValueError 
               << "incorrect number of nodal vel terms, "
               << "expected " << ndf 
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }

      double theDisp;
      Vector disp(ndf);
      for (int i = 0; i < ndf; ++i) {
        if (Tcl_GetDouble(interp, argv[currentArg++], &theDisp) != TCL_OK) {
          opserr << OpenSees::PromptValueError 
                 << "invalid nodal vel term at "
                 << " dof " << i + 1 
                 << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        }
        disp(i) = theDisp;
      }
      theNode->setTrialVel(disp);
      theNode->commitState();

    } else
      currentArg++;
  }


  //
  // Setup parameters for coordinates
  //
  for (int i=0; i<3; ++i) {
    if (coord_params[i] == nullptr)
      continue;
    char index[20];
    snprintf(index, sizeof(index), "%d", i + 1);
    std::string idx = std::to_string(i + 1);
    const char* args[2] = { "coord",  index }; 
    coord_params[i]->addComponent(theNode, args, 2);
  }

  //
  // add the node to the domain
  //
  if (theTclDomain->addNode(theNode) == false) {
    opserr << OpenSees::PromptValueError << "failed to add node to the domain\n";
    delete theNode;
    return TCL_ERROR;
  }

  return TCL_OK;
}


int
TclCommand_addNodalMass(ClientData context, 
                        Tcl_Interp *interp, 
                        ArgSize argc,
                        TCL_Char ** const argv)
{
  assert(context != nullptr);
  ModelRegistry *builder = static_cast<ModelRegistry*>(context);
  Domain *theTclDomain = builder->getDomain();


  // get the id of the node
  Xara::Tag nodeId;
  if (argc < 2) {
    opserr << OpenSees::PromptValueError
           << "Missing required argument tag"
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "invalid node tag: " << argv[1]
           << "\n";
    return TCL_ERROR;
  }

  Node* node = theTclDomain->getNode(nodeId);
  if (node == nullptr) {
    opserr << OpenSees::PromptValueError 
           << "node " << nodeId 
           << " does not exist in domain\n";
    return TCL_ERROR;
  }

  const int ndf = node->getNumberDOF();
  if (argc < 2 + 1) {
    opserr << OpenSees::PromptValueError 
           << "insufficient arguments, expected "  
           << ndf 
           << " mass values\n"; 
    return TCL_ERROR;
  }

  const int ndm = node->getCrds().Size();
  // check for mass terms
  bool equal_position = true;

  Matrix mass(ndf,ndf);
  double position_inertia[3];
  int argi=2;
  int n_mass = 0;
  while (argi < argc) {
    if (strcmp(argv[argi], "-position") == 0) {
      if (n_mass > 0) {
        opserr << OpenSees::PromptValueError 
               << "cannot specify -position option after mass terms\n";
        return TCL_ERROR;
      }
      equal_position = true;
      if (argi + 1 >= argc) {
        opserr << OpenSees::PromptValueError 
               << "missing nodal mass term after -position\n"
               << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[argi+1], &position_inertia[0]) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid nodal mass term\n";
          return TCL_ERROR;
      }
      mass(0,0) = position_inertia[0];
      for (int j=1; j<ndm; ++j) {
        position_inertia[j] = position_inertia[0];
        mass(j,j) = position_inertia[j];
      }
      argi += 2;
    } 
    else {
      double theMass;
      if (Tcl_GetDouble(interp, argv[argi], &theMass) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
                << "invalid nodal mass term "
                << argv[argi]
                << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      mass(n_mass,n_mass) = theMass;
      if (n_mass < ndm)
        position_inertia[n_mass] = theMass;
      n_mass++;
      argi++;
    }
  }

  node->setMass(mass);

  return TCL_OK;
}



int
TclCommand_getNDM(ClientData context, Tcl_Interp *interp, ArgSize argc, TCL_Char ** const argv)
{
  assert(context != nullptr);
  ModelRegistry* builder = static_cast<ModelRegistry*>(context);
  Domain *the_domain = builder->getDomain();

  int ndm;

  if (argc > 1) {
    int tag;
    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "ndm nodeTag? \n";
      return TCL_ERROR;
    }

    Node *theNode = the_domain->getNode(tag);
    if (theNode == nullptr) {
      opserr << OpenSees::PromptValueError 
             << "node with tag " << tag << " does not exist "
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }
    const Vector &coords = theNode->getCrds();
    ndm = coords.Size();

  } else {
    ndm = builder->getNDM();
  }

  Tcl_SetObjResult(interp, Tcl_NewIntObj(ndm));
  return TCL_OK;
}

int
TclCommand_getNDF(ClientData context, Tcl_Interp *interp, ArgSize argc, TCL_Char ** const argv)
{
  assert(context != nullptr);
  ModelRegistry* builder = static_cast<ModelRegistry*>(context);
  Domain *the_domain = builder->getDomain();
  int ndf;

  if (argc > 1) {
    int tag;
    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "ndf nodeTag? \n";
      return TCL_ERROR;
    }
    Node *theNode = the_domain->getNode(tag);
    if (theNode == nullptr) {
      opserr << OpenSees::PromptValueError 
             << "node with tag " << tag << " does not exist \n";
      return TCL_ERROR;
    }
    ndf = theNode->getNumberDOF();

  } else {
    ndf = builder->getNDF();
  }

  Tcl_SetObjResult(interp, Tcl_NewIntObj(ndf));
  return TCL_OK;
}
