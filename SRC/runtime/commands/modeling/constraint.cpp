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
#include <assert.h>
#include <Matrix.h>
#include <ID.h>
#include <Node.h>
#include <Domain.h>
#include <Parsing.h>
#include <Logging.h>
#include <ModelRegistry.h>
#include <vector>
#include <GroupSO3.h>

#include <runtimeAPI.h>
#include <Vector3D.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <MP_Constraint.h>

#include <LoadPattern.h>

#include <ImposedMotionSP.h>
#include <ImposedMotionSP1.h>
#include <MultiSupportPattern.h>


int
TclCommand_addHomogeneousBC(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc,
                            TCL_Char ** const argv)
{
  // fix tag <fixities>
  assert(clientData != nullptr);
  ModelRegistry *builder = static_cast<ModelRegistry*>(clientData);

  Domain *theTclDomain = builder->getDomain();

  if (argc < 3) {
    opserr << OpenSees::PromptValueError << "Missing required arguments\n";
    return TCL_ERROR;
  }


  // get the tag of the node
  int nodeId;
  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid tag\n";
    return TCL_ERROR;
  }

  // Alternate form:
  // 
  // fix $node -dof $dof <-value $value>
  //
  if (strcmp(argv[2], "-dof") == 0) {
    if (argc < 4) {
      opserr << OpenSees::PromptValueError << "missing required argument for -dof $dof\n";
      return TCL_ERROR;
    }
    int dof;
    if (Tcl_GetInt(interp, argv[3], &dof) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid dof\n";
      return TCL_ERROR;
    }
    // create a homogeneous constraint
    SP_Constraint *theSP = new SP_Constraint(nodeId, dof-1, 0.0, true, true);

    // add it to the domain
    if (theTclDomain->addSP_Constraint(theSP) == false) {
      opserr << OpenSees::PromptValueError << "could not add SP_Constraint to domain using fix "
                "command - node may already be constrained\n";
      delete theSP;
      return TCL_ERROR;
    }

    return TCL_OK;
  }

  int ndf = argc - 2;

  // check number of arguments
  if (argc < (2 + ndf)) {
    opserr << OpenSees::PromptValueError 
           << "bad command - want: fix nodeId " << ndf
           << " [0,1] conditions";
    return TCL_ERROR;
  }

  // get the fixity condition and add the constraint if fixed
  Tcl_Obj* list = Tcl_NewListObj(ndf, nullptr);
  for (int i = 0; i < ndf; ++i) {
    int theFixity;
    if (Tcl_GetInt(interp, argv[2 + i], &theFixity) != TCL_OK) {
      opserr << OpenSees::PromptValueError
             << "invalid fixity " << i + 1
             << " - load " << nodeId;
      opserr << " " << ndf << " fixities\n";
      return TCL_ERROR;

    } else {
      if (i+1 > builder->getNDF()) {
        opserr << OpenSees::PromptValueError
               << "dof " << i + 1 << " not allowed with NDF = "
               << builder->getNDF() << "\n";
        continue;
      }
      if (theFixity != 0) {
        // create a homogeneous constraint
        SP_Constraint *theSP = new SP_Constraint(nodeId, i, 0.0, true, true);

        // add it to the domain
        if (theTclDomain->addSP_Constraint(theSP) == false) {
          opserr << OpenSees::PromptValueError << "could not add SP_Constraint to domain using fix "
                    "command - node may already be constrained\n";
          delete theSP;
          return TCL_ERROR;

        } else {
//        Tcl(buffer, "%d ", theSP->getTag());
//        Tcl_AppendResult(interp, buffer, NULL);
          Tcl_ListObjAppendElement(interp, list, Tcl_NewDoubleObj(theSP->getTag()));
        }
      }
    }
  }

  Tcl_SetObjResult(interp, list);

  return TCL_OK;
}

int
TclCommand_addHomogeneousBC_X(ClientData clientData, Tcl_Interp *interp,
                                   Tcl_Size argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  ModelRegistry *builder = static_cast<ModelRegistry*>(clientData);

  int ndf = argc - 2;
  if (strcmp(argv[argc-2],"-tol") == 0)
    ndf -= 2;

  // check number of arguments
  if (argc < (2 + ndf)) {
    opserr << OpenSees::PromptValueError << "bad command - want: fixX xLoc " << ndf 
           << " [0,1] conditions" << "\n";
    return TCL_ERROR;
  }

  // get the xCrd of nodes to be constrained
  double xLoc;
  if (Tcl_GetDouble(interp, argv[1], &xLoc) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid xCrd - fixX xLoc " << ndf 
           << " [0,1] conditions" << "\n";
    return TCL_ERROR;
  }

  // read in the fixities
  ID fixity(ndf);
  for (int i=0; i<ndf; ++i) {
    if (Tcl_GetInt(interp, argv[2+i], &fixity(i)) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid fixity " << i+1 << "\n";
      return TCL_ERROR;
    }
  }

  // set the tolerance, the allowable difference in nodal coordinate and
  // what the value user specified to see if node is constrained or not
  double tol = 1.0e-10;
  if (argc >= (4 + ndf)) {
    if (strcmp(argv[2+ndf],"-tol") == 0)
    if (Tcl_GetDouble(interp, argv[3+ndf], &tol) != TCL_OK) {
      opserr << OpenSees::PromptValueError 
             << "invalid tol specified - fixX " << xLoc 
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }
  }


  // TODO: Why not add SP to Domain directly?
  builder->addSP_Constraint(0, xLoc, fixity, tol);

  return TCL_OK;
}

int
TclCommand_addHomogeneousBC_Y(ClientData clientData, Tcl_Interp *interp,
                                   Tcl_Size argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  ModelRegistry *builder = static_cast<ModelRegistry*>(clientData);

  int ndf = argc - 2;
  if (strcmp(argv[argc-2],"-tol") == 0)
    ndf -= 2;

  // check number of arguments
  if (argc < (2 + ndf)) {
    opserr << OpenSees::PromptValueError << "bad command - want: fixY yLoc " << ndf << " [0,1] conditions";
    return TCL_ERROR;
  }

  // get the yCrd of nodes to be constrained
  double yLoc;
  if (Tcl_GetDouble(interp, argv[1], &yLoc) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid yCrd - fixY yLoc " << ndf << " [0,1] conditions\n";
    return TCL_ERROR;
  }

  // read in the fixities
  ID fixity(ndf);
  for (int i=0; i<ndf; ++i) {
    if (Tcl_GetInt(interp, argv[2+i], &fixity(i)) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid fixity " << i+1 << " - fixY " << yLoc;
      opserr << " " << ndf << " fixities\n";
      return TCL_ERROR;
    }
  }

  // set the tolerance, the allowable difference in nodal coordinate and
  // what the value user specified to see if node is constrained or not
  double tol = 1.0e-10;
  if (argc >= (4 + ndf)) {
    if (strcmp(argv[2+ndf],"-tol") == 0)
      if (Tcl_GetDouble(interp, argv[3+ndf], &tol) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "invalid tol specified - fixY " 
               << yLoc << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
  }

  builder->addSP_Constraint(1, yLoc, fixity, tol);

  return TCL_OK;
}

int
TclCommand_addHomogeneousBC_Z(ClientData clientData, Tcl_Interp *interp,
                              Tcl_Size argc, TCL_Char ** const argv)
{

  assert(clientData != nullptr);
  ModelRegistry *builder = static_cast<ModelRegistry*>(clientData);

  int ndf = argc - 2;
  if (strcmp(argv[argc-2],"-tol") == 0)
    ndf -= 2;

  // check number of arguments
  if (argc < (2 + ndf)) {
    opserr << OpenSees::PromptValueError << "bad command - want: fixZ zLoc " << ndf << " [0,1] conditions";
    return TCL_ERROR;
  }

  // get the yCrd of nodes to be constrained
  double zLoc;
  if (Tcl_GetDouble(interp, argv[1], &zLoc) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid zCrd - fixZ zLoc " << ndf << " [0,1] conditions\n";
    return TCL_ERROR;
  }

  // read in the fixities
  ID fixity(ndf);
  for (int i=0; i<ndf; ++i) {
    if (Tcl_GetInt(interp, argv[2+i], &fixity(i)) != TCL_OK) {
      opserr << OpenSees::PromptValueError
             << "invalid fixity " << i+1 << " - fixZ " << zLoc;
      opserr << " " << ndf << " fixities\n";
      return TCL_ERROR;
    }
  }

  // set the tolerance, the allowable difference in nodal coordinate and
  // what the value user specified to see if node is constrained or not
  double tol = 1.0e-10;
  if (argc >= (4 + ndf)) {
    if (strcmp(argv[2+ndf],"-tol") == 0)
      if (Tcl_GetDouble(interp, argv[3+ndf], &tol) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "invalid tol specified - fixZ " << zLoc 
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
  }

  builder->addSP_Constraint(2, zLoc, fixity, tol);

  return TCL_OK;
}



int
TclCommand_addSP(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc,
                      TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  ModelRegistry* builder = static_cast<ModelRegistry*>(clientData);
  Domain *theTclDomain = builder->getDomain();

  if (argc > 1 && (strcmp(argv[1], "remove") == 0)) {
    if (argc < 3) {
      opserr << OpenSees::PromptValueError 
             << "want - remove sp spTag? -or- remove "
                "sp nodeTag? dofTag? <patternTag?>\n";
      return TCL_ERROR;
    }
    int tag;
    if (argc == 3) {
      if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "remove sp tag? failed to read tag: " << argv[2]
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      SP_Constraint *theSPconstraint = theTclDomain->removeSP_Constraint(tag);
      if (theSPconstraint != nullptr) {
        delete theSPconstraint;
      }
    } else {
      int nodeTag, dofTag;
      int patternTag = -1;

      if (Tcl_GetInt(interp, argv[2], &nodeTag) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "remove sp tag? failed to read node tag: " << argv[2]
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      if (Tcl_GetInt(interp, argv[3], &dofTag) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "remove sp tag? failed to read dof tag: " << argv[3]
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }

      if (argc == 5) {
        if (Tcl_GetInt(interp, argv[4], &patternTag) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "remove sp tag? failed to read pattern tag: "
                 << argv[4] 
                 << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        }
      }
      dofTag--; // one for C++ indexing of dof

      theTclDomain->removeSP_Constraint(nodeTag, dofTag, patternTag);

    }
    return TCL_OK;
  }

  LoadPattern *theTclLoadPattern = nullptr;
  if ((theTclLoadPattern = builder->getCurrentPattern<MultiSupportPattern>()))
    ;
  else if ((theTclLoadPattern = builder->getCurrentPattern<StaticPattern>()))
    ;
  else {
    opserr << OpenSees::PromptValueError
            << "no current load pattern supports single point constraints\n";
    return TCL_ERROR;
  }

  // check number of arguments
  if (argc < 4) {
    opserr << OpenSees::PromptValueError
           << "bad command - want: sp nodeId dofID value";
    return TCL_ERROR;
  }

  // get the nodeID, dofId and value of the constraint

  int nodeId, dofId;
  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid nodeId: " << argv[1] << " -  sp nodeId dofID value\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &dofId) != TCL_OK || dofId < 1) {
    opserr << OpenSees::PromptValueError << "invalid dofId: " << argv[2] << " -  sp ";
    opserr << nodeId << " dofID value\n";
    return TCL_ERROR;
  }

  // Decrement the DOF index by 1 to go to C/C++ 0-indexing
  dofId--; 

  double value;
  if (Tcl_GetDouble(interp, argv[3], &value) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid value: " << argv[3] << " -  sp ";
    opserr << nodeId << " dofID value\n";
    return TCL_ERROR;
  }

  bool isSpConst = false;
  bool retZeroInit = true;
  bool userSpecifiedPattern = false;
  int loadPatternTag = 0; // some pattern that will never be used!

  int endMarker = 4;
  while (endMarker != argc) {
    if (strcmp(argv[endMarker],"-const") == 0) {
      // allow user to specify const load
      isSpConst = true;
    }
    else if (strcmp(argv[endMarker], "-subtractInit") == 0) {
      retZeroInit = false;
    }
    else if (strcmp(argv[endMarker],"-pattern") == 0) {
      // allow user to specify load pattern other than current
      endMarker++;
      userSpecifiedPattern = true;
      if (endMarker == argc ||
          Tcl_GetInt(interp, argv[endMarker], &loadPatternTag) != TCL_OK) {

        opserr << OpenSees::PromptValueError << "invalid patternTag - load " << nodeId << "\n";
        return TCL_ERROR;
      }
    }
    endMarker++;
  }

  // if load pattern tag has not changed - get the pattern tag from current one
  if (userSpecifiedPattern == false) {
    if (theTclLoadPattern == nullptr) {
      opserr << OpenSees::PromptValueError << "no current pattern - sp " 
             << nodeId << " dofID value\n"; 
      return TCL_ERROR;
    } else {
      loadPatternTag = theTclLoadPattern->getTag();
    }
  }

  // LoadPattern *thePattern = theTclDomain->getLoadPattern(loadPatternTag);

  // create a homogeneous constraint
  SP_Constraint *theSP = new SP_Constraint(nodeId, dofId, value, isSpConst, retZeroInit);

  if (theTclDomain->addSP_Constraint(theSP, loadPatternTag) == false) {
    opserr << OpenSees::PromptValueError << "could not add SP_Constraint to domain ";
    delete theSP;
    return TCL_ERROR;
  }

  return TCL_OK;
}

int
TclCommand_addEqualDOF_MP(ClientData clientData, Tcl_Interp *interp,
                                Tcl_Size argc, TCL_Char ** const argv)
{
    ModelRegistry *builder = static_cast<ModelRegistry*>(clientData);
    Domain     *theTclDomain   = builder->getDomain();


    // Check number of arguments
    if (argc < 4) {
      opserr << OpenSees::PromptValueError 
             << "bad command - want: equalDOF RnodeID? CnodeID? DOF1? DOF2? ...";
      return TCL_ERROR;
    }

    // Read in the node IDs and the DOF
    int RnodeID, CnodeID, dofID;

    if (Tcl_GetInt(interp, argv[1], &RnodeID) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid RnodeID: " << argv[1]
             << " equalDOF RnodeID? CnodeID? DOF1? DOF2? ...";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[2], &CnodeID) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid CnodeID: " << argv[2]
             << " equalDOF RnodeID? CnodeID? DOF1? DOF2? ...";
      return TCL_ERROR;
    }

    // The number of DOF to be coupled
    int numDOF = argc - 3;

    // The constraint matrix ... U_c = C_cr * U_r
    Matrix Ccr (numDOF, numDOF);
    Ccr.Zero();

    // The vector containing the retained and constrained DOFs
    ID rcDOF (numDOF);

    int i, j;
    // Read the degrees of freedom which are to be coupled
    for (i = 3, j = 0; i < argc; i++, j++) {
      if (Tcl_GetInt(interp, argv[i], &dofID) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid dofID: " << argv[3]
               << " equalDOF RnodeID? CnodeID? DOF1? DOF2? ...";
        return TCL_ERROR;
      }

      dofID -= 1; // Decrement for C++ indexing
      if (dofID < 0) {
        opserr << OpenSees::PromptValueError << "invalid dofID: " << argv[i]
               << " must be >= 1";
        return TCL_ERROR;
      }
      rcDOF (j) = dofID;
      Ccr (j,j) = 1.0;
    }

    // Create the multi-point constraint
    MP_Constraint *theMP = new MP_Constraint(RnodeID, CnodeID, Ccr, rcDOF, rcDOF);

    // Add the multi-point constraint to the domain
    if (theTclDomain->addMP_Constraint(theMP) == false) {
      opserr << OpenSees::PromptValueError << "could not add equalDOF MP_Constraint to domain ";
      delete theMP;
      return TCL_ERROR;
    }
    Tcl_SetObjResult(interp, Tcl_NewIntObj(theMP->getTag()));
    return TCL_OK;
}


int
TclCommand_constrain(ClientData clientData, Tcl_Interp *interp,
                     Tcl_Size argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  ModelRegistry *builder = static_cast<ModelRegistry*>(clientData);
  Domain *domain = builder->getDomain();

  // Usage:
  //   constrain Rnode Cnode {d1 d2 ...} <-rotate {v1 v2 v3}>
  //
  // where v is a rotation vector (axis * angle in radians), R = Exp(Hat(v)),
  // and the enforced relation is Uc = R * Ur on the selected DOFs.

  Tcl_Size argi = 1;

  if (argc < 3) {
    opserr << OpenSees::PromptValueError
           << "insufficient arguments.\n"
           << "Usage: constrain Translation Rnode Cnode {d1 d2 ...} <-rotate {v1 v2 v3}>"
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  //
  // Read tags for the retained and constrained nodes
  //
  int RnodeID = 0, CnodeID = 0;
  if (Tcl_GetInt(interp, argv[argi + 0], &RnodeID) != TCL_OK) {
    opserr << OpenSees::PromptValueError
           << "invalid RnodeID: " << argv[argi + 0]
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[argi + 1], &CnodeID) != TCL_OK) {
    opserr << OpenSees::PromptValueError
           << "invalid CnodeID: " << argv[argi + 1]
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  Node *theRetainedNode    = domain->getNode(RnodeID);
  Node *theConstrainedNode = domain->getNode(CnodeID);
  if (theRetainedNode == nullptr || theConstrainedNode == nullptr) {
    opserr << OpenSees::PromptValueError
           << "Retained or Constrained node does not exist in the Domain: "
           << RnodeID << " " << CnodeID
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  // Parse DOF list argument. This is a Tcl list in one argv slot:
  // {d1 d2 ...}
  Tcl_Size listSize = 0;
  const char **listDOF = nullptr;
  if (Tcl_SplitList(interp, argv[argi + 2], &listSize, &listDOF) != TCL_OK) {
    opserr << OpenSees::PromptValueError
           << "invalid dof list: " << argv[argi + 2]
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }
  if (listSize <= 0) {
    Tcl_Free((char*)listDOF);
    opserr << OpenSees::PromptValueError
           << "empty dof list"
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  const int numDOF = static_cast<int>(listSize);

  // Prepare IDs and store 0-based indices of requested DOFs
  ID rDOF(numDOF);
  ID cDOF(numDOF);
  std::vector<int> dofIdx(numDOF, -1);

 
  
  // Check that nodes have valid spatial dimension
  int dim = theRetainedNode->getCrds().Size();
  if (dim < 1) {
    Tcl_Free((char*)listDOF);
    opserr << OpenSees::PromptValueError
          << "node coordinate dimension is invalid"
          << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }
  // Ensure both nodes have the same dimension
  if (dim != theConstrainedNode->getCrds().Size()) {
    Tcl_Free((char*)listDOF);
    opserr << OpenSees::PromptValueError
          << "retained and constrained node coordinate dimensions do not match"
          << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  // Parse and validate DOF indices
  for (int i = 0; i < numDOF; ++i) {
    int dof1 = 0;
    if (Tcl_GetInt(interp, listDOF[i], &dof1) != TCL_OK) {
      Tcl_Free((char*)listDOF);
      opserr << OpenSees::PromptValueError
             << "invalid dof value: " << listDOF[i]
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }
    if (dof1 <= 0) {
      Tcl_Free((char*)listDOF);
      opserr << OpenSees::PromptValueError
             << "dof values must be >= 1, got: " << dof1
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }

    // Convert to 0-based index
    const int d = dof1 - 1; // 0-based

    // Translation DOFs should live in [0, dim-1]
    if (d < 0 || d >= dim) {
      Tcl_Free((char*)listDOF);
      opserr << OpenSees::PromptValueError
             << "Translation constraint DOF out of range for node dimension. "
             << "Got dof=" << dof1 
             << " but node dimension is " << dim
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }

    // Check duplicates
    for (int k = 0; k < i; ++k) {
      if (dofIdx[k] == d) {
        Tcl_Free((char*)listDOF);
        opserr << OpenSees::PromptValueError
               << "duplicate dof in list: " << dof1
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
    }

    dofIdx[i] = d;
    cDOF(i)   = d;
    rDOF(i)   = d;
  }

  Tcl_Free((char*)listDOF);

  //
  // Parse optional -rotate {v1 v2 v3}
  //
  bool doRotate = false;
  Vector3D v{}; // rotation vector

  for (Tcl_Size a = argi + 3; a < argc; ++a) {
    if (strcmp(argv[a], "-rotate") == 0) {
      if (a + 1 >= argc) {
        opserr << OpenSees::PromptValueError
               << "missing vector after -rotate"
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      if (strcmp(argv[a + 1], "None") == 0) {
        argi++;
        continue; // no rotation
      }

      Tcl_Size nvec = 0;
      const char **vec = nullptr;
      if (Tcl_SplitList(interp, argv[a + 1], &nvec, &vec) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "invalid rotate vector list: " << argv[a + 1]
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }

      if (nvec != 3) {
        Tcl_Free((char*)vec);
        opserr << OpenSees::PromptValueError
               << "-rotate expects 3 components {v1 v2 v3}"
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, vec[0], &v[0]) != TCL_OK ||
          Tcl_GetDouble(interp, vec[1], &v[1]) != TCL_OK ||
          Tcl_GetDouble(interp, vec[2], &v[2]) != TCL_OK) {
        Tcl_Free((char*)vec);
        opserr << OpenSees::PromptValueError
               << "invalid rotate vector components: " << argv[a + 1]
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }

      Tcl_Free((char*)vec);
      doRotate = true;
      ++a; // skip the vector argument
    }
    else {
      opserr << OpenSees::PromptValueError
             << "unknown option: " << argv[a]
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }
  }

  //
  // Build constraint matrix Ccr such that Uc = Ccr * Ur
  //
  Matrix Ccr(numDOF, numDOF);
  Ccr.Zero();

  if (!doRotate) {
    // Identity -> equalDOF behavior on the selected DOFs
    for (int i = 0; i < numDOF; ++i)
      Ccr(i, i) = 1.0;
  }
  else {
    // Compute 3x3 rotation R = Exp(v) using Rodrigues’ formula
    Matrix3D R = ExpSO3(v);

    opserr << "R = " << Matrix(R);

    // Fill submatrix of R corresponding to requested DOFs
    for (int i = 0; i < numDOF; ++i) {
      const int ii = dofIdx[i]; // 0..(dim-1)
      for (int j = 0; j < numDOF; ++j) {
        const int jj = dofIdx[j];
        Ccr(i, j) = R(ii, jj);
      }
    }
  }

  opserr << "Ccr = " << Ccr;

  // Create and add MP (multi-point) constraint
  MP_Constraint *theMP = new MP_Constraint(RnodeID, CnodeID, Ccr, cDOF, rDOF);

  if (domain->addMP_Constraint(theMP) == false) {
    opserr << OpenSees::PromptValueError
           << "could not add MP_Constraint to domain"
           << OpenSees::SignalMessageEnd;
    delete theMP;
    return TCL_ERROR;
  }

  Tcl_SetObjResult(interp, Tcl_NewIntObj(theMP->getTag()));
  return TCL_OK;
}


#if 0
int
TclCommand_addEqualDOF_MP_Mixed(ClientData clientData, Tcl_Interp *interp,
                                Tcl_Size argc, TCL_Char ** const argv)
{
  // Ensure the destructor has not been called
  ModelRegistry *builder = static_cast<ModelRegistry*>(clientData);

  if (theTclBuilder == 0 || clientData == 0) {
    opserr << OpenSees::PromptValueError << "builder has been destroyed - equalDOF \n";
    return TCL_ERROR;
  }

  // Check number of arguments
  if (argc < 4) {
    opserr << OpenSees::PromptValueError << "bad command - want: equalDOFmixed RnodeID? CnodeID? numDOF? RDOF1? CDOF1? ... ...";
    return TCL_ERROR;
  }

  // Read in the node IDs and the DOF
  int RnodeID, CnodeID, dofIDR, dofIDC, numDOF;

  if (Tcl_GetInt(interp, argv[1], &RnodeID) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid RnodeID: " << argv[1]
          << " equalDOF RnodeID? CnodeID? numDOF? RDOF1? CDOF1? ...";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &CnodeID) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid CnodeID: " << argv[2]
          << " equalDOF RnodeID? CnodeID? numDOF? RDOF1? CDOF1? ...";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[3], &numDOF) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid numDOF: " << argv[2]
          << " equalDOF RnodeID? CnodeID? numDOF? RDOF1? CDOF1? ...";
    return TCL_ERROR;
  }

  // The number of DOF to be coupled
  //        int numDOF = argc - 3;

  // The constraint matrix ... U_c = C_cr * U_r
  Matrix Ccr (numDOF, numDOF);
  Ccr.Zero();

  // The vector containing the retained and constrained DOFs
  ID rDOF (numDOF);
  ID cDOF (numDOF);

  int i, j, k;
  // Read the degrees of freedom which are to be coupled
  for (i = 4, j = 5, k = 0; k < numDOF; i+=2, j+=2, k++) {
    if (Tcl_GetInt(interp, argv[i], &dofIDR) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid dofID: " << argv[3]
            << " equalDOF RnodeID? CnodeID? DOF1? DOF2? ...";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[j], &dofIDC) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid dofID: " << argv[3]
            << " equalDOF RnodeID? CnodeID? DOF1? DOF2? ...";
      return TCL_ERROR;
    }

    dofIDR -= 1; // Decrement for 0-based indexing
    dofIDC -= 1;
    if (dofIDC < 0 || dofIDR < 0) {
      opserr << OpenSees::PromptValueError << "invalid dofID: " << argv[i]
              << " must be >= 1";
      return TCL_ERROR;
    }
    rDOF(k) = dofIDR;
    cDOF(k) = dofIDC;
    Ccr(k,k) = 1.0;
  }

  // Create the multi-point constraint
  MP_Constraint *theMP = new MP_Constraint (RnodeID, CnodeID, Ccr, cDOF, rDOF); 

  // Add the multi-point constraint to the domain
  if (theTclDomain->addMP_Constraint (theMP) == false) {
    opserr << OpenSees::PromptValueError << "could not add equalDOF MP_Constraint to domain ";
    delete theMP;
    return TCL_ERROR;
  }

  Tcl_SetObjResult(interp, Tcl_NewIntObj(theMP->getTag()));
  return TCL_OK;
}
#endif


int
TclCommand_addImposedMotionSP(ClientData clientData, Tcl_Interp *interp,
                              Tcl_Size argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  ModelRegistry *builder = static_cast<ModelRegistry*>(clientData);
  Domain* domain = builder->getDomain();


  // check number of arguments
  if (argc < 4) {
    opserr << OpenSees::PromptValueError 
           << "bad command - want: imposedMotion nodeId dofID gMotionID\n";
    return TCL_ERROR;
  }

  // get the nodeID, dofId and value of the constraint
  int nodeId, dofId, gMotionID;

  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "invalid nodeId: " << argv[1]
           << " - imposedMotion nodeId dofID gMotionID" 
           << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2], &dofId) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid dofId: " << argv[2] 
           << " -  imposedMotion " << nodeId << " dofID gMotionID"
           << "\n";
    return TCL_ERROR;
  }
  dofId--; // DECREMENT THE DOF VALUE BY 1 TO GO TO OUR C++ INDEXING

  if (Tcl_GetInt(interp, argv[3], &gMotionID) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid gMotionID: " << argv[3] << " -  imposedMotion ";
    opserr << nodeId << " dofID gMotionID\n";
    return TCL_ERROR;
  }

  bool alt = false;
  if (argc == 5) {
    if (strcmp(argv[4],"-other") == 0)
      alt = true;
  }

  //
  // check valid node & dof
  //
  Node *theNode = domain->getNode(nodeId);
  if (theNode == nullptr) {
    opserr << OpenSees::PromptValueError 
           << "invalid node " << argv[2] << " node not found\n ";
    return -1;
  }

  int nDof = theNode->getNumberDOF();
  if (dofId < 0 || dofId >= nDof) {
    opserr << OpenSees::PromptValueError 
           << "invalid dofId: " << argv[2] << " dof specified cannot be <= 0 or greater than num dof at nod\n "; 
    return -2;
  }

  MultiSupportPattern *thePattern = builder->getCurrentPattern<MultiSupportPattern>();
  if (thePattern == nullptr) {
    opserr << "ERROR no multi-support pattern found\n";
    return TCL_ERROR;
  }

  int loadPatternTag = thePattern->getTag();

  // create a new ImposedMotionSP
  SP_Constraint *theSP;
  if (alt == true) {
    theSP = new ImposedMotionSP1(nodeId, dofId, loadPatternTag, gMotionID);
  } else {
    theSP = new ImposedMotionSP(nodeId, dofId, loadPatternTag, gMotionID);
  }

  if (thePattern->addSP_Constraint(theSP) == false) {
    opserr << OpenSees::PromptValueError 
           << "could not add SP_Constraint to pattern "
           << OpenSees::SignalMessageEnd;
    delete theSP;
    return TCL_ERROR;
  }

  return TCL_OK;
}

