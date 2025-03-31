/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
//
// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 09/07
// Revision: A
//
// Description: This file contains the function to parse the TCL input
// for the actuator element.
//
#include <BasicModelBuilder.h>
#include <tcl.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <Domain.h>
#include <Logging.h>
#include <Actuator.h>
#include <Adapter.h>
#include <ActuatorCorot.h>


int
TclCommand_addActuator(ClientData clientData, Tcl_Interp *interp, int argc,
                            TCL_Char ** const argv)
{
  constexpr static int eleArgStart = 1;
  // ensure the destructor has not been called
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;

  // check the number of arguments is correct
  if ((argc - eleArgStart) < 6) {
    opserr << OpenSees::PromptValueError << "insufficient arguments\n";
    opserr << "Want: element actuator eleTag iNode jNode EA ipPort "
              "<-doRayleigh> <-rho rho>\n";
    return TCL_ERROR;
  }

  Element *theElement = 0;
  int ndm = builder->getNDM();

  // get the id and end nodes
  int tag, iNode, jNode;
  double EA;
  int ipPort;
  int doRayleigh = 0;
  double rho = 0.0;

  if (Tcl_GetInt(interp, argv[1 + eleArgStart], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid element tag" << "\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2 + eleArgStart], &iNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid iNode\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3 + eleArgStart], &jNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid jNode\n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[4 + eleArgStart], &EA) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid EA\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[5 + eleArgStart], &ipPort) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid ipPort\n";
    return TCL_ERROR;
  }
  for (int i = 6 + eleArgStart; i < argc; ++i) {
    if (strcmp(argv[i], "-doRayleigh") == 0)
      doRayleigh = 1;
  }
  for (int i = 6 + eleArgStart; i < argc; ++i) {
    if (i + 1 < argc && strcmp(argv[i], "-rho") == 0) {
      if (Tcl_GetDouble(interp, argv[i + 1], &rho) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid rho\n";
        return TCL_ERROR;
      }
    }
  }

  // now create the actuator and add it to the Domain
  theElement =
      new Actuator(tag, ndm, iNode, jNode, EA, ipPort, doRayleigh, rho);


  Domain* theTclDomain = builder->getDomain();
  if (theTclDomain->addElement(theElement) == false) {
    opserr << OpenSees::PromptValueError << "could not add element to the domain\n";
    delete theElement;
    return TCL_ERROR;
  }

  // if get here we have successfully created the actuator and added it to the
  // domain
  return TCL_OK;
}

int
TclCommand_addActuatorCorot(ClientData clientData, Tcl_Interp *interp,
                                 int argc, TCL_Char ** const argv)
{
  constexpr static int eleArgStart = 1;
  assert(clientData != nullptr);
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;


  // check the number of arguments is correct
  if ((argc - eleArgStart) < 6) {
    opserr << OpenSees::PromptValueError << "insufficient arguments\n";
    opserr << "Want: element corotActuator eleTag iNode jNode EA ipPort "
              "<-doRayleigh> <-rho rho>\n";
    return TCL_ERROR;
  }

  Element *theElement = nullptr;
  int ndm = builder->getNDM();

  // get the id and end nodes
  int tag, iNode, jNode;
  double EA;
  int ipPort;
  int doRayleigh = 0;
  double rho = 0.0;

  if (Tcl_GetInt(interp, argv[1 + eleArgStart], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid corotActuator eleTag" << "\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2 + eleArgStart], &iNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid iNode\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3 + eleArgStart], &jNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid jNode\n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[4 + eleArgStart], &EA) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid EA\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[5 + eleArgStart], &ipPort) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid ipPort\n";
    return TCL_ERROR;
  }
  for (int i = 6 + eleArgStart; i < argc; ++i) {
    if (strcmp(argv[i], "-doRayleigh") == 0)
      doRayleigh = 1;
  }
  for (int i = 6 + eleArgStart; i < argc; ++i) {
    if (i + 1 < argc && strcmp(argv[i], "-rho") == 0) {
      if (Tcl_GetDouble(interp, argv[i + 1], &rho) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid rho\n";
        opserr << "corotActuator element: " << tag << "\n";
        return TCL_ERROR;
      }
    }
  }

  // now create the corotActuator and add it to the Domain
  theElement =
      new ActuatorCorot(tag, ndm, iNode, jNode, EA, ipPort, doRayleigh, rho);


  Domain* theTclDomain = builder->getDomain();
  if (theTclDomain->addElement(theElement) == false) {
    opserr << OpenSees::PromptValueError << "could not add element to the domain\n";
    delete theElement;
    return TCL_ERROR;
  }

  // if get here we have successfully created the corotActuator and added it to
  // the domain
  return TCL_OK;
}


int
TclCommand_addAdapter(ClientData clientData, Tcl_Interp *interp, int argc,
                           TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;

  constexpr static int eleArgStart = 1;
  // check the number of arguments is correct
  if ((argc - eleArgStart) < 8) {
    opserr << OpenSees::PromptValueError << "insufficient arguments\n";
    opserr << "Want: element adapter eleTag -node Ndi Ndj ... -dof dofNdi -dof "
              "dofNdj ... -stif Kij ipPort <-doRayleigh> <-mass Mij>\n";
    return TCL_ERROR;
  }

  // get the id and end nodes
  int tag, node, dof, ipPort, argi, i, j, k;
  int numNodes = 0, 
      numDOFj  = 0, 
      numDOF   = 0;
  int doRayleigh = 0;
  Matrix *mass   = nullptr;
  Element *theElement = nullptr;

  if (Tcl_GetInt(interp, argv[1 + eleArgStart], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid adapter eleTag" << "\n";
    return TCL_ERROR;
  }
  // read the number of nodes
  if (strcmp(argv[2 + eleArgStart], "-node") != 0) {
    opserr << OpenSees::PromptValueError << "expecting -node flag\n";
    return TCL_ERROR;
  }
  argi = 3 + eleArgStart;
  i = argi;
  while (strcmp(argv[i], "-dof") != 0 && i < argc) {
    numNodes++;
    i++;
  }
  if (numNodes == 0) {
    opserr << OpenSees::PromptValueError << "no nodes specified\n";
    return TCL_ERROR;
  }
  // create the ID arrays to hold the nodes and dofs
  ID nodes(numNodes);
  ID *dofs = new ID[numNodes];

  // fill in the nodes ID
  for (i = 0; i < numNodes; ++i) {
    if (Tcl_GetInt(interp, argv[argi], &node) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid node\n";
      return TCL_ERROR;
    }
    nodes(i) = node;
    argi++;
  }
  for (j = 0; j < numNodes; j++) {
    // read the number of dofs per node j
    numDOFj = 0;
    if (strcmp(argv[argi], "-dof") != 0) {
      opserr << OpenSees::PromptValueError << "expect -dof\n";
      return TCL_ERROR;
    }
    argi++;
    i = argi;
    while (strcmp(argv[i], "-dof") != 0 && strcmp(argv[i], "-stif") != 0 &&
           i < argc) {
      numDOFj++;
      numDOF++;
      i++;
    }
    // fill in the dofs ID array
    ID dofsj(numDOFj);
    for (i = 0; i < numDOFj; ++i) {
      if (Tcl_GetInt(interp, argv[argi], &dof) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid dof\n";
        return TCL_ERROR;
      }
      dofsj(i) = dof - 1;
      argi++;
    }
    dofs[j] = dofsj;
  }
  // get stiffness matrix
  Matrix kb(numDOF, numDOF);
  if (strcmp(argv[argi], "-stif") != 0) {
    opserr << OpenSees::PromptValueError << "expecting -stif flag\n";
    return TCL_ERROR;
  }
  argi++;
  if (argc - 1 < argi + numDOF * numDOF) {
    opserr << OpenSees::PromptValueError << "incorrect number of stiffness terms\n";
    return TCL_ERROR;
  }
  double stif;
  for (j = 0; j < numDOF; j++) {
    for (k = 0; k < numDOF; k++) {
      if (Tcl_GetDouble(interp, argv[argi], &stif) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid stiffness term\n";
        return TCL_ERROR;
      }
      kb(j, k) = stif;
      argi++;
    }
  }
  // get ip-port
  if (Tcl_GetInt(interp, argv[argi], &ipPort) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid ipPort\n";
    return TCL_ERROR;
  }
  argi++;
  // get optional rayleigh flag
  for (int i = argi; i < argc; ++i) {
    if (strcmp(argv[i], "-doRayleigh") == 0)
      doRayleigh = 1;
  }
  // get optional mass matrix
  for (int i = argi; i < argc; ++i) {
    if (strcmp(argv[i], "-mass") == 0) {
      if (argc - 1 < i + numDOF * numDOF) {
        opserr << OpenSees::PromptValueError << "incorrect number of mass terms\n";
        return TCL_ERROR;
      }
      mass = new Matrix(numDOF, numDOF);
      double m;
      for (j = 0; j < numDOF; j++) {
        for (k = 0; k < numDOF; k++) {
          if (Tcl_GetDouble(interp, argv[i + 1 + numDOF * j + k], &m) != TCL_OK) {
            opserr << OpenSees::PromptValueError << "invalid mass term\n";
            return TCL_ERROR;
          }
          (*mass)(j, k) = m;
        }
      }
    }
  }

  // now create the adapter and add it to the Domain
  if (mass == 0)
    theElement = new Adapter(tag, nodes, dofs, kb, ipPort, 0, 0, doRayleigh);
  else
    theElement =
        new Adapter(tag, nodes, dofs, kb, ipPort, 0, 0, doRayleigh, mass);

  // cleanup dynamic memory
  if (dofs != 0)
    delete[] dofs;


  Domain* theTclDomain = builder->getDomain();
  if (theTclDomain->addElement(theElement) == false) {
    opserr << OpenSees::PromptValueError << "could not add element to the domain\n";
    delete theElement;
    return TCL_ERROR;
  }

  // if get here we have successfully created the adapter and added it to the
  // domain
  return TCL_OK;
}
