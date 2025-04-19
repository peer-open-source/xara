/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */
//
// Description: This file contains the function to parse the TCL input
// for the genericClient element.
//
// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 11/06
//
#include <tcl.h>
#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <BasicModelBuilder.h>
#include <GenericClient.h>
#include <GenericCopy.h>

int
TclBasicBuilder_addGenericClient(ClientData clientData, Tcl_Interp *interp,
                                 int argc, TCL_Char ** const argv)
{
  const int eleArgStart = 1;

  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;
  Domain* theTclDomain = builder->getDomain();

  if (builder == 0 || clientData == 0) {
    opserr << "WARNING builder has been destroyed - genericClient\n";
    return TCL_ERROR;
  }

  // check the number of arguments is correct
  if ((argc - eleArgStart) < 8) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element genericClient eleTag -node Ndi Ndj ... -dof "
              "dofNdi -dof dofNdj ... -server ipPort <ipAddr> <-ssl> <-udp> "
              "<-dataSize size> <-noRayleigh>\n";
    return TCL_ERROR;
  }

  Element *theElement = nullptr;

  // get the id and end nodes
  int tag, node, dof, ipPort;
  char *ipAddr = 0;
  int ssl = 0, udp = 0;
  int dataSize = 256;
  int doRayleigh = 1;

  if (Tcl_GetInt(interp, argv[1 + eleArgStart], &tag) != TCL_OK) {
    opserr << "WARNING invalid genericClient eleTag\n";
    return TCL_ERROR;
  }
  // read the number of nodes
  if (strcmp(argv[2 + eleArgStart], "-node") != 0) {
    opserr << "WARNING expecting -node flag\n";
    return TCL_ERROR;
  }
  int numNodes = 0;
  int argi = 3 + eleArgStart;
  for (int i=argi; strcmp(argv[i], "-dof") != 0 && i < argc;) {
    numNodes++;
    i++;
  }
  if (numNodes == 0) {
    opserr << "WARNING no nodes specified\n";
    return TCL_ERROR;
  }
  // create the ID arrays to hold the nodes and dofs
  ID nodes(numNodes);
  ID *dofs = new ID[numNodes];

  // fill in the nodes ID
  for (int i = 0; i < numNodes; ++i) {
    if (Tcl_GetInt(interp, argv[argi], &node) != TCL_OK) {
      opserr << "WARNING invalid node\n";
      opserr << "genericClient element: " << tag << endln;
      return TCL_ERROR;
    }
    nodes(i) = node;
    argi++;
  }
  for (int j = 0; j < numNodes; j++) {
    // read the number of dofs per node j
    int numDOFj = 0;
    if (strcmp(argv[argi], "-dof") != 0) {
      opserr << "WARNING expect -dof\n";
      opserr << "genericClient element: " << tag << endln;
      return TCL_ERROR;
    }
    argi++;

    for (int i=argi; 
         strcmp(argv[i], "-dof") != 0 && strcmp(argv[i], "-server") != 0 &&
         strcmp(argv[i], "-doRayleigh") != 0 &&
         strcmp(argv[i], "-noRayleigh") != 0 && i < argc; 
         i++) {
      numDOFj++;
    }
    // fill in the dofs ID array
    ID dofsj(numDOFj);
    for (int i = 0; i < numDOFj; ++i) {
      if (Tcl_GetInt(interp, argv[argi], &dof) != TCL_OK) {
        opserr << "WARNING invalid dof\n";
        opserr << "genericClient element: " << tag << endln;
        return TCL_ERROR;
      }
      dofsj(i) = dof - 1;
      argi++;
    }
    dofs[j] = dofsj;
  }
  if (strcmp(argv[argi], "-server") == 0) {
    argi++;
    if (Tcl_GetInt(interp, argv[argi], &ipPort) != TCL_OK) {
      opserr << "WARNING invalid ipPort\n";
      opserr << "genericClient element: " << tag << endln;
      return TCL_ERROR;
    }
    argi++;
    if (argi < argc && strcmp(argv[argi], "-doRayleigh") != 0 &&
        strcmp(argv[argi], "-noRayleigh") != 0 &&
        strcmp(argv[argi], "-dataSize") != 0 &&
        strcmp(argv[argi], "-ssl") != 0 && strcmp(argv[argi], "-udp") != 0) {
      ipAddr = new char[strlen(argv[argi]) + 1];
      strcpy(ipAddr, argv[argi]);
      argi++;
    } else {
      ipAddr = new char[9 + 1];
      strcpy(ipAddr, "127.0.0.1");
    }
    for (int i = argi; i < argc; ++i) {
      if (strcmp(argv[i], "-ssl") == 0) {
        ssl = 1;
        udp = 0;
      } else if (strcmp(argv[i], "-udp") == 0) {
        udp = 1;
        ssl = 0;
      } else if (strcmp(argv[i], "-dataSize") == 0) {
        if (Tcl_GetInt(interp, argv[i + 1], &dataSize) != TCL_OK) {
          opserr << "WARNING invalid dataSize\n";
          opserr << "genericClient element: " << tag << endln;
          return TCL_ERROR;
        }
      }
    }
  } else {
    opserr << "WARNING expecting -server string but got ";
    opserr << argv[argi] << endln;
    opserr << "genericClient element: " << tag << endln;
    return TCL_ERROR;
  }
  for (int i = argi; i < argc; ++i) {
    if (strcmp(argv[i], "-doRayleigh") == 0) {
      doRayleigh = 1;
    } else if (strcmp(argv[i], "-noRayleigh") == 0) {
      doRayleigh = 0;
    }
  }

  // now create the GenericClient
  theElement = new GenericClient(tag, nodes, dofs, ipPort, ipAddr, ssl, udp,
                                 dataSize, doRayleigh);

  // cleanup dynamic memory
  if (dofs != 0)
    delete[] dofs;


  // then add the GenericClient to the domain
  if (theTclDomain->addElement(theElement) == false) {
    opserr << "WARNING could not add element to the domain\n";
    opserr << "genericClient element: " << tag << endln;
    delete theElement;
    return TCL_ERROR;
  }

  return TCL_OK;
}

int
TclBasicBuilder_addGenericCopy(ClientData clientData, Tcl_Interp *interp, int argc,
                               TCL_Char ** const argv)
{
  const int eleArgStart = 1;
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;
  Domain* theTclDomain = builder->getDomain();

  // check the number of arguments is correct
  if ((argc - eleArgStart) < 6) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: expElement genericCopy eleTag -node Ndi ... -src srcTag\n";
    return TCL_ERROR;
  }


  // get the id and end nodes
  int tag, node, srcTag;
  int numNodes = 0;

  if (Tcl_GetInt(interp, argv[1 + eleArgStart], &tag) != TCL_OK) {
    opserr << "WARNING invalid genericCopy eleTag\n";
    return TCL_ERROR;
  }
  // read the number of nodes
  if (strcmp(argv[2 + eleArgStart], "-node") != 0) {
    opserr << "WARNING expecting -node flag\n";
    return TCL_ERROR;
  }
  int argi = 3 + eleArgStart;

  for (int i=argi; strcmp(argv[i], "-src") != 0 && i < argc;) {
    numNodes++;
    i++;
  }
  if (numNodes == 0) {
    opserr << "WARNING no nodes specified\n";
    return TCL_ERROR;
  }
  // create and fill in the ID array to hold the nodes
  ID nodes(numNodes);
  for (int i = 0; i < numNodes; ++i) {
    if (Tcl_GetInt(interp, argv[argi], &node) != TCL_OK) {
      opserr << "WARNING invalid node\n";
      return TCL_ERROR;
    }
    nodes(i) = node;
    argi++;
  }
  if (strcmp(argv[argi], "-src") != 0) {
    opserr << "WARNING expect -src\n";
    return TCL_ERROR;
  }
  argi++;
  if (Tcl_GetInt(interp, argv[argi], &srcTag) != TCL_OK) {
    opserr << "WARNING invalid srcTag\n";
    return TCL_ERROR;
  }

  // now create the GenericCopy
  Element* theElement = new GenericCopy(tag, nodes, srcTag);

  // then add the GenericCopy to the domain
  if (theTclDomain->addElement(theElement) == false) {
    opserr << "WARNING could not add element to the domain\n";
    opserr << "genericCopy element: " << tag << endln;
    delete theElement;
    return TCL_ERROR;
  }

  return TCL_OK;
}
