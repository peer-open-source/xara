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
// File: ~/element/TclElementCommands.C
//
// Written: fmk
// Created: 07/99
// Revision: A
//
// Description: This file contains the implementation of the TclElementCommands.
// The file contains the routine TclElementCommands which is invoked by the
// TclBasicBuilder.
//
// What: "@(#) TclBasicBuilder.C, revA"
#include <Parsing.h>
#include <Logging.h>
#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <fElmt02.h>
#include <ModelRegistry.h>



extern "C" int elmt04_(double *d, double *ul, double *xl, int *ix, double *tl, 
                       double *s, double *r, int *ndf, int *ndm, int *nst, int *isw, 
                       double *dm, int *nen, int *n, int *nh1, int *nh2, int *nh3, 
                       double *h, double *ctan, int *ior, int *iow);



extern "C" int elmt05_(double *d, double *ul, double *xl, int *ix, double *tl, 
                       double *s, double *r, int *ndf, int *ndm, int *nst, int *isw, 
                       double *dm, int *nen, int *n, int *nh1, int *nh2, int *nh3, 
                       double *h, double *ctan, int *ior, int *iow);


int
TclBasicBuilder_addFeapTruss(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc,
                             TCL_Char ** const argv, int eleArgStart)
{
  // Ensure the destructor has not been called -
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;

  int ndm = builder->getNDM();
  int ndf = builder->getNDF();

  if (ndm != 2 && ndf != 2) {
    opserr
        << "WARNING - fTruss eleTag? iNode? jNode? A? E? needs ndm=2, ndf=2\n";
    return TCL_ERROR;
  }

  // check the number of arguments is correct
  if ((argc - eleArgStart) < 5) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element fTruss eleTag? iNode? jNode? A? E?\n";

    return TCL_ERROR;
  }

  // Get the id and end nodes
  int trussId, iNode, jNode;
  double A, E;
  if (Tcl_GetInt(interp, argv[1 + eleArgStart], &trussId) != TCL_OK) {
    opserr << "WARNING invalid truss eleTag" << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2 + eleArgStart], &iNode) != TCL_OK) {
    opserr << "WARNING invalid iNode" << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3 + eleArgStart], &jNode) != TCL_OK) {
    opserr << "WARNING invalid jNode" << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[4 + eleArgStart], &A) != TCL_OK) {
    opserr << "WARNING invalid A" << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[5 + eleArgStart], &E) != TCL_OK) {
    opserr << "WARNING invalid E" << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  // now create the truss and add it to the Domain
  fElmt02 *theTruss = new fElmt02(trussId, iNode, jNode, A, E);

  if (builder->getDomain()->addElement(theTruss) == false) {
    opserr << "WARNING could not add element to the domain" << OpenSees::SignalMessageEnd;
    delete theTruss;
    return TCL_ERROR;
  }
  return TCL_OK;
}
