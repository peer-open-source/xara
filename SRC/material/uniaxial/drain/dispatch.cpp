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
// Written: MHS
// Created: Aug 2001
//
// Description: This file contains the implementation of the
// TclBasicBuilder_addDrainMaterial() function.
#include <assert.h>
#include <Logging.h>
#include <Parsing.h>
#include <ModelRegistry.h>
#include <DrainHardeningMaterial.h>
#include <DrainBilinearMaterial.h>
#include <DrainClough1Material.h>
#include <DrainClough2Material.h>
#include <DrainPinch1Material.h>

#include <Vector.h>
#include <string.h>


int
TclBasicBuilder_addDrainMaterial(ClientData clientData, 
                                 Tcl_Interp *interp,
                                 int argc, 
                                 TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  ModelRegistry *model = static_cast<ModelRegistry*>(clientData);

  if (argc < 3) {
    opserr << "WARNING insufficient number of arguments" 
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  int tag;
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial tag" 
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  UniaxialMaterial *theMaterial = nullptr;

  if (strcmp(argv[1], "Hardening2") == 0 ||
      strcmp(argv[1], "Hardening02") == 0) {
    if (argc < 7) {
      opserr << "WARNING invalid number of arguments" << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }

    double E, sigY, Hiso, Hkin;

    if (Tcl_GetDouble(interp, argv[3], &E) != TCL_OK) {
      opserr << "WARNING invalid E" << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[4], &sigY) != TCL_OK) {
      opserr << "WARNING invalid sigY" << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[5], &Hiso) != TCL_OK) {
      opserr << "WARNING invalid Hiso" << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[6], &Hkin) != TCL_OK) {
      opserr << "WARNING invalid Hkin" << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }

    theMaterial = new DrainHardeningMaterial(tag, E, sigY, Hiso, Hkin);
  }

  else if (strcmp(argv[1], "BiLinear") == 0) {
    if (argc < 19) {
      opserr << "WARNING insufficient arguments" << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }

    Vector input(16);
    double temp;

    for (int i = 3, j = 0; j < 16; i++, j++) {
      if (Tcl_GetDouble(interp, argv[i], &temp) != TCL_OK) {
        opserr << "WARNING invalid input, data " << i << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      input(j) = temp;
    }

    theMaterial = new DrainBilinearMaterial(tag, input);
  }

  else if (strcmp(argv[1], "Clough1") == 0) {
    if (argc < 19) {
      opserr << "WARNING insufficient arguments" << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }

    Vector input(16);
    double temp;

    for (int i = 3, j = 0; j < 16; i++, j++) {
      if (Tcl_GetDouble(interp, argv[i], &temp) != TCL_OK) {
        opserr << "WARNING invalid input, data " << i << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      input(j) = temp;
    }

    theMaterial = new DrainClough1Material(tag, input);
  }

  else if (strcmp(argv[1], "Clough2") == 0) {
    if (argc < 19) {
      opserr << "WARNING insufficient arguments" << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }

    Vector input(16);
    double temp;

    for (int i = 3, j = 0; j < 16; i++, j++) {
      if (Tcl_GetDouble(interp, argv[i], &temp) != TCL_OK) {
        opserr << "WARNING invalid input, data " << i << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      input(j) = temp;
    }

    theMaterial = new DrainClough2Material(tag, input);
  }

  else if (strcmp(argv[1], "Pinch1") == 0) {
    if (argc < 22) {
      opserr << "WARNING insufficient arguments" << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }

    Vector input(19);
    double temp;

    for (int i = 3, j = 0; j < 19; i++, j++) {
      if (Tcl_GetDouble(interp, argv[i], &temp) != TCL_OK) {
        opserr << "WARNING invalid input, data " << i << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      input(j) = temp;
    }

    theMaterial = new DrainPinch1Material(tag, input);
  }

  if (theMaterial == nullptr) {
    opserr << "WARNING unknown material type: " << argv[1] << OpenSees::SignalMessageEnd;
    return TCL_ERROR;
  }

  if (model->addTaggedObject<UniaxialMaterial>(*theMaterial) != TCL_OK) {
    delete theMaterial;
    return TCL_ERROR;
  }
  return TCL_OK;
}
