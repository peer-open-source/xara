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

#include "HyperbolicGapMaterial.h"
#include <BasicModelBuilder.h>
#include <Logging.h>
#include <Parsing.h>
#include <string.h>
#include <tcl.h>

int
TclCommand_HyperbolicGapMaterial(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{
//
// Written: MD
// Created: April 2008
//
  int tag;
  double Kmax, Kur, Rf, Fult, gap;

  if (argc < 8) {
    opserr << "WARNING insufficient number of arguments\n";
    return TCL_ERROR;
  }
  
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial tag\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[3], &Kmax) != TCL_OK) {
    opserr << "WARNING invalid Kmax\n";
    return TCL_ERROR;	
  }

  if (Tcl_GetDouble(interp, argv[4], &Kur) != TCL_OK) {
    opserr << "WARNING invalid Kur\n";
    return TCL_ERROR;	
  }

  if (Tcl_GetDouble(interp, argv[5], &Rf) != TCL_OK) {
    opserr << "WARNING invalid Rf\n";
    return TCL_ERROR;	
  }

  if (Tcl_GetDouble(interp, argv[6], &Fult) != TCL_OK) {
    opserr << "WARNING invalid Fult\n";
    return TCL_ERROR;	
  }

  if (Tcl_GetDouble(interp, argv[7], &gap) != TCL_OK) {
    opserr << "WARNING invalid gap\n";
    return TCL_ERROR;	
  }
  
  if (gap>=0) {
    opserr << "Initial gap size must be negative for compression-only material, setting to negative\n";

    gap = -gap;
  }
  if (Fult>0) {
      opserr << "Fult must be negative for compression-only material, setting to negative\n";
    Fult = -Fult;
  }
  if (Kmax == 0.0) {
      opserr << "Kmax is zero, continuing with Kmax = Fult/0.002\n";
      if (Fult != 0.0)
          Kmax = fabs(Fult)/0.002;
      else {
          opserr << "Kmax and Fult are zero\n";
          return TCL_ERROR;
      }
  }

  UniaxialMaterial *theMaterial = new HyperbolicGapMaterial(tag, Kmax, Kur, Rf, Fult, gap);

  BasicModelBuilder* builder = (BasicModelBuilder*)(clientData);
  if (builder->addTaggedObject<UniaxialMaterial>(*theMaterial) == TCL_OK)
    return TCL_OK;
  else
    return TCL_ERROR;

  return TCL_ERROR;
}
