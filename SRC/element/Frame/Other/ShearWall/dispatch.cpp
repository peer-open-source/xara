//
// Description: This file contains the 
// implementation of TclBasicBuilder_addDispBeamColumnInt().
// Based on TclBasicBuilder_addDispBeamColumn().
//
// Created: 07/04 
// Modified by: LMS 
//
#include <tcl.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <Domain.h>
#include <FrameTransform.h>

#include "DispBeamColumn2dInt.h"
#include <runtimeAPI.h>
#include <BasicModelBuilder.h>
#include <SectionForceDeformation.h>


int
TclBasicBuilder_addDispBeamColumnInt(ClientData clientData, Tcl_Interp *interp,
                                     int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;

  int ndm = builder->getNDM();
  int ndf = builder->getNDF();

  int ok = 0;
  if (ndm == 2 && ndf == 3)
    ok = 1;

  if (ok == 0) {
    opserr << "WARNING -- NDM = " << ndm << " and NDF = " << ndf
           << " not compatible with dispBeamColumn element" << endln;
    return TCL_ERROR;
  }

  if (argc < 9) { // 8
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element dispBeamColumn eleTag? iNode? jNode? nIP? secTag? "
              "transfTag? C1? t1? NStrip1? t2? NStrip2? t3? NStrip3?\n";
    return TCL_ERROR;
  }

  // get the id and end nodes
  int eleTag, iNode, jNode, nIP, transfTag;
  double C1;
  int secTag[10]; // Max size of integration rule ... can change if needed
  int argi = 2;

  if (Tcl_GetInt(interp, argv[argi++], &eleTag) != TCL_OK) {
    opserr << "WARNING invalid dispBeamColumn eleTag" << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[argi++], &iNode) != TCL_OK) {
    opserr << "WARNING invalid iNode ";
    opserr << "dispBeamColumn element: " << eleTag << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[argi++], &jNode) != TCL_OK) {
    opserr << "WARNING invalid jNode ";
    opserr << "dispBeamColumn element: " << eleTag << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[argi++], &nIP) != TCL_OK) {
    opserr << "WARNING invalid nIP ";
    opserr << "dispBeamColumn element: " << eleTag << endln;
    return TCL_ERROR;
  }

  if (strcmp(argv[argi], "-sections") == 0) {
    argi++;
    if (argi + nIP > argc) {
      opserr
          << "WARNING insufficient number of section tags - element "
             "dispBeamColumn eleTag? iNode? jNode? nIP? secTag? transfTag?\n";
      return TCL_ERROR;
    }
    int section;
    for (int i = 0; i < nIP; i++) {
      if (Tcl_GetInt(interp, argv[argi + i], &section) != TCL_OK) {
        opserr << "WARNING invalid secTag - element dispBeamColumn eleTag? "
                  "iNode? jNode? nIP? secTag? transfTag?\n";
        return TCL_ERROR;
      }
      secTag[i] = section;
    }
    argi += nIP;
  }

  else {
    int section;
    if (Tcl_GetInt(interp, argv[argi++], &section) != TCL_OK) {
      opserr << "WARNING invalid secTag - element dispBeamColumn eleTag? "
                "iNode? jNode? nIP? secTag? transfTag?\n";
      return TCL_ERROR;
    }
    for (int i = 0; i < nIP; i++)
      secTag[i] = section;
  }

  if (argi >= argc || Tcl_GetInt(interp, argv[argi++], &transfTag) != TCL_OK) {
    opserr << "WARNING invalid transfTag? - element dispBeamColumn eleTag? "
              "iNode? jNode? nIP? secTag? transfTag?\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[argi++], &C1) != TCL_OK) {
    opserr << "WARNING invalid dispBeamColumn C1" << endln;
    return TCL_ERROR;
  }

  double massDens = 0.0;

  while (argi != argc) {
    if (strcmp(argv[argi++], "-mass") == 0 && argi < argc) {
      if (Tcl_GetDouble(interp, argv[argi++], &massDens) != TCL_OK) {
        opserr << "WARNING invalid massDens - element dispBeamColumn eleTag? "
                  "iNode? jNode? nIP? secTag? transfTag? C1? t? NStrip?\n";
        return TCL_ERROR;
      }
    }
  }

  SectionForceDeformation **sections = new SectionForceDeformation *[nIP];

  for (int j = 0; j < nIP; j++) {
    SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(secTag[j]);
    if (theSection == nullptr) {
      delete[] sections;
      return TCL_ERROR;
    }

    sections[j] = theSection;
  }

  Element *theElement = nullptr;

  if (ndm == 2) {

    FrameTransform2d *theTransf = builder->getTypedObject<FrameTransform2d>(transfTag);
    if (theTransf == nullptr) {
      delete[] sections;
      return TCL_ERROR;
    }

    // now create the DispBeamColumn and add it to the Domain
    theElement = new DispBeamColumn2dInt(eleTag, iNode, jNode, nIP, sections,
                                         *theTransf, C1, massDens);

    delete[] sections;
  }

  if (theElement == nullptr) {
    opserr << "WARNING ran out of memory creating element\n";
    opserr << "dispBeamColumn element: " << eleTag << endln;
    return TCL_ERROR;
  }

  if (builder->getDomain()->addElement(theElement) == false) {
    opserr << "WARNING could not add element to the domain\n";
    opserr << "dispBeamColumn element: " << eleTag << endln;
    delete theElement;
    return TCL_ERROR;
  }

  // if get here we have successfully created the element and added it to the
  // domain
  return TCL_OK;
}
