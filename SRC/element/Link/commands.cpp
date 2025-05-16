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
// Created: 08/08
// Revision: A
//
// Description: This file contains the function to parse the TCL input
// for the TwoNodeLink element.

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <tcl.h>
#include <Domain.h>
#include <Logging.h>
#include <Parsing.h>
#include <ID.h>
#include <Vector.h>
#include <UniaxialMaterial.h>
#include <SectionForceDeformation.h>
#include <BasicModelBuilder.h>
#include <TwoNodeLink.h>
#include <TwoNodeLinkSection.h>
#include <VectorND.h>

using OpenSees::VectorND;

template <int ndm>
int
SetupLink(const VectorND<ndm> &xi, const VectorND<ndm>& xj,
          const VectorND<ndm>& x, const VectorND<ndm>& y,
          const ID &direction, const ID &dir,
          Matrix &trans, double &L, bool onP0)
{
#if 0
  const Vector &end1Crd = theNodes[0]->getCrds();
  const Vector &end2Crd = theNodes[1]->getCrds();	
  Vector xp = end2Crd - end1Crd;
  L = xp.Norm();

  // setup x and y orientation vectors
  if (L > DBL_EPSILON)  {
      if (x.Size() == 0)  {
          x.resize(3);
          x.Zero();
          x(0) = xp(0);
          if (xp.Size() > 1)
              x(1) = xp(1);
          if (xp.Size() > 2)
              x(2) = xp(2);

      } else if (onP0)  {
          opserr << "WARNING TwoNodeLink::setUp() - " 
              << "element: " << this->getTag() << endln
              << "ignoring nodes and using specified "
              << "local x vector to determine orientation\n";
      }
      if (y.Size() == 0)  {
          y.resize(3);
          y.Zero();
          y(0) = -xp(1);
          if (xp.Size() > 1)
              y(1) = xp(0);
          if (xp.Size() > 2)
              opserr << "WARNING TwoNodeLink::setUp() - " 
                  << "element: " << this->getTag() << endln
                  << "no local y vector specified\n";
      }
  } else  {
      if (x.Size() == 0)  {
          x.resize(3);
          x(0) = 1.0; x(1) = 0.0; x(2) = 0.0;
      }
      if (y.Size() == 0)  {
          y.resize(3);
          y(0) = 0.0; y(1) = 1.0; y(2) = 0.0;
      }
  }

  // check that vectors for orientation are of correct size
  if (x.Size() != 3 || y.Size() != 3)  {
      opserr << "TwoNodeLink::setUp() - "
          << "element: " << this->getTag() << endln
          << "incorrect dimension of orientation vectors\n";
      exit(-1);
  }

  // establish orientation of element for the transformation matrix
  // z = x cross yp
  static Vector z(3);
  z(0) = x(1)*y(2) - x(2)*y(1);
  z(1) = x(2)*y(0) - x(0)*y(2);
  z(2) = x(0)*y(1) - x(1)*y(0);

  // y = z cross x
  y(0) = z(1)*x(2) - z(2)*x(1);
  y(1) = z(2)*x(0) - z(0)*x(2);
  y(2) = z(0)*x(1) - z(1)*x(0);

  // compute length(norm) of vectors
  double xn = x.Norm();
  double yn = y.Norm();
  double zn = z.Norm();

  // check valid x and y vectors, i.e. not parallel and of zero length
  if (xn == 0 || yn == 0 || zn == 0)  {
      opserr << "TwoNodeLink::setUp() - "
          << "element: " << this->getTag() << endln
          << "invalid orientation vectors\n";
      exit(-1);
  }

  // create transformation matrix of direction cosines
  for (int i=0; i<3; i++)  {
      trans(0,i) = x(i)/xn;
      trans(1,i) = y(i)/yn;
      trans(2,i) = z(i)/zn;
  }
#endif
  return 0;
}


int
TclCommand_addTwoNodeLink(ClientData clientData, Tcl_Interp *interp, int argc,
                               TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;
  constexpr static int eleArgStart = 1;


  Element *theElement = nullptr;
  int ndm = builder->getNDM();
  int ndf = builder->getNDF();

  // check the number of arguments is correct
  if ((argc - eleArgStart) < 8) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: twoNodeLink eleTag iNode jNode -mat matTags -dir dirs "
              "<-orient <x1 x2 x3> y1 y2 y3> <-pDelta Mratios> <-shearDist "
              "sDratios> <-doRayleigh> <-mass m>\n";
    return TCL_ERROR;
  }

  // get the id and end nodes
  int tag, iNode, jNode, numMat, matTag, numDir, dirID;
  int argi = eleArgStart + 1;
  Vector Mratio(0), shearDistI(0);
  int doRayleigh = 0;
  double mass = 0.0;

  if (Tcl_GetInt(interp, argv[argi], &tag) != TCL_OK) {
    opserr << "WARNING invalid eleTag\n";
    return TCL_ERROR;
  }
  argi++;
  if (Tcl_GetInt(interp, argv[argi], &iNode) != TCL_OK) {
    opserr << "WARNING invalid iNode\n";
    return TCL_ERROR;
  }
  argi++;
  if (Tcl_GetInt(interp, argv[argi], &jNode) != TCL_OK) {
    opserr << "WARNING invalid jNode\n";
    return TCL_ERROR;
  }
  argi++;
  // read the number of materials
  numMat = 0;
  if (strcmp(argv[argi], "-mat") != 0) {
    opserr << "WARNING expecting -mat flag\n";
    return TCL_ERROR;
  }
  argi++;
  int i = argi;
  while (i < argc && strcmp(argv[i], "-dir") != 0) {
    numMat++;
    i++;
  }
  if (numMat == 0) {
    opserr << "WARNING no directions specified\n";
    return TCL_ERROR;
  }
  // create array of uniaxial materials
  UniaxialMaterial **theMaterials = new UniaxialMaterial *[numMat];
  for (int i = 0; i < numMat; i++) {
    theMaterials[i] = nullptr;
    if (Tcl_GetInt(interp, argv[argi], &matTag) != TCL_OK) {
      opserr << "WARNING invalid matTag\n";
      opserr << "twoNodeLink element: " << tag << endln;
      return TCL_ERROR;
    }
    theMaterials[i] = builder->getTypedObject<UniaxialMaterial>(matTag);
    if (theMaterials[i] == 0) {
      return TCL_ERROR;
    }
    argi++;
  }
  // read the number of directions
  numDir = 0;
  if (strcmp(argv[argi], "-dir") != 0) {
    opserr << "WARNING expecting -dir flag\n";
    return TCL_ERROR;
  }
  argi++;
  i = argi;
  while (i < argc && strcmp(argv[i], "-orient") != 0 &&
         strcmp(argv[i], "-pDelta") != 0 &&
         strcmp(argv[i], "-shearDist") != 0 &&
         strcmp(argv[i], "-doRayleigh") != 0 && strcmp(argv[i], "-mass") != 0) {
    numDir++;
    i++;
  }
  if (numDir != numMat) {
    opserr << "WARNING wrong number of directions specified\n";
    return TCL_ERROR;
  }
  // create the ID array to hold the direction IDs
  ID theDirIDs(numDir);
  // fill in the directions
  for (i = 0; i < numDir; i++) {
    if (Tcl_GetInt(interp, argv[argi], &dirID) != TCL_OK) {
      opserr << "WARNING invalid direction ID\n";
      opserr << "twoNodeLink element: " << tag << endln;
      return TCL_ERROR;
    }
    if (dirID < 1 || dirID > ndf) {
      opserr << "WARNING invalid direction ID: ";
      opserr << "dir = " << dirID << " > ndf = " << ndf;
      opserr << "\ntwoNodeLink element: " << tag << endln;
      return TCL_ERROR;
    }
    theDirIDs(i) = dirID - 1;
    argi++;
  }

  // check for optional arguments
  Vector x(0), y(0);
  for (i = argi; i < argc; i++) {
    if (strcmp(argv[i], "-orient") == 0) {
      int j = i + 1;
      int numOrient = 0;
      while (j < argc && strcmp(argv[j], "-pDelta") != 0 &&
             strcmp(argv[j], "-shearDist") != 0 &&
             strcmp(argv[j], "-doRayleigh") != 0 &&
             strcmp(argv[j], "-mass") != 0) {
        numOrient++;
        j++;
      }
      if (numOrient == 3) {
        y.resize(3);
        double value;
        // read the y values
        for (int j = 0; j < 3; j++) {
          if (Tcl_GetDouble(interp, argv[i + 1 + j], &value) != TCL_OK) {
            opserr << "WARNING invalid -orient value\n";
            opserr << "twoNodeLink element: " << tag << endln;
            return TCL_ERROR;
          } else {
            y(j) = value;
          }
        }
      } else if (numOrient == 6) {
        x.resize(3);
        y.resize(3);
        double value;
        // read the x values
        for (int j = 0; j < 3; j++) {
          if (Tcl_GetDouble(interp, argv[i + 1 + j], &value) != TCL_OK) {
            opserr << "WARNING invalid -orient value\n";
            opserr << "twoNodeLink element: " << tag << endln;
            return TCL_ERROR;
          } else {
            x(j) = value;
          }
        }
        // read the y values
        for (int j = 0; j < 3; j++) {
          if (Tcl_GetDouble(interp, argv[i + 4 + j], &value) != TCL_OK) {
            opserr << "WARNING invalid -orient value\n";
            opserr << "twoNodeLink element: " << tag << endln;
            return TCL_ERROR;
          } else {
            y(j) = value;
          }
        }
      } else {
        opserr << "WARNING insufficient arguments after -orient flag\n";
        opserr << "twoNodeLink element: " << tag << endln;
        return TCL_ERROR;
      }
    }
  }
  for (int i = argi; i < argc; i++) {
    if (i + 1 < argc && strcmp(argv[i], "-pDelta") == 0) {
      double Mr;
      Mratio.resize(4);
      if (ndm == 2) {
        Mratio.Zero();
        for (int j = 0; j < 2; j++) {
          if (Tcl_GetDouble(interp, argv[i + 1 + j], &Mr) != TCL_OK) {
            opserr << "WARNING invalid -pDelta value\n";
            opserr << "twoNodeLink element: " << tag << endln;
            return TCL_ERROR;
          }
          Mratio(2 + j) = Mr;
        }
      } else if (ndm == 3) {
        for (int j = 0; j < 4; j++) {
          if (Tcl_GetDouble(interp, argv[i + 1 + j], &Mr) != TCL_OK) {
            opserr << "WARNING invalid -pDelta value\n";
            opserr << "twoNodeLink element: " << tag << endln;
            return TCL_ERROR;
          }
          Mratio(j) = Mr;
        }
      }
    }
  }
  for (i = argi; i < argc; i++) {
    if (i + 1 < argc && strcmp(argv[i], "-shearDist") == 0) {
      double sDI;
      shearDistI.resize(2);
      if (ndm == 2) {
        if (Tcl_GetDouble(interp, argv[i + 1], &sDI) != TCL_OK) {
          opserr << "WARNING invalid -shearDist value\n";
          opserr << "twoNodeLink element: " << tag << endln;
          return TCL_ERROR;
        }
        shearDistI(0) = sDI;
        shearDistI(1) = 0.5;
      } else if (ndm == 3) {
        for (int j = 0; j < 2; j++) {
          if (Tcl_GetDouble(interp, argv[i + 1 + j], &sDI) != TCL_OK) {
            opserr << "WARNING invalid -shearDist value\n";
            opserr << "twoNodeLink element: " << tag << endln;
            return TCL_ERROR;
          }
          shearDistI(j) = sDI;
        }
      }
    }
  }
  for (i = argi; i < argc; i++) {
    if (strcmp(argv[i], "-doRayleigh") == 0)
      doRayleigh = 1;
  }
  for (i = argi; i < argc; i++) {
    if (i + 1 < argc && strcmp(argv[i], "-mass") == 0) {
      if (Tcl_GetDouble(interp, argv[i + 1], &mass) != TCL_OK) {
        opserr << "WARNING invalid -mass value\n";
        opserr << "twoNodeLink element: " << tag << endln;
        return TCL_ERROR;
      }
    }
  }

  // now create the twoNodeLink
  theElement = new TwoNodeLink(tag, ndm, iNode, jNode, theDirIDs, theMaterials,
                               y, x, Mratio, shearDistI, doRayleigh, mass);

  // cleanup dynamic memory
  if (theMaterials != 0)
    delete[] theMaterials;


  // then add the twoNodeLink to the domain
  if (builder->getDomain()->addElement(theElement) == false) {
    opserr << "WARNING could not add element to the domain\n";
    delete theElement;
    return TCL_ERROR;
  }

  return TCL_OK;
}


int
TclCommand_addTwoNodeLinkSection(ClientData clientData, Tcl_Interp *interp, int argc,
                               TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;
  constexpr static int eleArgStart = 1;


  int ndm = builder->getNDM();

  // Check the number of arguments is correct
  if (argc < 5) {
      opserr << "WARNING insufficient arguments\n";
      //               0       1           2      3     4     5
      opserr << "Want: element twoNodeLink eleTag iNode jNode secTag <-orient <x1 x2 x3> y1 y2 y3> <-pDelta Mratios> <-shearDist sDratios> <-doRayleigh> <-mass m>\n";
    return TCL_ERROR;
  }

  int argi = eleArgStart + 1;
  Vector Mratio(0), shearDistI(0), x(0), y(0);
  int doRayleigh = 0;
  double mass = 0.0;

  // get the id and end nodes
  int tag, iNode, jNode, stag;
  if (Tcl_GetInt(interp, argv[argi], &tag) != TCL_OK) {
    opserr << "WARNING invalid eleTag\n";
    return TCL_ERROR;
  }
  argi++;
  if (Tcl_GetInt(interp, argv[argi], &iNode) != TCL_OK) {
    opserr << "WARNING invalid iNode\n";
    return TCL_ERROR;
  }
  argi++;
  if (Tcl_GetInt(interp, argv[argi], &jNode) != TCL_OK) {
    opserr << "WARNING invalid jNode\n";
    return TCL_ERROR;
  }
  argi++;
  if (Tcl_GetInt(interp, argv[argi], &stag) != TCL_OK) {
    opserr << "WARNING invalid section tag\n";
    return TCL_ERROR;
  }
  SectionForceDeformation* theSection = builder->getTypedObject<SectionForceDeformation>(stag);
  if (theSection == nullptr)
    return TCL_ERROR;
  argi++;

  //
  int numdata;
  for (int i=argi; i<argc; i++) {

      if (strcmp(argv[i], "-orient") == 0) {
          if (OPS_GetNumRemainingInputArgs() < 3) {
              opserr << "WARNING: insufficient arguments after -orient\n";
              return 0;
          }
          numdata = 3;
          x.resize(3);
          if (OPS_GetDoubleInput(&numdata, &x(0)) < 0) {
              opserr << "WARNING: invalid -orient values\n";
              return 0;
          }
          if (OPS_GetNumRemainingInputArgs() < 3) {
              y = x;
              x = Vector();
              continue;
          }
          y.resize(3);
          if (OPS_GetDoubleInput(&numdata, &y(0)) < 0) {
              y = x;
              x = Vector();
              continue;
          }
      }

      else if (strcmp(argv[i], "-pDelta") == 0) {
          Mratio.resize(4);
          Mratio.Zero();
          numdata = 4;
          double* ptr = &Mratio(0);
          if (ndm == 2) {
              numdata = 2;
              ptr += 2;
          }
          if (OPS_GetNumRemainingInputArgs() < numdata) {
              opserr << "WARNING: insufficient data for -pDelta\n";
              return 0;
          }
          if (OPS_GetDoubleInput(&numdata, ptr) < 0) {
              opserr << "WARNING: invalid -pDelta value\n";
              return 0;
          }
      }
      else if (strcmp(argv[i], "-shearDist") == 0) {
          shearDistI.resize(2);
          numdata = 2;
          if (ndm == 2) {
              numdata = 1;
              shearDistI(1) = 0.5;
          }
          if (OPS_GetNumRemainingInputArgs() < numdata) {
              opserr << "WARNING: insufficient data for -shearDist\n";
              return 0;
          }
          if (OPS_GetDoubleInput(&numdata, &shearDistI(0)) < 0) {
              opserr << "WARNING: invalid -shearDist value\n";
              return 0;
          }
      }

      else if (strcmp(argv[i], "-doRayleigh") == 0) {
          doRayleigh = 1;
      }

      else if (strcmp(argv[i], "-mass") == 0) {
        if (argc < i + 2) {
          opserr << OpenSees::PromptValueError << "not enough arguments, expected -mass $mass\n";
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[i + 1], &mass) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid mass\n";
          return TCL_ERROR;
        }
        i += 1;
      }
  }

  // Check for optional arguments
  for (int i = argi; i < argc; i++) {
    if (strcmp(argv[i], "-orient") == 0) {
      int j = i + 1;
      int numOrient = 0;
      while (j < argc && strcmp(argv[j], "-pDelta") != 0 &&
             strcmp(argv[j], "-shearDist") != 0 &&
             strcmp(argv[j], "-doRayleigh") != 0 &&
             strcmp(argv[j], "-mass") != 0) {
        numOrient++;
        j++;
      }
      if (numOrient == 3) {
        y.resize(3);
        double value;
        // read the y values
        for (int j = 0; j < 3; j++) {
          if (Tcl_GetDouble(interp, argv[i + 1 + j], &value) != TCL_OK) {
            opserr << "WARNING invalid -orient value\n";
            opserr << "twoNodeLinkSection element: " << tag << endln;
            return TCL_ERROR;
          } else {
            y(j) = value;
          }
        }
      } else if (numOrient == 6) {
        x.resize(3);
        y.resize(3);
        double value;
        // read the x values
        for (int j = 0; j < 3; j++) {
          if (Tcl_GetDouble(interp, argv[i + 1 + j], &value) != TCL_OK) {
            opserr << "WARNING invalid -orient value\n";
            opserr << "twoNodeLinkSection element: " << tag << endln;
            return TCL_ERROR;
          } else {
            x(j) = value;
          }
        }
        // read the y values
        for (int j = 0; j < 3; j++) {
          if (Tcl_GetDouble(interp, argv[i + 4 + j], &value) != TCL_OK) {
            opserr << "WARNING invalid -orient value\n";
            opserr << "twoNodeLinkSection element: " << tag << endln;
            return TCL_ERROR;
          } else {
            y(j) = value;
          }
        }
      } else {
        opserr << "WARNING insufficient arguments after -orient flag\n";
        opserr << "twoNodeLinkSection element: " << tag << endln;
        return TCL_ERROR;
      }

    }

    else if (i + 1 < argc && strcmp(argv[i], "-pDelta") == 0) {
      double Mr;
      Mratio.resize(4);
      if (ndm == 2) {
        Mratio.Zero();
        for (int j = 0; j < 2; j++) {
          if (Tcl_GetDouble(interp, argv[i + 1 + j], &Mr) != TCL_OK) {
            opserr << "WARNING invalid -pDelta value\n";
            opserr << "twoNodeLinkSection element: " << tag << endln;
            return TCL_ERROR;
          }
          Mratio(2 + j) = Mr;
        }
      } else if (ndm == 3) {
        for (int j = 0; j < 4; j++) {
          if (Tcl_GetDouble(interp, argv[i + 1 + j], &Mr) != TCL_OK) {
            opserr << "WARNING invalid -pDelta value\n";
            opserr << "twoNodeLinkSection element: " << tag << endln;
            return TCL_ERROR;
          }
          Mratio(j) = Mr;
        }
      }
    }
    //
    else if (i + 1 < argc && strcmp(argv[i], "-shearDist") == 0) {
      double sDI;
      shearDistI.resize(2);
      if (ndm == 2) {
        if (Tcl_GetDouble(interp, argv[i + 1], &sDI) != TCL_OK) {
          opserr << "WARNING invalid -shearDist value\n";
          opserr << "twoNodeLinkSection element: " << tag << endln;
          return TCL_ERROR;
        }
        shearDistI(0) = sDI;
        shearDistI(1) = 0.5;
      } else if (ndm == 3) {
        for (int j = 0; j < 2; j++) {
          if (Tcl_GetDouble(interp, argv[i + 1 + j], &sDI) != TCL_OK) {
            opserr << "WARNING invalid -shearDist value\n";
            opserr << "twoNodeLinkSection element: " << tag << endln;
            return TCL_ERROR;
          }
          shearDistI(j) = sDI;
        }
      }
    }

    else if (strcmp(argv[i], "-doRayleigh") == 0)
      doRayleigh = 1;

    else if (i + 1 < argc && strcmp(argv[i], "-mass") == 0) {
      if (Tcl_GetDouble(interp, argv[i + 1], &mass) != TCL_OK) {
        opserr << "WARNING invalid -mass value\n";
        opserr << "twoNodeLinkSection element: " << tag << endln;
        return TCL_ERROR;
      }
    }
  }

  // now create the twoNodeLink
  Element* theElement = new TwoNodeLinkSection(tag, ndm, iNode, jNode,
        *theSection, y, x, Mratio, shearDistI, doRayleigh, mass);


  // then add the twoNodeLinkSection to the domain
  if (builder->getDomain()->addElement(theElement) == false) {
    opserr << "WARNING could not add element to the domain\n";
    opserr << "twoNodeLinkSection element: " << tag << endln;
    delete theElement;
    return TCL_ERROR;
  }

  // if get here we have successfully created the twoNodeLinkSection and added it to
  // the domain
  return TCL_OK;
}
