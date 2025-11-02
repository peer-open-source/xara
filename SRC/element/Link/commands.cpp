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
#include <float.h>
#include <string.h>
#include <cassert>
#include <Node.h>
#include <Domain.h>
#include <Logging.h>
#include <Parsing.h>
#include <ArgumentTracker.h>
#include <ID.h>
#include <Vector.h>
#include <UniaxialMaterial.h>
#include <FrameSection.h>
#include <ModelRegistry.h>
#include <TwoNodeLink.h>
#include <TwoNodeLinkSection.h>
#include <VectorND.h>

using OpenSees::VectorND;


int
TclCommand_addTwoNodeLink(ClientData clientData, Tcl_Interp *interp,
                          Tcl_Size argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  ModelRegistry *builder = (ModelRegistry*)clientData;

  // element TwoNodeLink eleTag iNode jNode -mat {m1 ...}|m1 m2 ...
  //                             -dir {d1 ...}|d1 d2 ...
  //   <-orient <x1 x2 x3> <y1 x2 y3> | -orient <x1 x2 x3> (NDM=1,2)
  //   | -orient <y1 y2 y3> (NDM=3 with x from nodes)>
  //   <-pDelta Mratios> <-shearDist sDratios> <-doRayleigh> <-mass m>


  const int ndm = builder->getNDM();

  if (argc < 8) {
    opserr << OpenSees::PromptValueError
           << "insufficient arguments for element " << argv[1] << "\n";
    return TCL_ERROR;
  }

  int tag, iNode, jNode;
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid element tag " << argv[2] << "\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid iNode " << argv[3] << "\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid jNode " << argv[4] << "\n";
    return TCL_ERROR;
  }

  int argi = 5;

  // Defaults
  Vector3D x{};
  Vector3D y{0, 1, 0};

  {
    Domain *theDomain = builder->getDomain();
    Node *ndI = theDomain->getNode(iNode);
    Node *ndJ = theDomain->getNode(jNode);
    if (ndI == nullptr || ndJ == nullptr) {
      opserr << OpenSees::PromptValueError << "invalid node tag(s)\n";
      return TCL_ERROR;
    }
    const Vector &end1 = ndI->getCrds();
    const Vector &end2 = ndJ->getCrds();
    for (int i=0; i<ndm; ++i)
      x(i) = end2(i) - end1(i);

    if (x.norm() < DBL_EPSILON) {
      x(0) = 1.0;
      x(1) = 0.0;
      x(2) = 0.0;
    }
    else {
      y(0) = -x(1);
      y(1) =  x(0);
      y(2) =  0.0;
    }
  }

  std::vector<UniaxialMaterial*> mats;
  std::vector<int> dirs_vec;
  Vector Mratio;
  Vector sDistI;
  int doRayleigh = 0;
  double mass = 0.0;

  //
  // Keyword loop
  //
  while (argi < argc) {
    if (strcmp(argv[argi], "-mat") == 0) {
      if (argi + 1 >= argc) {
        opserr << OpenSees::PromptValueError
               << "-mat requires arguments\n";
        return TCL_ERROR;
      }
      // Back-compat detection:
      // If argv[argi+2] exists and is an int, treat as unpacked sequence starting at argi+1.
      // Otherwise treat argv[argi+1] as a Tcl list.
      bool unpacked = false;
      if (argi + 2 < argc) {
        int tmp;
        if (Tcl_GetInt(interp, argv[argi+2], &tmp) == TCL_OK)
          unpacked = true;
      }

      if (unpacked) {
        int j = argi + 1;
        int mt;
        while (j < argc && Tcl_GetInt(interp, argv[j], &mt) == TCL_OK) {
          UniaxialMaterial *umat = builder->getTypedObject<UniaxialMaterial>(mt);
          if (umat == nullptr) {
            return TCL_ERROR;
          }
          mats.push_back(umat);
          j++;
        }
        if (mats.empty()) {
          opserr << OpenSees::PromptValueError
                 << "no materials provided after -mat\n";
          return TCL_ERROR;
        }
        argi = j; // advance past all parsed ints
      }
      else {
        int list_argc;
        TCL_Char **list_argv;
        if (Tcl_SplitList(interp, argv[argi+1], &list_argc, &list_argv) != TCL_OK || list_argc < 1) {
          opserr << OpenSees::PromptValueError << "invalid material list after -mat\n";
          return TCL_ERROR;
        }
        for (int i=0; i<list_argc; ++i) {
          int mt;
          if (Tcl_GetInt(interp, list_argv[i], &mt) != TCL_OK) {
            opserr << OpenSees::PromptValueError << "invalid material tag "
                   << list_argv[i]
                   << "\n";
            Tcl_Free((char*)list_argv);
            return TCL_ERROR;
          }
          UniaxialMaterial *umat = builder->getTypedObject<UniaxialMaterial>(mt);
          if (umat == nullptr) {
            Tcl_Free((char*)list_argv);
            return TCL_ERROR;
          }
          mats.push_back(umat);
        }
        Tcl_Free((char*)list_argv);
        argi += 2; // "-mat" + list token
      }
    }

    else if (strcmp(argv[argi], "-dir") == 0 || strcmp(argv[argi], "-dof") == 0) {
      if (argi + 1 >= argc) {
        opserr << OpenSees::PromptValueError
               << "-dir requires arguments\n";
        return TCL_ERROR;
      }
      bool unpacked = false;
      if (argi + 2 < argc) {
        int tmp;
        if (Tcl_GetInt(interp, argv[argi+2], &tmp) == TCL_OK)
          unpacked = true;
      }

      if (unpacked) {
        int j = argi + 1;
        int d;
        while (j < argc && Tcl_GetInt(interp, argv[j], &d) == TCL_OK) {
          dirs_vec.push_back(d);
          j++;
        }
        if (dirs_vec.empty()) {
          opserr << OpenSees::PromptValueError
                 << "no directions provided after -dir\n";
          return TCL_ERROR;
        }
        argi = j; // advance past all parsed ints
      } else {
        int list_argc;
        TCL_Char **list_argv;
        if (Tcl_SplitList(interp, argv[argi+1], &list_argc, &list_argv) != TCL_OK || list_argc < 1) {
          opserr << OpenSees::PromptValueError
                 << "invalid direction list after -dir\n";
          return TCL_ERROR;
        }
        for (int i=0; i<list_argc; ++i) {
          int d;
          if (Tcl_GetInt(interp, list_argv[i], &d) != TCL_OK) {
            opserr << OpenSees::PromptValueError
                   << "invalid direction " << list_argv[i] << "\n";
            Tcl_Free((char*)list_argv);
            return TCL_ERROR;
          }
          dirs_vec.push_back(d);
        }
        Tcl_Free((char*)list_argv);
        argi += 2; // "-dir" + list token
      }
    }

    else if (strcmp(argv[argi], "-orient") == 0) {
      argi++;
      int remaining = argc - argi;
      if (remaining < 3) {
        opserr << OpenSees::PromptValueError
               << "insufficient arguments after -orient\n";
        return TCL_ERROR;
      }
      if (remaining >= 6) {
        for (int k=0; k<3; ++k) {
          if (Tcl_GetDouble(interp, argv[argi++], &x(k)) != TCL_OK) {
            opserr << OpenSees::PromptValueError
                   << "invalid -orient x component\n";
            return TCL_ERROR;
          }
        }
        for (int k=0; k<3; ++k) {
          if (Tcl_GetDouble(interp, argv[argi++], &y(k)) != TCL_OK) {
            opserr << OpenSees::PromptValueError
                   << "invalid -orient y component\n";
            return TCL_ERROR;
          }
        }
      } else {
        if (ndm == 1 || ndm == 2) {
          for (int k=0; k<3; ++k) {
            if (Tcl_GetDouble(interp, argv[argi++], &x(k)) != TCL_OK) {
              opserr << OpenSees::PromptValueError << "invalid -orient x component\n";
              return TCL_ERROR;
            }
          }
        } else {
          for (int k=0; k<3; ++k) {
            if (Tcl_GetDouble(interp, argv[argi++], &y(k)) != TCL_OK) {
              opserr << OpenSees::PromptValueError << "invalid -orient y component\n";
              return TCL_ERROR;
            }
          }
        }
      }
    }

    else if (strcmp(argv[argi], "-pDelta") == 0) {
      argi++;
      Mratio.resize(4); Mratio.Zero();
      if (ndm == 2) {
        if (argi + 2 > argc) {
          opserr << OpenSees::PromptValueError << "insufficient data for -pDelta\n";
          return TCL_ERROR;
        }
        for (int k=2; k<4; ++k) {
          if (Tcl_GetDouble(interp, argv[argi++], &Mratio(k)) != TCL_OK) {
            opserr << OpenSees::PromptValueError << "invalid -pDelta value\n";
            return TCL_ERROR;
          }
        }
      } else {
        if (argi + 4 > argc) {
          opserr << OpenSees::PromptValueError << "insufficient data for -pDelta\n";
          return TCL_ERROR;
        }
        for (int k=0; k<4; ++k) {
          if (Tcl_GetDouble(interp, argv[argi++], &Mratio(k)) != TCL_OK) {
            opserr << OpenSees::PromptValueError << "invalid -pDelta value\n";
            return TCL_ERROR;
          }
        }
      }
    }

    else if (strcmp(argv[argi], "-shearDist") == 0) {
      argi++;
      sDistI.resize(2); sDistI(0)=0.0; sDistI(1)=0.5;
      if (ndm == 2) {
        if (argi + 1 > argc) {
          opserr << OpenSees::PromptValueError << "insufficient data for -shearDist\n";
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[argi++], &sDistI(0)) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid -shearDist value\n";
          return TCL_ERROR;
        }
      } else {
        if (argi + 2 > argc) {
          opserr << OpenSees::PromptValueError << "insufficient data for -shearDist\n";
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[argi++], &sDistI(0)) != TCL_OK ||
            Tcl_GetDouble(interp, argv[argi++], &sDistI(1)) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid -shearDist value\n";
          return TCL_ERROR;
        }
      }
    }

    else if (strcmp(argv[argi], "-doRayleigh") == 0) {
      argi++;
      doRayleigh = 1;
    }

    else if (strcmp(argv[argi], "-mass") == 0) {
      argi++;
      if (argi >= argc) {
        opserr << OpenSees::PromptValueError << "insufficient mass value\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[argi++], &mass) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid -mass value\n";
        return TCL_ERROR;
      }
    }

    else {
      opserr << OpenSees::PromptParseError << "unexpected option " << argv[argi] << "\n";
      return TCL_ERROR;
    }
  }

  //
  // Validate -mat / -dir
  //
  if (mats.empty()) {
    opserr << OpenSees::PromptValueError << "no materials provided\n";
    return TCL_ERROR;
  }
  if (dirs_vec.size() != mats.size()) {
    opserr << OpenSees::PromptValueError
           << "wrong number of directions, expected " << static_cast<int>(mats.size()) << "\n";
    return TCL_ERROR;
  }

  ID dirs((int)mats.size());
  for (int i=0; i<(int)mats.size(); ++i)
    dirs(i) = dirs_vec[i] - 1; // 1-based to 0-based

  //
  // Create and add element
  //
  Element *theEle = new TwoNodeLink(tag, ndm, iNode, jNode,
                                    dirs, &mats[0],
                                    y, x,
                                    Mratio, sDistI,
                                    doRayleigh, mass);

  if (!builder->getDomain()->addElement(theEle)) {
    opserr << OpenSees::PromptValueError
           << "could not add element to the domain\n";
    delete theEle;
    return TCL_ERROR;
  }

  return TCL_OK;
}


int
TclCommand_addTwoNodeLinkSection(ClientData clientData, Tcl_Interp *interp,
                                  Tcl_Size argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  ModelRegistry *builder = (ModelRegistry*)clientData;

  // Usage:
  // element TwoNodeLinkSection eleTag iNode jNode secTag
  //   <-orient <x1 x2 x3> <y1 y2 y3> | -orient <x1 x2 x3> (NDM=1,2)
  //   | -orient <y1 y2 y3> (NDM=3 with x from nodes)>
  //   <-pDelta Mratios> <-shearDist sDratios> <-doRayleigh> <-mass m>

  const int ndm = builder->getNDM();

  if (argc < 6) {
    opserr << OpenSees::PromptValueError
           << "insufficient arguments for element " << argv[1] << "\n";
    return TCL_ERROR;
  }

  int tag, iNode, jNode, secTag;
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid element tag " << argv[2] << "\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid iNode " << argv[3] << "\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid jNode " << argv[4] << "\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[5], &secTag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid section tag " << argv[5] << "\n";
    return TCL_ERROR;
  }

  // Fetch section
  SectionForceDeformation *theSection =
      builder->getTypedObject<FrameSection>(secTag);
  if (theSection == nullptr) {
    opserr << OpenSees::PromptValueError
           << "no SectionForceDeformation with tag " << secTag << "\n";
    return TCL_ERROR;
  }

  // Default orientation: x = end2 - end1; y = (0,1,0)
  Vector x(3); x.Zero();
  Vector y(3); y(0)=0.0; y(1)=1.0; y(2)=0.0;

  {
    Domain *theDomain = builder->getDomain();
    Node *ndI = theDomain->getNode(iNode);
    Node *ndJ = theDomain->getNode(jNode);
    if (ndI == nullptr || ndJ == nullptr) {
      opserr << OpenSees::PromptValueError << "invalid node tag(s)\n";
      return TCL_ERROR;
    }
    const Vector &end1 = ndI->getCrds();
    const Vector &end2 = ndJ->getCrds();
    for (int i=0; i<ndm; ++i)
      x(i) = end2(i) - end1(i);
  }

  // Options
  Vector Mratio;   // will be resized on demand
  Vector sDistI;   // will be resized on demand
  int doRayleigh = 0;
  double mass = 0.0;

  // If only required args are present, construct minimal element
  if (argc == 6) {
    Element *theEle = new TwoNodeLinkSection(tag, ndm, iNode, jNode, *theSection);
    if (!builder->getDomain()->addElement(theEle)) {
      opserr << OpenSees::PromptValueError << "could not add element to the domain\n";
      delete theEle;
      return TCL_ERROR;
    }
    return TCL_OK;
  }

  // Parse keyword options
  int i = 6;
  while (i < argc) {
    if (strcmp(argv[i], "-orient") == 0) {
      i++;
      // Read 6 doubles (x then y) if available; otherwise follow NDM logic
      int remaining = argc - i;
      if (remaining >= 6) {
        // x
        for (int k=0; k<3; ++k) {
          if (Tcl_GetDouble(interp, argv[i++], &x(k)) != TCL_OK) {
            opserr << OpenSees::PromptValueError
                   << "invalid orient x component"
                   << OpenSees::SignalMessageEnd;
            return TCL_ERROR;
          }
        }
        // y
        for (int k=0; k<3; ++k) {
          if (Tcl_GetDouble(interp, argv[i++], &y[k]) != TCL_OK) {
            opserr << OpenSees::PromptValueError
                   << "invalid orient y component"
                   << OpenSees::SignalMessageEnd;
            return TCL_ERROR;
          }
        }
      } else {
        // Only one vector provided: if NDM=1,2 treat as x; if NDM=3 treat as y.
        if (remaining < 3) {
          opserr << OpenSees::PromptValueError
                 << "insufficient arguments for orient"
                 << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        }
        if (ndm == 1 || ndm == 2) {
          for (int k=0; k<3; ++k) {
            if (Tcl_GetDouble(interp, argv[i++], &x[k]) != TCL_OK) {
              opserr << OpenSees::PromptValueError
                     << "invalid orient x component"
                     << OpenSees::SignalMessageEnd;
              return TCL_ERROR;
            }
          }
        } else {
          for (int k=0; k<3; ++k) {
            if (Tcl_GetDouble(interp, argv[i++], &y(k)) != TCL_OK) {
              opserr << OpenSees::PromptValueError << "invalid -orient y component\n";
              return TCL_ERROR;
            }
          }
        }
      }
    }

    else if (strcmp(argv[i], "-pDelta") == 0) {
      i++;
      Mratio.resize(4); Mratio.Zero();
      if (ndm == 2) {
        // Read 2 values into the last two entries (matches your pointer offset logic)
        if (i + 2 > argc) {
          opserr << OpenSees::PromptValueError << "insufficient data for -pDelta\n";
          return TCL_ERROR;
        }
        for (int k=2; k<4; ++k) {
          if (Tcl_GetDouble(interp, argv[i++], &Mratio(k)) != TCL_OK) {
            opserr << OpenSees::PromptValueError << "invalid -pDelta value\n";
            return TCL_ERROR;
          }
        }
      } else {
        if (i + 4 > argc) {
          opserr << OpenSees::PromptValueError << "insufficient data for -pDelta\n";
          return TCL_ERROR;
        }
        for (int k=0; k<4; ++k) {
          if (Tcl_GetDouble(interp, argv[i++], &Mratio(k)) != TCL_OK) {
            opserr << OpenSees::PromptValueError << "invalid -pDelta value\n";
            return TCL_ERROR;
          }
        }
      }
    }

    else if (strcmp(argv[i], "-shearDist") == 0) {
      i++;
      sDistI.resize(2); sDistI(0) = 0.0; sDistI(1) = 0.5; // default for 2D path matches your code
      if (ndm == 2) {
        if (i + 1 > argc) {
          opserr << OpenSees::PromptValueError 
                 << "insufficient data for -shearDist"
                 << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[i++], &sDistI(0)) != TCL_OK) {
          opserr << OpenSees::PromptValueError 
                 << "invalid shearDist value"
                 << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        }
        // sDistI(1) stays 0.5 as in your original
      } else {
        if (i + 2 > argc) {
          opserr << OpenSees::PromptValueError 
                 << "insufficient data for shearDist"
                 << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[i++], &sDistI(0)) != TCL_OK ||
            Tcl_GetDouble(interp, argv[i++], &sDistI(1)) != TCL_OK) {
          opserr << OpenSees::PromptValueError
                 << "invalid shearDist"
                 << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        }
      }
    }

    else if (strcmp(argv[i], "-doRayleigh") == 0) {
      i++;
      doRayleigh = 1;
    }

    else if (strcmp(argv[i], "-mass") == 0) {
      i++;
      if (i >= argc) {
        opserr << OpenSees::PromptValueError
               << "insufficient mass values"
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i++], &mass) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "invalid -mass value"
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      }
    }

    else {
      opserr << OpenSees::PromptParseError
             << "unexpected option " << argv[i]
             << OpenSees::SignalMessageEnd;
      return TCL_ERROR;
    }
  }

  // Create element (full signature)
  Element *theEle = new TwoNodeLinkSection(tag, ndm, iNode, jNode,
                                           *theSection, y, x, Mratio, sDistI,
                                           doRayleigh, mass);

  if (!builder->getDomain()->addElement(theEle)) {
    opserr << OpenSees::PromptValueError
           << "could not add element to the domain"
           << OpenSees::SignalMessageEnd;
    delete theEle;
    return TCL_ERROR;
  }

  return TCL_OK;
}
