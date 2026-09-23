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
//
#include <string.h>
#ifdef _MSC_VER 
#  include <string.h>
#  define strcasecmp _stricmp
#else
#  include <strings.h>
#endif
#include <set>
#include <array>
#include <cstdlib>

#include <Parsing.h>
#include <Logging.h>
#include <Domain.h>
#include <ModelRegistry.h>
#include <ArgumentTracker.h>

#include <Brick.h>
#include <Brick02.h>
// #define XARA_HAVE_H8E12
#ifdef XARA_HAVE_H8E12
# include <H8E12.h>
#endif
#include <BbarBrick.h>
#include <BbarBrickWithSensitivity.h>
#include <Twenty_Node_Brick.h>
#include <FourNodeTetrahedron.h>
#include <TenNodeTetrahedron.h>
#include <element/community/UWelements/SSPbrick.h>

class SolidElement {
public:
  enum class ElementType {
    Brick,
    Brick02,
    BbarBrick,
    SSPBrick,
    BbarBrickWithSensitivity,
    H8E12,
    TwentyNodeBrick,
    TwentySevenNodeBrick,
    Tet4,
    Tet10
  };

  enum class CellType {
    H8,
    H20,
    H27,
    T4,
    T10
  };

  SolidElement(ElementType type) : element_type(type) {}

  ElementType type() const { return element_type; }

  CellType cell() const {
    switch (element_type) {
      case ElementType::Brick:
      case ElementType::Brick02:
      case ElementType::BbarBrick:
      case ElementType::BbarBrickWithSensitivity:
      case ElementType::H8E12:
      case ElementType::SSPBrick:
        return CellType::H8;
      case ElementType::TwentyNodeBrick:
        return CellType::H20;
      case ElementType::TwentySevenNodeBrick:
        return CellType::H27;
      case ElementType::Tet4:
        return CellType::T4;
      case ElementType::Tet10:
        return CellType::T10;
    }
  }

  std::size_t nodes() const {
    switch (cell()) {
      case CellType::H8:  return  8;
      case CellType::H20: return 20;
      case CellType::H27: return 27;
      case CellType::T4:  return  4;
      case CellType::T10: return 10;
    }
  }

private:
  const ElementType element_type;
};


template <std::size_t nen>
static int
CreateSolidElement(ClientData clientData, 
                   Tcl_Interp *interp,
                   ArgSize argc,
                   TCL_Char **const argv,
                   const SolidElement& solid_element)
{

  ModelRegistry* builder = (ModelRegistry*)clientData;
  Domain* theTclDomain = builder->getDomain();

  enum class Position : int {
    Material,
    EndRequired,
      B1, B2, B3,
    End
  };
  ArgumentTracker<Position> tracker;
  std::set<int> positional;


  // get the element tag
  int tag;
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid element tag" << "\n";
    return TCL_ERROR;
  }

  // Parse nodes
  std::array<int, nen> node_tags{};
  int argi = 3;
  {
    int list_argc;
    TCL_Char **list_argv;
    if (Tcl_SplitList(interp, argv[argi], &list_argc, &list_argv) == TCL_OK && list_argc == nen) {
      for (int i = 0; i < list_argc; ++i) {
        int node;
        if (Tcl_GetInt(interp, list_argv[i], &node) != TCL_OK) {
          opserr << OpenSees::PromptValueError 
                 << "invalid node " << list_argv[i]
                 << "\n";
          return TCL_ERROR;
        }
        node_tags[i] = node;
      }
      Tcl_Free((char *)list_argv);
      argi += 1;
    }
    else {
      if (argi + nen > argc) {
        opserr << OpenSees::PromptValueError 
               << "expected " << nen << " nodes for element type " << argv[1] 
               << "\n";
        return TCL_ERROR;
      }
      for (int i=0; i<nen; i++) {
        int node;
        if (Tcl_GetInt(interp, argv[argi++], &node) != TCL_OK) {
          opserr << OpenSees::PromptValueError 
                 << "invalid node tag " << argv[argi-1] 
                 << "\n";
          return TCL_ERROR;
        }
        node_tags[i] = node;
      }
    }
  }

  //
  // Positional/Keyword arguments
  //

  int matID;
  std::array<double, 3> body_force = {0.0, 0.0, 0.0};
  // keyword pass
  for (int i=argi; i<argc; i++) {
    if (strcasecmp(argv[i], "-material") == 0) {
      if (i+1 >= argc) {
        opserr << OpenSees::PromptValueError << "expected material tag after -material\n";
        return TCL_ERROR;
      }
      if (Tcl_GetInt(interp, argv[i+1], &matID) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid material tag " << argv[i+1] << "\n";
        return TCL_ERROR;
      }
      i += 1;
      tracker.consume(Position::Material);
    }
    else if (strcasecmp(argv[i], "-b") == 0) {
      if (i+1 >= argc) {
        opserr << OpenSees::PromptValueError << "expected 3 body force components after -b\n";
        return TCL_ERROR;
      }
      int argc_b;
      TCL_Char **argv_b;
      if (Tcl_SplitList(interp, argv[i+1], &argc_b, &argv_b) != TCL_OK || argc_b != 3) {
        opserr << OpenSees::PromptValueError << "expected 3 body force components after -b\n";
        return TCL_ERROR;
      }
      for (int j=0; j<3; j++) {
        if (Tcl_GetDouble(interp, argv_b[j], &body_force[j]) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid body force component " << argv_b[j] << "\n";
          return TCL_ERROR;
        }
      }
      Tcl_Free((char *)argv_b);
      i += 1;
      tracker.consume(Position::B1);
      tracker.consume(Position::B2);
      tracker.consume(Position::B3);
    }
    else {
      positional.insert(i);
    }
  }

  //
  // Positional arguments
  //
  for (int i: positional) {
    if (tracker.current() == Position::EndRequired)
      continue;
    switch (tracker.current()) {
      case Position::Material:
        if (Tcl_GetInt(interp, argv[i], &matID) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid material tag " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::Material);
        break;
      case Position::B1:
        if (Tcl_GetDouble(interp, argv[i], &body_force[0]) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid body force component " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::B1);
        break;
      case Position::B2:
        if (Tcl_GetDouble(interp, argv[i], &body_force[1]) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid body force component " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::B2);
        break;
      case Position::B3:
        if (Tcl_GetDouble(interp, argv[i], &body_force[2]) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid body force component " << argv[i] << "\n";
          return TCL_ERROR;
        }
        tracker.consume(Position::B3);
        break;
      default:
        break;
    }
  }

  //
  // Finalize/validate arguments
  //
  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == nullptr)
    return TCL_ERROR;

  double b1 = body_force[0],
         b2 = body_force[1],
         b3 = body_force[2];
  // now create the element and add it to the Domain
  Element *theBrick = nullptr;
  if constexpr (nen == 8) {
    switch (solid_element.type()) {
      case SolidElement::ElementType::Brick:
        theBrick = new Brick(tag, node_tags, *theMaterial, b1, b2, b3);
        break;
      case SolidElement::ElementType::Brick02:
        theBrick = new Brick02(tag, node_tags, *theMaterial, b1, b2, b3);
        break;
      case SolidElement::ElementType::BbarBrick:
        theBrick = new BbarBrick(tag, node_tags, *theMaterial, b1, b2, b3);
        break;
      case SolidElement::ElementType::BbarBrickWithSensitivity:
        theBrick = new BbarBrickWithSensitivity(tag, node_tags, *theMaterial, b1, b2, b3);
        break;
      case SolidElement::ElementType::SSPBrick:
        theBrick = new SSPbrick(tag, node_tags, *theMaterial, b1, b2, b3);
        break;
#ifdef XARA_HAVE_H8E12
      case SolidElement::ElementType::H8E12:
        theBrick = new H8E12(tag, node_tags, *theMaterial, b1, b2, b3);
        break;
#endif
      default:
        opserr << OpenSees::PromptValueError << "invalid element type for 8-node solid element\n";
        return TCL_ERROR;
    }
  }
  else if constexpr (nen == 20) {
    theBrick = new Twenty_Node_Brick(tag, node_tags, *theMaterial, b1, b2, b3);
  }
  else if constexpr (nen == 27) {
    return TCL_ERROR;
  }
  else if constexpr (nen == 4) {
    theBrick = new FourNodeTetrahedron(tag, node_tags, *theMaterial, b1, b2, b3);
  }
  else if constexpr (nen == 10) {
    theBrick = new TenNodeTetrahedron(tag, node_tags, *theMaterial, b1, b2, b3);
  }
  else {
    static_assert(false, "invalid number of nodes for solid element");
  }


  if (theTclDomain->addElement(theBrick) == false) {
    opserr << "WARNING could not add element to the domain\n";
    delete theBrick;
    return TCL_ERROR;
  }

  return TCL_OK;
}


int
XaraElemCmd_SolidElement(ClientData clientData, 
                         Tcl_Interp *interp,
                         ArgSize argc,
                         TCL_Char **const argv)
{
  //
  SolidElement::ElementType element_type;
  if ((strcasecmp(argv[1], "Brick") == 0) || 
      (strcasecmp(argv[1], "stdBrick") == 0) || 
      (strcasecmp(argv[1], "H8") == 0)) {
    element_type = SolidElement::ElementType::Brick;
  }
  else if (strcasecmp(argv[1], "Brick02") == 0) {
    element_type = SolidElement::ElementType::Brick02;
  }
  else if (strcasecmp(argv[1], "bbarBrick") == 0) {
    element_type = SolidElement::ElementType::BbarBrick;
  }
  else if (strcasecmp(argv[1], "bbarBrickWithSensitivity") == 0) {
    element_type = SolidElement::ElementType::BbarBrickWithSensitivity;
  }
  else if (strcasecmp(argv[1], "SSPbrick") == 0) {
    element_type = SolidElement::ElementType::SSPBrick;
  }
#ifdef XARA_HAVE_H8E12
  else if (strcasecmp(argv[1], "H8E12") == 0) {
    element_type = SolidElement::ElementType::H8E12;
  }
#endif
  else if ((strcasecmp(argv[1], "20NodeBrick") == 0) || 
           (strcasecmp(argv[1], "H20") == 0)) {
    element_type = SolidElement::ElementType::TwentyNodeBrick;
  }
  else if ((strcasecmp(argv[1], "27NodeBrick") == 0) || 
           (strcasecmp(argv[1], "H27") == 0)) {
    element_type = SolidElement::ElementType::TwentySevenNodeBrick;
  }
  else if ((strcasecmp(argv[1], "Tet4") == 0) || 
           (strcasecmp(argv[1], "FourNodeTetrahedron") == 0) ||
           (strcasecmp(argv[1], "T4") == 0)) {
    element_type = SolidElement::ElementType::Tet4;
  }
  else if ((strcasecmp(argv[1], "Tet10") == 0) || 
           (strcasecmp(argv[1], "TenNodeTetrahedron") == 0) ||
           (strcasecmp(argv[1], "T10") == 0)) {
    element_type = SolidElement::ElementType::Tet10;
  }
  else {
    opserr << OpenSees::PromptValueError << "invalid element type\n";
    return TCL_ERROR;
  }
  SolidElement solid_element(element_type);

  switch (solid_element.nodes()) {
    case 8:
      return CreateSolidElement<8>(clientData, interp, argc, argv, solid_element);
    case 20:
      return CreateSolidElement<20>(clientData, interp, argc, argv, solid_element);
    case 27:
      return CreateSolidElement<27>(clientData, interp, argc, argv, solid_element);
    case 4:
      return CreateSolidElement<4>(clientData, interp, argc, argv, solid_element);
    case 10:
      return CreateSolidElement<10>(clientData, interp, argc, argv, solid_element);
    default:
      opserr << OpenSees::PromptValueError << "invalid number of nodes for element type\n";
      return TCL_ERROR;
  }
}

//
// Description: This file contains the implementation of
//    XaraElemCmd_H8UP() ,
//    TclBasicBuilder_addTwentyEightNodeBrickUP(),
//    TclBasicBuilder_addBBarBrickUP()
//
//
// Jinchi Lu and Zhaohui Yang (May 2004)
//
#include <BrickUP.h>
#include <BBarBrickUP.h>
#include <Twenty_Eight_Node_BrickUP.h>

int
XaraElemCmd_H8UP(ClientData clientData, Tcl_Interp *interp, ArgSize argc,
                           TCL_Char ** const argv)
{
  ModelRegistry *builder = (ModelRegistry*)clientData;
  Domain* theTclDomain = builder->getDomain();

  if (builder == 0 || clientData == 0) {
    opserr << "WARNING builder has been destroyed\n";
    return TCL_ERROR;
  }

  if (builder->getNDM() != 3 || builder->getNDF() != 4) {
    opserr << "WARNING -- model dimensions and/or nodal DOF not compatible "
              "with QuadUP element\n";
    return TCL_ERROR;
  }

  // check the number of arguments is correct
  int argStart = 2;

  if ((argc - argStart) < 15) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element brickUP eleTag? N1? N2? N3? N4? N5? N6? N7? N8? "
              "matTag? bulk? rhof? perm_x? perm_y? perm_z? <b1? b2? b3?>\n";
    return TCL_ERROR;
  }

  // get the id and end nodes
  int brickUPId, Nod[8], matID;
  double bk, r, perm1, perm2, perm3;
  double b1 = 0.0;
  double b2 = 0.0;
  double b3 = 0.0;

  if (Tcl_GetInt(interp, argv[argStart], &brickUPId) != TCL_OK) {
    opserr << "WARNING invalid brickUP eleTag" << "\n";
    return TCL_ERROR;
  }

  for (int i = 0; i < 8; i++)
    if (Tcl_GetInt(interp, argv[1 + argStart + i], &(Nod[i])) != TCL_OK) {
      opserr << "WARNING invalid Node number\n";
      return TCL_ERROR;
    }

  if (Tcl_GetInt(interp, argv[9 + argStart], &matID) != TCL_OK) {
    opserr << "WARNING invalid matID\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[10 + argStart], &bk) != TCL_OK) {
    opserr << "WARNING invalid fluid bulk modulus\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[11 + argStart], &r) != TCL_OK) {
    opserr << "WARNING invalid fluid mass density\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[12 + argStart], &perm1) != TCL_OK) {
    opserr << "WARNING invalid permeability_x\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[13 + argStart], &perm2) != TCL_OK) {
    opserr << "WARNING invalid permeability_y\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[14 + argStart], &perm3) != TCL_OK) {
    opserr << "WARNING invalid permeability_z\n";
    return TCL_ERROR;
  }

  if ((argc - argStart) >= 16) {
    if (Tcl_GetDouble(interp, argv[15 + argStart], &b1) != TCL_OK) {
      opserr << "WARNING invalid b1\n";
      return TCL_ERROR;
    }
  }
  if ((argc - argStart) >= 17) {
    if (Tcl_GetDouble(interp, argv[16 + argStart], &b2) != TCL_OK) {
      opserr << "WARNING invalid b2\n";
      return TCL_ERROR;
    }
  }
  if ((argc - argStart) >= 18) {
    if (Tcl_GetDouble(interp, argv[17 + argStart], &b3) != TCL_OK) {
      opserr << "WARNING invalid b3\n";
      return TCL_ERROR;
    }
  }

  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == nullptr) {
    return TCL_ERROR;
  }

  // now create the brickUP and add it to the Domain
  BrickUP *theBrickUP = new BrickUP(
      brickUPId, Nod[0], Nod[1], Nod[2], Nod[3], Nod[4], Nod[5], Nod[6], Nod[7],
      *theMaterial, bk, r, perm1, perm2, perm3, b1, b2, b3);


  if (theTclDomain->addElement(theBrickUP) == false) {
    opserr << "WARNING could not add element to the domain\n";
    delete theBrickUP;
    return TCL_ERROR;
  }

  return TCL_OK;
}

int
TclBasicBuilder_addTwentyEightNodeBrickUP(ClientData clientData, Tcl_Interp *interp,
                                          int argc, TCL_Char ** const argv)

{
  ModelRegistry *builder = (ModelRegistry*)clientData;
  Domain* theTclDomain = builder->getDomain();


  if (builder->getNDM() != 3) {
    opserr << "WARNING -- model dimensions and/or nodal DOF not compatible "
              "with 20_8_BrickUP element\n";
    return TCL_ERROR;
  }

  // check the number of arguments is correct
  int argStart = 2;

  if ((argc - argStart) < 27) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element 20_8_BrickUP eleTag? N1? N2? N3? N4? N5? N6? N7? "
              "N8? N9? N10? N11? N12? N13? N14? N15? N16? N17? N18? N19? N20? "
              "matTag? bulk? rhof? perm_x? perm_y? perm_z? <b1? b2? b3?>\n";
    return TCL_ERROR;
  }

  // get the id and end nodes
  int brickUPId, Nod[20], matID;
  double bk, r, perm1, perm2, perm3;
  double b1 = 0.0;
  double b2 = 0.0;
  double b3 = 0.0;

  if (Tcl_GetInt(interp, argv[argStart], &brickUPId) != TCL_OK) {
    opserr << "WARNING invalid 20_8_BrickUP eleTag" << "\n";
    return TCL_ERROR;
  }

  for (int i = 0; i < 20; i++)
    if (Tcl_GetInt(interp, argv[1 + argStart + i], &(Nod[i])) != TCL_OK) {
      opserr << "WARNING invalid Node number\n";
      return TCL_ERROR;
    }

  if (Tcl_GetInt(interp, argv[21 + argStart], &matID) != TCL_OK) {
    opserr << "WARNING invalid matID\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[22 + argStart], &bk) != TCL_OK) {
    opserr << "WARNING invalid fluid bulk modulus\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[23 + argStart], &r) != TCL_OK) {
    opserr << "WARNING invalid fluid mass density\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[24 + argStart], &perm1) != TCL_OK) {
    opserr << "WARNING invalid permeability_x\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[25 + argStart], &perm2) != TCL_OK) {
    opserr << "WARNING invalid permeability_y\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[26 + argStart], &perm3) != TCL_OK) {
    opserr << "WARNING invalid permeability_z\n";
    return TCL_ERROR;
  }

  if ((argc - argStart) >= 28) {
    if (Tcl_GetDouble(interp, argv[27 + argStart], &b1) != TCL_OK) {
      opserr << "WARNING invalid b1\n";
      return TCL_ERROR;
    }
  }
  if ((argc - argStart) >= 29) {
    if (Tcl_GetDouble(interp, argv[28 + argStart], &b2) != TCL_OK) {
      opserr << "WARNING invalid b2\n";
      return TCL_ERROR;
    }
  }
  if ((argc - argStart) >= 30) {
    if (Tcl_GetDouble(interp, argv[29 + argStart], &b3) != TCL_OK) {
      opserr << "WARNING invalid b3\n";
      return TCL_ERROR;
    }
  }

  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);

  if (theMaterial == nullptr) {
    return TCL_ERROR;
  }

  // now create the brickUP and add it to the Domain
  TwentyEightNodeBrickUP *theTwentyEightNodeBrickUP =
      new TwentyEightNodeBrickUP(
          brickUPId, Nod[0], Nod[1], Nod[2], Nod[3], Nod[4], Nod[5], Nod[6],
          Nod[7], Nod[8], Nod[9], Nod[10], Nod[11], Nod[12], Nod[13], Nod[14],
          Nod[15], Nod[16], Nod[17], Nod[18], Nod[19], *theMaterial, bk, r,
          perm1, perm2, perm3, b1, b2, b3);

  if (theTclDomain->addElement(theTwentyEightNodeBrickUP) == false) {
    opserr << "WARNING could not add element to the domain\n";
    delete theTwentyEightNodeBrickUP;
    return TCL_ERROR;
  }

  // if get here we have successfully created the element and added it to the
  // domain
  return TCL_OK;
}



/*  *****************************************************************************
    BBAR  BRICK  U_P
    *****************************************************************************
 */

int
TclBasicBuilder_addBBarBrickUP(ClientData clientData, Tcl_Interp *interp, ArgSize argc,
                               TCL_Char ** const argv)
{
  ModelRegistry *builder = (ModelRegistry*)clientData;
  Domain* theTclDomain = builder->getDomain();


  if (builder == nullptr || clientData == nullptr) {
    opserr << "WARNING builder has been destroyed\n";
    return TCL_ERROR;
  }

  if (builder->getNDM() != 3 || builder->getNDF() != 4) {
    opserr << "WARNING -- model dimensions and/or nodal DOF not compatible "
              "with QuadUP element\n";
    return TCL_ERROR;
  }

  // check the number of arguments is correct
  int argStart = 2;

  if ((argc - argStart) < 15) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: element BBarBrickUP eleTag? N1? N2? N3? N4? N5? N6? N7? "
              "N8? matTag? bulk? rhof? perm_x? perm_y? perm_z? <b1? b2? b3?>\n";
    return TCL_ERROR;
  }

  // get the id and end nodes
  int BBarBrickUPId, Nod[8], matID;
  double bk, r, perm1, perm2, perm3;
  double b1 = 0.0;
  double b2 = 0.0;
  double b3 = 0.0;

  if (Tcl_GetInt(interp, argv[argStart], &BBarBrickUPId) != TCL_OK) {
    opserr << "WARNING invalid eleTag" << "\n";
    return TCL_ERROR;
  }

  for (int i = 0; i < 8; i++)
    if (Tcl_GetInt(interp, argv[1 + argStart + i], &(Nod[i])) != TCL_OK) {
      opserr << "WARNING invalid Node number\n";
      return TCL_ERROR;
    }

  if (Tcl_GetInt(interp, argv[9 + argStart], &matID) != TCL_OK) {
    opserr << "WARNING invalid matID\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[10 + argStart], &bk) != TCL_OK) {
    opserr << "WARNING invalid fluid bulk modulus\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[11 + argStart], &r) != TCL_OK) {
    opserr << "WARNING invalid fluid mass density\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[12 + argStart], &perm1) != TCL_OK) {
    opserr << "WARNING invalid permeability_x\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[13 + argStart], &perm2) != TCL_OK) {
    opserr << "WARNING invalid permeability_y\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[14 + argStart], &perm3) != TCL_OK) {
    opserr << "WARNING invalid permeability_z\n";
    return TCL_ERROR;
  }

  if ((argc - argStart) >= 16) {
    if (Tcl_GetDouble(interp, argv[15 + argStart], &b1) != TCL_OK) {
      opserr << "WARNING invalid b1\n";
      return TCL_ERROR;
    }
  }
  if ((argc - argStart) >= 17) {
    if (Tcl_GetDouble(interp, argv[16 + argStart], &b2) != TCL_OK) {
      opserr << "WARNING invalid b2\n";
      return TCL_ERROR;
    }
  }
  if ((argc - argStart) >= 18) {
    if (Tcl_GetDouble(interp, argv[17 + argStart], &b3) != TCL_OK) {
      opserr << "WARNING invalid b3\n";
      return TCL_ERROR;
    }
  }

  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == nullptr) {
    return TCL_ERROR;
  }

  // now create the BBarBrickUP and add it to the Domain
  BBarBrickUP *theBBarBrickUP = new BBarBrickUP(
      BBarBrickUPId, Nod[0], Nod[1], Nod[2], Nod[3], Nod[4], Nod[5], Nod[6],
      Nod[7], *theMaterial, bk, r, perm1, perm2, perm3, b1, b2, b3);


  if (theTclDomain->addElement(theBBarBrickUP) == false) {
    opserr << "WARNING could not add element to the domain\n";
    delete theBBarBrickUP;
    return TCL_ERROR;
  }

  // if get here we have successfully created the element and added it to the
  // domain
  return TCL_OK;
}
