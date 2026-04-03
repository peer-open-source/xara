//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the implementation of the commands used
// to add zero length elements to the model.
//
// Written: fmk
// Created: 01/00
//
#include <array>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include <Parsing.h>
#include <Logging.h>

#include <ZeroLength.h>
#include <ZeroLengthND.h>
#include <ZeroLengthSection.h>
#include <ZeroLengthContact2D.h>
#include <ZeroLengthRocking.h>
#include <ZeroLengthContact3D.h>
#include <ModelRegistry.h>
#include <ID.h>
#include <Node.h>
#include <Vector.h>
#include <Domain.h>
#include <UniaxialMaterial.h>
#include <NDMaterial.h>
#include <FrameSection.h>
#include <SectionForceDeformation.h>
#include <Matrix3D.h>
#include <Vector3D.h>


using namespace OpenSees;

int
ParseLinkAxes(Tcl_Interp *interp, Tcl_Size argc, TCL_Char ** const argv, 
              int ndm, const Vector3D& dx, Matrix3D& R);
            
int
TclCommand_addZeroLength(ClientData clientData, Tcl_Interp *interp, Tcl_Size argc,
                              TCL_Char ** const argv)
{
  ModelRegistry* builder = (ModelRegistry*)clientData;
  Domain* domain = builder->getDomain();

  // Syntax:
  // element ZeroLength eleTag? iNode? jNode? 
  //          -mat matID1? ... -dir dirMat1? .. 
  //          <-orient x1? x2? x3? y1? y2? y3?>

  // the spatial dimension of the problem
  const int ndm = builder->getNDM();

  //
  // first scan the command line to obtain eleID, iNode, jNode, material ID's
  // and their directions, and the orientation of ele xPrime and yPrime not
  // along the global x and y axis
  //

  int eleTag;
  std::array<int, 2> nodes;

  // a quick check on number of args
  if (argc < 9) {
    opserr << "WARNING too few arguments "
           << "\n";
    return TCL_ERROR;
  }

  // get the ele tag
  if (Tcl_GetInt(interp, argv[2], &eleTag) != TCL_OK) {
    opserr << "WARNING invalied eleTag " << argv[2]
           << "\n";
    return TCL_ERROR;
  }

  // get the two end nodes
  if (Tcl_GetInt(interp, argv[3], &nodes[0]) != TCL_OK) {
    opserr << "WARNING invalied iNode " << argv[3]
           << "\n";

    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[4], &nodes[1]) != TCL_OK) {
    opserr << "WARNING invalid jNode " << argv[4]
           << "\n";
    return TCL_ERROR;
  }

  Matrix3D T{};
  {
    Vector3D dx{};
    Domain *theDomain = builder->getDomain();
    Node *ndI = theDomain->getNode(nodes[0]);
    Node *ndJ = theDomain->getNode(nodes[1]);
    if (ndI == nullptr || ndJ == nullptr) {
      opserr << OpenSees::PromptValueError << "invalid node tag(s)\n";
      return TCL_ERROR;
    }
    const Vector &end1 = ndI->getCrds();
    const Vector &end2 = ndJ->getCrds();
    for (int i=0; i<ndm; ++i)
      dx(i) = end2(i) - end1(i);
    
    if (ParseLinkAxes(interp, argc-5, argv+5, ndm, dx, T) != TCL_OK)
      return TCL_ERROR;
  }


  // create an array of material pointers, to do this first count
  // the materials to create the array then get matID's and from ModelBuilder
  // obtain pointers to the material objects

  // read the number of materials
  // int numMat = 0;
  // int numDir = 0;
  int argi = 5;

  // if (strcmp(argv[5], "-mat") != 0) {
  //   opserr << "WARNING expecting -mat flag %s %s %s %s\n"
  //          << argv[5] << "- element ZeroLength eleTag? iNode? jNode? "
  //          << "-mat matID1? ... -dir dirMat1? .. "
  //          << "<-orient x1? x2? x3? y1? y2? y3?>\n";

  //   return TCL_ERROR;
  // }

  // while ((argi < argc) && (strcmp(argv[argi], "-dir") != 0)) {
  //   numMat++;
  //   argi++;
  // }

  // if (argi == argc) { // check we encountered the -dirn flag
  //   opserr << "WARNING no -dirn flag encountered "
  //          << "\n";
  //   return TCL_ERROR;
  // }

  // if (numMat == 0) {
  //   opserr << "WARNING no materials specified "
  //          << "\n";
  //   return TCL_ERROR;
  // }

  // create the array
  // UniaxialMaterial **theMats = new UniaxialMaterial *[numMat];
  // UniaxialMaterial **theDampMats = new UniaxialMaterial *[numMat];
  std::vector<int> mat_tags;
  std::vector<int> dof_tags;
  std::vector<int> dmp_tags;

  // fill in the material array
  // for (int i = 0; i < numMat; i++) {
  //   theDampMats[i] = 0;

  //   int matID;

  //   // read the material tag
  //   if (Tcl_GetInt(interp, argv[argi], &matID) != TCL_OK) {
  //     opserr << "WARNING invalid matID " << argv[argi]
  //            << "- element ZeroLength eleTag? iNode? jNode? "
  //            << "-mat matID1? ... -dir dirMat1? .. "
  //            << "<-orient x1? x2? x3? y1? y2? y3?>\n";
  //     delete[] theMats;
  //     return TCL_ERROR;
  //   } 
  //   else {

  //     // get a pointer to the material from the modelbuilder
  //     argi++;
  //     UniaxialMaterial *theMat = builder->getTypedObject<UniaxialMaterial>(matID);
  //     if (theMat == nullptr) {
  //       delete[] theMats;
  //       return TCL_ERROR;
  //     }
  //     // add the material to the array
  //     theMats[i] = theMat;
  //   }
  // }

  // // now read the dirn ID's for the materials added
  // argi = 6 + numMat;
  // if (strcmp(argv[argi], "-dir") != 0) {
  //   opserr << "WARNING expecting -dirn flag " << argv[argi]
  //          << "\n";
  //   delete[] theMats;
  //   return TCL_ERROR;
  // }
  // if ((argi + numMat) > argc) {
  //   opserr << "WARNING not enough directions provided for element " << eleTag
  //          << "\n";

  //   delete[] theMats;
  //   return TCL_ERROR;
  // }

  // create an ID to hold the directions
  // ID theDirns(numMat);
  // argi++;
  // int dirnID;

  // const int ndf = builder->getNDF();
  // // read the dirn identifiers
  // for (int j = 0; j < numMat; j++) {
  //   if ((Tcl_GetInt(interp, argv[argi], &dirnID) != TCL_OK) || (dirnID < 1) || (dirnID > ndf)) {
  //     opserr << OpenSees::PromptValueError 
  //            << "invalid directiion " << argv[argi]
  //            << "\n";

  //     delete[] theMats;
  //     return TCL_ERROR;
  //   }
  //   else {
  //     theDirns[j] = dirnID - 1; // the minus g3 to C++
  //     argi++;
  //   }
  // }

  // finally check the command line to see if user specified orientation
  int doRayleighDamping = 0;

  while (argi < argc) {
    if (strcmp(argv[argi], "-doRayleigh") == 0) {
      doRayleighDamping = 1;
      argi++;
      if (argi < argc)
        if ((Tcl_GetInt(interp, argv[argi], &doRayleighDamping) == TCL_OK))
          argi++;
    }
    else if ((strcmp(argv[argi], "-dir") == 0) || 
             (strcmp(argv[argi], "-dof") == 0)) {
      argi++;
      while (argi < argc && argv[argi][0] != '-') {
        int dofID;
        if (Tcl_GetInt(interp, argv[argi], &dofID) != TCL_OK) {
          break;
        }
        dof_tags.push_back(dofID);
        argi++;
      }
    }
    else if ((strcmp(argv[argi], "-mat") == 0) || 
             (strcmp(argv[argi], "-mats") == 0)) {
      argi++;
      while (argi < argc && argv[argi][0] != '-') {
        int matID;
        if (Tcl_GetInt(interp, argv[argi], &matID) != TCL_OK) {
          break;
        }
        mat_tags.push_back(matID);
        argi++;
      }
    }
    else if (strcmp(argv[argi], "-dampMats") == 0) {
      doRayleighDamping = 2;
      argi++;
      while (argi < argc && argv[argi][0] != '-') {
        int matID;
        if (Tcl_GetInt(interp, argv[argi], &matID) != TCL_OK) {
          break;
        }
        dmp_tags.push_back(matID);
        argi++;
      }

    } else
      argi++;
  }

  if (dof_tags.size() == 0) {
    opserr << OpenSees::PromptValueError 
           << "no directions specified for element " << eleTag
           << "\n";
    return TCL_ERROR;
  }
  if (mat_tags.size() == 0) {
    opserr << OpenSees::PromptValueError 
           << "no materials specified for element " << eleTag
           << "\n";
    return TCL_ERROR;
  }
  if (dof_tags.size() != mat_tags.size()) {
    opserr << OpenSees::PromptValueError 
           << "number of materials, directions and damping materials do not match\n";
    return TCL_ERROR;
  }
  if (doRayleighDamping == 2 && dmp_tags.size() != mat_tags.size()) {
    opserr << OpenSees::PromptValueError 
           << "number of materials, directions and damping materials do not match\n";
    return TCL_ERROR;
  }

  //
  // Get the material pointers from the model builder
  //
  int numMat = mat_tags.size();
  UniaxialMaterial **theMats = new UniaxialMaterial *[numMat];
  UniaxialMaterial **theDampMats = new UniaxialMaterial *[numMat];
  ID theDirns(numMat);
  for (int i=0; i<numMat; ++i) {
    UniaxialMaterial *theMat = builder->getTypedObject<UniaxialMaterial>(mat_tags[i]);
    if (theMat == nullptr) {
      delete[] theMats;
      delete[] theDampMats;
      return TCL_ERROR;
    }
    theMats[i] = theMat;

    if (doRayleighDamping == 2) {
      UniaxialMaterial *theDampMat = builder->getTypedObject<UniaxialMaterial>(dmp_tags[i]);
      if (theDampMat == nullptr) {
        delete[] theMats;
        delete[] theDampMats;
        return TCL_ERROR;
      }
      theDampMats[i] = theDampMat;
    }

    theDirns(i) = dof_tags[i] - 1;
  }

  //
  // now we create the element and add it to the domain
  //

  Element *theEle = nullptr;

  if (doRayleighDamping != 2)
    theEle = new ZeroLength(eleTag, ndm, nodes[0], nodes[1],
                            T,
                            numMat, theMats,
                            theDirns, doRayleighDamping);
  else
    theEle = new ZeroLength(eleTag, ndm, nodes[0], nodes[1], 
                            T,
                            numMat, theMats,
                            theDampMats, theDirns, doRayleighDamping);


  if (domain->addElement(theEle) == false) {
    delete[] theMats;
    return TCL_ERROR;
  }

  // return the memory we stole and return OK
  delete[] theMats;
  delete[] theDampMats;
  return TCL_OK;
}

int
TclCommand_addZeroLengthSection(ClientData clientData, Tcl_Interp *interp,
                                     int argc, TCL_Char ** const argv)
{

  ModelRegistry* builder = (ModelRegistry*)clientData;
  Domain* domain = builder->getDomain();
  int ndm = builder->getNDM(); // the spatial dimension of the problem

  //
  // first scan the command line to obtain eleID, iNode, jNode, material ID's
  // and their directions, and the orientation of ele xPrime and yPrime not
  // along the global x and y axis
  //

  int eleTag, iNode, jNode;

  // a quick check on number of args
  if (argc < 6) {
    opserr << "WARNING too few arguments "
           << "want - element zeroLengthSection eleTag? iNode? jNode? "
           << "secTag? "
           << "<-orient x1? x2? x3? y1? y2? y3?>\n";
    return TCL_ERROR;
  }

  // get the ele tag
  if (Tcl_GetInt(interp, argv[2], &eleTag) != TCL_OK) {
    opserr << "WARNING invalied eleTag " << argv[2]
           << "\n";
    return TCL_ERROR;
  }

  // get the two end nodes
  if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK) {
    opserr << "WARNING invalied iNode " << argv[3]
           << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK) {
    opserr << "WARNING invalid jNode " << argv[4]
           << "\n";
    return TCL_ERROR;
  }

  int secTag;
  if (Tcl_GetInt(interp, argv[5], &secTag) != TCL_OK) {
    opserr << "WARNING invalid secTag " << argv[5]
           << "\n";
    return TCL_ERROR;
  }

  // create the vectors for the element orientation
  Vector x(3);
  x(0) = 1.0;
  x(1) = 0.0;
  x(2) = 0.0;
  Vector y(3);
  y(0) = 0.0;
  y(1) = 1.0;
  y(2) = 0.0;

  int argi = 6;
  int doRayleighDamping = 1;

  // finally check the command line to see if user specified orientation
  while (argi < argc) {
    if (strcmp(argv[argi], "-orient") == 0) {
      if (argc < (argi + 7)) {
        opserr << "WARNING not enough parameters after -orient flag for ele "
               << eleTag
               << "- element zeroLengthSection eleTag? iNode? jNode? secTag? "
               << "<-orient x1? x2? x3? y1? y2? y3?>\n";
        return TCL_ERROR;
      } else {
        argi++;
        double value;
        // read the x values
        for (int i = 0; i < 3; i++) {
          if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK) {
            opserr
                << "WARNING invalid -orient value for ele  " << eleTag
                << argv[argi]
                << "- element zeroLengthSection eleTag? iNode? jNode secTag? "
                << "<-orient x1? x2? x3? y1? y2? y3?>\n";
            return TCL_ERROR;
          } else {
            argi++;
            x(i) = value;
          }
        }
        // read the y values
        for (int j = 0; j < 3; j++) {
          if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK) {
            opserr
                << "WARNING invalid -orient value for ele  " << eleTag
                << argv[argi]
                << "- element zeroLengthSection eleTag? iNode? jNode? secTag? "
                << "<-orient x1? x2? x3? y1? y2? y3?>\n";
            return TCL_ERROR;
          } else {
            argi++;
            y(j) = value;
          }
        }
      }
    } else if (strcmp(argv[argi], "-doRayleigh") == 0) {
      doRayleighDamping = 1;
      argi++;
      if (argi < argc)
        if ((Tcl_GetInt(interp, argv[argi], &doRayleighDamping) == TCL_OK))
          argi++;
    } else
      argi++;
  }

  //
  // now we create the element and add it to the domain
  //
  SectionForceDeformation *theSection = builder->getTypedObject<FrameSection>(secTag);
  if (theSection == nullptr)
    return TCL_ERROR;

  Element *theEle = new ZeroLengthSection(eleTag, ndm, iNode, jNode, x, y,
                                          *theSection, doRayleighDamping);

  if (domain->addElement(theEle) == false)
    return TCL_ERROR;

  return TCL_OK;
}

// add by Gang Wang for Contact Element
// Tcl parse of zeroLengthContact2D command
//
// Command:
//
// element zeroLengthContact2D $tag $secondaryNd $primaryNd $Kn $Kt $fs
// $ContactDirction
//
//
int
TclCommand_addZeroLengthContact2D(ClientData clientData, Tcl_Interp *interp,
                                       int argc, TCL_Char ** const argv)
{
  ModelRegistry* builder = (ModelRegistry*)clientData;
  Domain* domain = builder->getDomain();

  //
  // first scan the command line to obtain eleID, SecondaryNode, PrimaryNode,
  //
  int eleTag, iNode, jNode;

  // opserr << argc;

  // a quick check on number of args
  if (argc < 11) {
    opserr << "ZeroLengthContact2D::WARNING too few arguments "
           << "want - element ZeroLengthContact2D eleTag? iNode? jNode? Kn? "
              "Kt? fs? -normal Nx? Ny?";
    return TCL_ERROR;
  }

  // get the ele tag
  if (Tcl_GetInt(interp, argv[2], &eleTag) != TCL_OK) {
    opserr << "ZeroLengthContact2D::WARNING invalied eleTag " << argv[2]
           << "\n";
    return TCL_ERROR;
  }

  // get the two end nodes
  if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK) {
    opserr << "ZeroLengthContact2D::WARNING invalied iNode " << argv[3] << "\n";

    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK) {
    opserr << "ZeroLengthContact2D::WARNING invalid jNode " << argv[4] << "\n";
    return TCL_ERROR;
  }

  double Kn, Kt, fs;

  // read the material properties
  if (Tcl_GetDouble(interp, argv[5], &Kn) != TCL_OK) {
    opserr << "ZeroLengthContact2D::WARNING invalid Kn " << argv[5] << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[6], &Kt) != TCL_OK) {
    opserr << "ZeroLengthContact2D::WARNING invalid Kt " << argv[6] << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[7], &fs) != TCL_OK) {
    opserr << "ZeroLengthContact2D::WARNING invalid fs " << argv[7] << "\n";
    return TCL_ERROR;
  }

  // changed to specify any contact normal direction
  if (strcmp(argv[8], "-normal") != 0) {
    opserr << "ZeroLengthContact2D:: expecting "
           << "- element ZeroLengthContact2D eleTag? iNode? jNode? Kn? Kt? fs? "
              "-normal Nx? Ny? \n";
    return TCL_ERROR;
  }

  Vector NormalDir(2);

  int argi = 9;

  double value;
  // read the NormalDir values
  for (int i = 0; i < 2; i++) {
    if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK) {
      opserr << "ZeroLengthContact2D:: invalid -normal value for ele " << eleTag
             << "- element ZeroLengthContact2D eleTag? iNode? jNode? Kn? Kt? "
                "fs? -normal Nx? Ny? \n";

      return TCL_ERROR;
    } else {
      argi++;
      NormalDir(i) = value;
    }
  }

  //
  // now we create the element and add it to the domain
  //

  Element *theEle;

  theEle = new ZeroLengthContact2D(eleTag, iNode, jNode, Kn, Kt, fs, NormalDir);

  if (theEle == 0) {
    return TCL_ERROR;
  }

  if (domain->addElement(theEle) == false) {
    return TCL_ERROR;
  }

  // return the memory we stole and return OK
  return TCL_OK;
}

// add by Gang Wang for Contact Element
// Tcl parse of zeroLengthContact3D command

// Command:
//
// element zeroLengthContact3D $tag $secondaryNd $primaryNd $Kn $Kt $fs $c
// $ContactDir
//
//

int
TclCommand_addZeroLengthContact3D(ClientData clientData, Tcl_Interp *interp,
                                       int argc, TCL_Char ** const argv)
{
  ModelRegistry* builder = (ModelRegistry*)clientData;
  Domain* domain = builder->getDomain();

  //
  // first scan the command line to obtain eleID, SecondaryNode, PrimaryNode,

  int eleTag, iNode, jNode;

  // a quick check on number of args
  if (argc < 10) {
    opserr << "ZeroLengthContact3D::WARNING too few arguments "
           << "want - element ZeroLengthContact3D eleTag? iNode? jNode? Kn? "
              "Kt? fs? c? dir?";
    return TCL_ERROR;
  }

  // get the ele tag
  if (Tcl_GetInt(interp, argv[2], &eleTag) != TCL_OK) {
    opserr << "ZeroLengthContact3D::WARNING invalied eleTag " << argv[2]
           << "\n";
    return TCL_ERROR;
  }

  // get the two end nodes
  if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK) {
    opserr << "ZeroLengthContact3D::WARNING invalied iNode " << argv[3] << "\n";

    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK) {
    opserr << "ZeroLengthContact3D::WARNING invalid jNode " << argv[4] << "\n";
    return TCL_ERROR;
  }

  double Kn, Kt, fs, c;

  // read the material properties
  if (Tcl_GetDouble(interp, argv[5], &Kn) != TCL_OK) {
    opserr << "ZeroLengthContact3D::WARNING invalid Kn " << argv[5] << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[6], &Kt) != TCL_OK) {
    opserr << "ZeroLengthContact3D::WARNING invalid Kt " << argv[6] << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[7], &fs) != TCL_OK) {
    opserr << "ZeroLengthContact3D::WARNING invalid fs " << argv[7] << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[8], &c) != TCL_OK) {
    opserr << "ZeroLengthContact3D::WARNING invalid c " << argv[8] << "\n";
    return TCL_ERROR;
  }

  int dir;

  if (Tcl_GetInt(interp, argv[9], &dir) != TCL_OK) {
    opserr << "ZeroLengthContact3D::WARNING invalid direction " << argv[9]
           << "\n";
    return TCL_ERROR;
  }

  //
  // now we create the element and add it to the domain
  //

  Element *theEle;

  double originX, originY;
  originX = 0;
  originY = 0;

  if (dir == 0) {
    if (argc == 12) {
      if (Tcl_GetDouble(interp, argv[10], &originX) != TCL_OK) {
        opserr << "ZeroLengthContact3D::WARNING invalid originX " << argv[9]
               << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[11], &originY) != TCL_OK) {
        opserr << "ZeroLengthContact3D::WARNING invalid originY " << argv[10]
               << "\n";
        return TCL_ERROR;
      }
    }
  }

  theEle = new ZeroLengthContact3D(eleTag, iNode, jNode, dir, Kn, Kt, fs, c,
                                   originX, originY);

  if (theEle == 0) {
    return TCL_ERROR;
  }

  if (domain->addElement(theEle) == false) {
    return TCL_ERROR;
  }

  // return the memory we stole and return OK
  return TCL_OK;
}

int
TclCommand_addZeroLengthND(ClientData clientData, Tcl_Interp *interp, int argc,
                                TCL_Char ** const argv)
{
  ModelRegistry* builder = (ModelRegistry*)clientData;
  Domain* domain = builder->getDomain();
  int ndm = builder->getNDM(); // the spatial dimension of the problem

  //
  // first scan the command line to obtain eleID, iNode, jNode, material ID's
  // and their directions, and the orientation of ele xPrime and yPrime not
  // along the global x and y axis
  //

  int eleTag, iNode, jNode;

  // a quick check on number of args
  if (argc < 6) {
    opserr << "WARNING too few arguments %s %s %s\n"
           << "want - element zeroLengthND eleTag? iNode? jNode? "
           << "NDTag? <1DTag?>"
           << "<-orient x1? x2? x3? y1? y2? y3?>\n";
    return TCL_ERROR;
  }

  // get the ele tag
  if (Tcl_GetInt(interp, argv[2], &eleTag) != TCL_OK) {
    opserr << "WARNING invalied eleTag " << argv[2]
           << " - element zeroLengthND eleTag? iNode? jNode? NDTag? <1DTag?> "
              "<-orient x1? x2? x3? y1? y2? y3?>\n";

    return TCL_ERROR;
  }

  // get the two end nodes
  if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK) {
    opserr << "WARNING invalied iNode " << argv[3]
           << "- element zeroLengthND eleTag? iNode? jNode? "
           << "NDTag? <1DTag?>"
           << "<-orient x1? x2? x3? y1? y2? y3?>\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK) {
    opserr << "WARNING invalid jNode " << argv[4]
           << "- element zeroLengthND eleTag? iNode? jNode? "
           << "NDTag? <1DTag?> <-orient x1? x2? x3? y1? y2? y3?>\n";
    return TCL_ERROR;
  }

  int NDTag;

  if (Tcl_GetInt(interp, argv[5], &NDTag) != TCL_OK) {
    opserr << "WARNING invalid NDTag "
           << argv[5]
           << "\n";

    return TCL_ERROR;
  }

  UniaxialMaterial *the1DMat = 0;

  int argi = 6;

  if (argc > 6 && strcmp(argv[6], "-orient") != 0) {

    int uniTag;

    if (Tcl_GetInt(interp, argv[6], &uniTag) != TCL_OK) {
      opserr << "WARNING invalid NDTag " << argv[5]
             << "- element zeroLengthND eleTag? iNode? jNode? "
             << "NDTag? <1DTag?> <-orient x1? x2? x3? y1? y2? y3?>\n";
      return TCL_ERROR;
    }

    the1DMat = builder->getTypedObject<UniaxialMaterial>(uniTag);

    if (the1DMat == nullptr)
      opserr << "WARNING UniaxialMaterial " << uniTag
             << " not found in model, proceeding without\n";

    argi = 7;
  }

  // create the vectors for the element orientation
  Vector x(3);
  x(0) = 1.0;
  x(1) = 0.0;
  x(2) = 0.0;
  Vector y(3);
  y(0) = 0.0;
  y(1) = 1.0;
  y(2) = 0.0;

  // finally check the command line to see if user specified orientation
  if (argi < argc) {
    if (strcmp(argv[argi], "-orient") == 0) {
      if (argc < (argi + 7)) {
        opserr << "WARNING not enough parameters after -orient flag for ele "
               << eleTag << "- element zeroLengthND eleTag? iNode? jNode? "
               << "NDTag? <1DTag?> <-orient x1? x2? x3? y1? y2? y3?>\n";
        return TCL_ERROR;
      } else {
        argi++;
        double value;
        // read the x values
        for (int i = 0; i < 3; i++) {
          if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK) {
            opserr << "WARNING invalid -orient value for ele  " << eleTag
                   << argv[argi]
                   << "- element zeroLengthND eleTag? iNode? jNode? "
                   << "NDTag? <1DTag?> <-orient x1? x2? x3? y1? y2? y3?>\n";

            return TCL_ERROR;
          } else {
            argi++;
            x(i) = value;
          }
        }
        // read the y values
        for (int j = 0; j < 3; j++) {
          if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK) {
            opserr << "WARNING invalid -orient value for ele  " << eleTag << " "
                   << argv[argi]
                   << "- element zeroLengthND eleTag? iNode? jNode? "
                   << "NDTag? <1DTag?> <-orient x1? x2? x3? y1? y2? y3?>\n";
            return TCL_ERROR;
          } else {
            argi++;
            y(j) = value;
          }
        }
      }
    }
  }

  //
  // now we create the element and add it to the domain
  //

  NDMaterial *theNDMat = builder->getTypedObject<NDMaterial>(NDTag); // OPS_getNDMaterial(NDTag);
  if (theNDMat == nullptr)
    return TCL_ERROR;


  Element *theEle = nullptr;

  if (the1DMat == nullptr)
    theEle = new ZeroLengthND(eleTag, ndm, iNode, jNode, x, y, *theNDMat);
  else
    theEle =
        new ZeroLengthND(eleTag, ndm, iNode, jNode, x, y, *theNDMat, *the1DMat);

  if (theEle == 0)
    return TCL_ERROR;

  if (domain->addElement(theEle) == false)
    return TCL_ERROR;

  return TCL_OK;
}

int
TclCommand_addZeroLengthRocking(ClientData clientData, Tcl_Interp *interp,
                                     int argc, TCL_Char ** const argv)
{
  ModelRegistry* builder = (ModelRegistry*)clientData;
  Domain* domain = builder->getDomain();
  int ndm = builder->getNDM(); // the spatial dimension of the problem

  //
  // first scan the command line to obtain eleID, iNode, jNode, and the
  // orientation of ele xPrime and yPrime not along the global x and y axis
  //

  int eleTag, iNode, jNode;

  // a quick check on number of args
  if (argc < 9) {
    opserr << "WARNING too few arguments "
           << "want - element ZeroLengthRocking eleTag? iNode? jNode? "
           << "kr? radius? theta0? kappa? <-orient x1? x2? x3? y1? y2? y3?>\n";

    return TCL_ERROR;
  }

  // get the ele tag
  if (Tcl_GetInt(interp, argv[2], &eleTag) != TCL_OK) {
    opserr << "WARNING invalied eleTag " << argv[2]
           << "- element ZeroLengthRocking eleTag? iNode? jNode? "
           << "kr? radius? theta0? kappa? <-orient x1? x2? x3? y1? y2? y3?>\n";
    return TCL_ERROR;
  }

  // get the two end nodes
  if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK) {
    opserr << "WARNING invalied iNode " << argv[3]
           << "- element ZeroLengthRocking eleTag? iNode? jNode? "
           << "kr? radius? theta0? kappa? <-orient x1? x2? x3? y1? y2? y3?>\n";

    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK) {
    opserr << "WARNING invalid jNode " << argv[4]
           << "- element ZeroLengthRocking eleTag? iNode? jNode? "
           << "kr? radius? theta0? kappa? <-orient x1? x2? x3? y1? y2? y3?>\n";
    return TCL_ERROR;
  }

  // look for rocking required inputs
  double kr = 0;
  double R = 0;
  double theta = 0;
  double kap = 1.0e12;

  // rotational stiffness
  if (Tcl_GetDouble(interp, argv[5], &kr) != TCL_OK) {
    opserr << "WARNING invalid kr " << argv[5]
           << "- element ZeroLengthRocking eleTag? iNode? jNode? "
           << "kr? radius? theta0? kappa? <-orient x1? x2? x3? y1? y2? y3?>\n";
    return TCL_ERROR;
  }

  // rocking radius
  if (Tcl_GetDouble(interp, argv[6], &R) != TCL_OK) {
    opserr << "WARNING invalid radius " << argv[6]
           << "- element ZeroLengthRocking eleTag? iNode? jNode? "
           << "kr? radius? theta0? kappa? <-orient x1? x2? x3? y1? y2? y3?>\n";
    return TCL_ERROR;
  }

  // theta0
  if (Tcl_GetDouble(interp, argv[7], &theta) != TCL_OK) {
    opserr << "WARNING invalid theta0 " << argv[7]
           << "- element ZeroLengthRocking eleTag? iNode? jNode? "
           << "kr? radius? theta0? kappa? <-orient x1? x2? x3? y1? y2? y3?>\n";
    return TCL_ERROR;
  }

  // kappa
  if (Tcl_GetDouble(interp, argv[8], &kap) != TCL_OK) {
    opserr << "WARNING invalid kappa " << argv[8]
           << "- element ZeroLengthRocking eleTag? iNode? jNode? "
           << "kr? radius? theta0? kappa? <-orient x1? x2? x3? y1? y2? y3?>\n";
    return TCL_ERROR;
  }

  // create the vectors for the element orientation
  Vector x(3);
  x(0) = 1.0;
  x(1) = 0.0;
  x(2) = 0.0;
  Vector y(3);
  y(0) = 0.0;
  y(1) = 1.0;
  y(2) = 0.0;
  double xi = 1.0e-8;
  double dTol = 1.0e-7;
  double vTol = 1.0e-7;
  int argi = 9;

  while (argi < argc) {
    if (strcmp(argv[argi], "-orient") == 0) {
      if (argc < (argi + 7)) {
        opserr
            << "WARNING not enough parameters after -orient flag for ele "
            << eleTag << "- element ZeroLengthRocking eleTag? iNode? jNode? "
            << "kr? radius? theta0? kappa? <-orient x1? x2? x3? y1? y2? y3?>\n";
        return TCL_ERROR;

      } else {
        argi++;
        double value;
        // read the x values
        for (int i = 0; i < 3; i++) {
          if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK) {
            opserr << "WARNING invalid -orient value for ele  " << eleTag
                   << argv[i] << "- element ZeroLength eleTag? iNode? jNode? "
                   << "kr? radius? theta0? kappa? <-orient x1? x2? x3? y1? y2? "
                      "y3?>\n";
            return TCL_ERROR;
          } else {
            argi++;
            x(i) = value;
          }
        }
        // read the y values
        for (int j = 0; j < 3; j++) {
          if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK) {
            opserr << "WARNING invalid -orient value for ele  " << eleTag
                   << argv[argi]
                   << "- element ZeroLength eleTag? iNode? jNode? "
                   << "kr? radius? theta0? kappa? <-orient x1? x2? x3? y1? y2? "
                      "y3?>\n";
            return TCL_ERROR;
          } else {
            argi++;
            y(j) = value;
          }
        }
      }

    } else if (strcmp(argv[argi], "-xi") == 0) {
      if (argc < (argi + 2)) {
        opserr << "WARNING not enough parameters after -xi flag for ele "
               << eleTag << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      } else {
        argi++;
        if (Tcl_GetDouble(interp, argv[argi], &xi) != TCL_OK) {
          opserr << "WARNING invalid -xi value for ele  " << eleTag << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        } else
          argi++;
      }

    } else if (strcmp(argv[argi], "-dTol") == 0) {
      if (argc < (argi + 2)) {
        opserr << "WARNING not enough parameters after -dTol flag for ele "
               << eleTag << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      } else {
        argi++;
        if (Tcl_GetDouble(interp, argv[argi], &dTol) != TCL_OK) {
          opserr << "WARNING invalid -dTol value for ele  " << eleTag << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        } else
          argi++;
      }

    } else if (strcmp(argv[argi], "-vTol") == 0) {
      if (argc < (argi + 2)) {
        opserr << "WARNING not enough parameters after -vTol flag for ele "
               << eleTag << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
      } else {
        argi++;
        if (Tcl_GetDouble(interp, argv[argi], &vTol) != TCL_OK) {
          opserr << "WARNING invalid -vTol value for ele  " << eleTag << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        } else
          argi++;
      }

    } else
      argi++;
  }

  //
  // now we create the element and add it to the domain
  //
  Element *theEle;
  theEle = new ZeroLengthRocking(eleTag, ndm, iNode, jNode, x, y, kr, R, theta,
                                 kap, xi, dTol, vTol);
  if (theEle == 0)
    return TCL_ERROR;

  if (domain->addElement(theEle) == false)
    return TCL_ERROR;

  return TCL_OK;
}
