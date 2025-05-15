#include <tcl.h>
#include <assert.h>
#include <Domain.h>
#include <BasicModelBuilder.h>
#include <elementAPI.h>
#include <element/community/UWelements/SSPquadUP.h>
#include <element/community/UWelements/SSPbrick.h>
#include <NDMaterial.h>

#ifdef _MSC_VER 
#  include <string.h>
#  define strcasecmp _stricmp
#else
#  include <strings.h>
#endif
#define strcmp strcasecmp

static Element *TclDispatch_SSPbrick(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv);
static Element *TclDispatch_SSPbrickUP(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv);
static int      TclCommand_addSSPquad(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv);
static Element *TclDispatch_SSPquadUP(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv);

int
TclCommand_SSP_Element(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  Element* theEle = nullptr;
  assert(clientData != nullptr);

  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  Domain* domain = builder->getDomain();

  if (strcasecmp(argv[1], "SSPquad")==0) {
    return TclCommand_addSSPquad(clientData, interp, argc, argv);
  }
  else if (strcasecmp(argv[1], "SSPquadUP")==0) {
    theEle = TclDispatch_SSPquadUP(clientData, interp, argc, argv);
  }
  else if (strcasecmp(argv[1], "SSPbrick")==0) {
    theEle = TclDispatch_SSPbrick(clientData, interp, argc, argv);
  }
  else if (strcasecmp(argv[1], "SSPbrickUP")==0) {
    theEle = TclDispatch_SSPbrickUP(clientData, interp, argc, argv);
  }

  if (theEle && domain->addElement(theEle))
    return TCL_OK;

  return TCL_ERROR;
}

static Element*
TclDispatch_SSPbrick(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  static int num_SSPbrick;
  if (num_SSPbrick == 0) {
    num_SSPbrick++;
    opslog << "SSPbrick element - Written: C.McGann, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to an element that will be returned
  Element *theElement = 0;


  if (argc < 10) {
    opserr
        << "Invalid #args, want: element SSPbrick eleTag? iNode? jNode? kNode? "
           "lNode? mNode? nNode? pNode? qNode? matTag? <b1? b2? b3?>\n";
    return 0;
  }

  int iData[10];
  double dData[3];
  dData[0] = 0.0;
  dData[1] = 0.0;
  dData[2] = 0.0;

  int numData = 10;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element SSPbrick " << iData[0]
           << endln;
    return 0;
  }

  int matID = iData[9];
  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == nullptr) {
    return nullptr;
  }

  if (argc == 13) {
    numData = 3;
    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "WARNING invalid optional data: element SSPbrick " << iData[0]
             << endln;
      return 0;
    }
  }

  // parsing was successful, allocate the element
  theElement = new SSPbrick(iData[0], iData[1], iData[2], iData[3], iData[4],
                            iData[5], iData[6], iData[7], iData[8],
                            *theMaterial, dData[0], dData[1], dData[2]);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type SSPbrick\n";
    return 0;
  }

  return theElement;
}



#include <element/community/UWelements/SSPbrickUP.h>
static Element*
TclDispatch_SSPbrickUP(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  static int num_SSPbrickUP;
  if (num_SSPbrickUP == 0) {
    num_SSPbrickUP++;
    opslog << "SSPbrickUP element - Written: C.McGann, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to an element that will be returned
  Element *theElement = 0;

  if (argc < 17) {
    opserr << "Invalid #args, want: element SSPbrickUP eleTag? iNode? jNode? "
              "kNode? lNode? mNode? nNode? pNode? qNode? matTag? fBulk? fDen? "
              "k1? k2? k3? e? alpha? <b1? b2? b3?>\n";
    return 0;
  }

  int iData[10];
  double dData[10];
  dData[7] = 0.0;
  dData[8] = 0.0;
  dData[9] = 0.0;

  int numData = 10;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element SSPbrickUP " << iData[0]
           << endln;
    return 0;
  }

  int matID = iData[9];
  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == nullptr) {
    return nullptr;
  }

  numData = 7;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid double data: element SSPbrickUP " << iData[0]
           << endln;
    return 0;
  }

  if (argc == 20) {
    numData = 3;
    if (OPS_GetDoubleInput(&numData, &dData[7]) != 0) {
      opserr << "WARNING invalid optional data: element SSPbrickUP " << iData[0]
             << endln;
      return 0;
    }
  }

  // parsing was successful, allocate the element
  theElement = new SSPbrickUP(
      iData[0], iData[1], iData[2], iData[3], iData[4], iData[5], iData[6],
      iData[7], iData[8], *theMaterial, dData[0], dData[1], dData[2], dData[3],
      dData[4], dData[5], dData[6], dData[7], dData[8], dData[9]);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type SSPbrickUP\n";
    return 0;
  }

  return theElement;
}




#include <element/community/UWelements/SSPquad.h>
static int
TclCommand_addSSPquad(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  if (builder->getNDM() != 2 || builder->getNDF() != 2) {
    opserr << "WARNING -- model dimensions and/or nodal DOF not compatible "
              "with quad element\n";
    return TCL_ERROR;
  }

  // check the number of arguments is correct
//int argStart = 2;

//if ((argc - argStart) < 8) {
//  opserr << "WARNING insufficient arguments\n";
//  opserr << "Want: element FourNodeQuad eleTag? iNode? jNode? kNode? lNode? "
//            "thk? type? matTag? <pressure? rho? b1? b2?>\n";
//  return TCL_ERROR;
//}

  // get the id and end nodes
  int tag, iNode, jNode, kNode, lNode, matID;
  double thickness = 1.0;
  double b1 = 0.0;
  double b2 = 0.0;

  static int num_SSPquad;
  if (num_SSPquad == 0) {
    num_SSPquad++;
    opslog << "SSPquad element - Written: C.McGann, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  if (argc < 10) {
    opserr << "Invalid #args, want: element SSPquad eleTag? iNode? jNode? "
              "kNode? lNode? matTag? type? thickness? <b1? b2?>?\n";
    return TCL_ERROR;
  }

#if 1
  int argi = 1;
  if (Tcl_GetInt(interp, argv[++argi], &tag) != TCL_OK) {   // 2
    opserr << "WARNING invalid SSPquad eleTag" << endln;
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[++argi], &iNode) != TCL_OK) { // 3
    opserr << "WARNING invalid iNode\n";
    opserr << "SSPquad element: " << tag << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[++argi], &jNode) != TCL_OK) { // 4
    opserr << "WARNING invalid jNode\n";
    opserr << "SSPquad element: " << tag << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[++argi], &kNode) != TCL_OK) {
    opserr << "WARNING invalid kNode\n";
    opserr << "SSPquad element: " << tag << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[++argi], &lNode) != TCL_OK) {
    opserr << "WARNING invalid lNode\n";
    opserr << "SSPquad element: " << tag << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[++argi], &thickness) != TCL_OK) { // 7
    opserr << "WARNING invalid thickness\n";
    opserr << "SSPquad element: " << tag << endln;
    return TCL_ERROR;
  }

  TCL_Char *type = argv[++argi];

  if (Tcl_GetInt(interp, argv[++argi], &matID) != TCL_OK) {
    opserr << "WARNING invalid matID\n";
    opserr << "SSPquad element: " << tag << endln;
    return TCL_ERROR;
  }

  if (argi < argc-1) {
    if (Tcl_GetDouble(interp, argv[++argi], &b1) != TCL_OK) {
      opserr << "WARNING invalid b1\n";
      opserr << "SSPquad element: " << tag << endln;
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[++argi], &b2) != TCL_OK) {
      opserr << "WARNING invalid b2\n";
      opserr << "SSPquad element: " << tag << endln;
      return TCL_ERROR;
    }
  }

  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == nullptr)
    return TCL_ERROR;


  Element *theElem =
      new SSPquad(tag, iNode, jNode, kNode, lNode, *theMaterial,
                       type, thickness, b1, b2);

  if (builder->getDomain()->addElement(theElem) == false) {
    opserr << "WARNING could not add element to the domain\n";
    opserr << "FourNodeQuad element: " << tag << endln;
    delete theElem;
    return TCL_ERROR;
  }

  return TCL_OK;
#else
  int iData[6];
  const char *theType;
  double dData[3] = {1.0, 0.0, 0.0};

  int numData = 6;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element SSPquad " << iData[0]
           << endln;
    return 0;
  }

  theType = OPS_GetString();

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid thickness data: element SSPquad " << iData[0]
           << endln;
    return 0;
  }

  int matID = iData[5];
  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);

  if (theMaterial == 0) {
    opserr << "WARNING element SSPquad " << iData[0] << endln;
    opserr << " Material: " << matID << "not found\n";
    return 0;
  }

  if (argc == 10) {
    numData = 2;
    if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
      opserr << "WARNING invalid optional data: element SSPquad " << iData[0]
             << endln;
      return 0;
    }
  }

  // parsing was successful, allocate the element
  theElement = new SSPquad(iData[0], iData[1], iData[2], iData[3], iData[4],
                           *theMaterial, theType, dData[0], dData[1], dData[2]);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type SSPquad\n";
    return 0;
  }
  return theElement;
#endif
}

static Element*
TclDispatch_SSPquadUP(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  static int num_SSPquadUP;
  if (num_SSPquadUP == 0) {
    num_SSPquadUP++;
    opslog << "SSPquadUP element - Written: C.McGann, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to an element that will be returned
  Element *theElement = nullptr;

  // LM change
  if (argc < 13) {
    opserr << "Invalid #args, want: element SSPquadUP eleTag? iNode? jNode? "
              "kNode? lNode? matTag? t? fBulk? fDen? k1? k2? e? alpha? <b1? "
              "b2?> <Pup? Plow? Pleft? Pright?>?\n";
    return 0;
  }

  int iData[6];
  double dData[13];
  dData[7] = 0.0;
  dData[8] = 0.0;
  dData[9] = 0.0;
  dData[10] = 0.0;
  dData[11] = 0.0;
  dData[12] = 0.0;
  // LM change

  int numData = 6;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element SSPquadUP " << iData[0]
           << endln;
    return 0;
  }

  numData = 7;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid double data: element SSPquadUP " << iData[0]
           << endln;
    return 0;
  }

  int matID = iData[5];
  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == nullptr)
    return nullptr;


  // LM change
  if (argc == 15) {
    numData = 2;
    if (OPS_GetDoubleInput(&numData, &dData[7]) != 0) {
      opserr << "WARNING invalid optional data: element SSPquadUP " << iData[0]
             << endln;
      return 0;
    }
  } else if (argc == 19) {
    numData = 6;
    if (OPS_GetDoubleInput(&numData, &dData[7]) != 0) {
      opserr << "WARNING invalid optional data: element SSPquadUP " << iData[0]
             << endln;
      return 0;
    }
  }

  // parsing was successful, allocate the element
  theElement = new SSPquadUP(
      iData[0], iData[1], iData[2], iData[3], iData[4], *theMaterial, dData[0],
      dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7],
      dData[8], dData[9], dData[10], dData[11], dData[12]);
  // LM change

  if (theElement == 0) {
    opserr << "WARNING could not create element of type SSPquadUP\n";
    return 0;
  }

  return theElement;
}


#if 0
#include <element/community/UWelements/BeamContact2D.h>

static Element*
TclDispatch_BeamContact2D(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  if (num_BeamContact2D == 0) {
    num_BeamContact2D++;
    opslog << "BeamContact2D element - Written: C.McGann, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to an element that will be returned
  Element *theElement = 0;

  if (argc < 9) {
    opserr << "Invalid #args, want: element BeamContact2D eleTag? iNode? "
              "jNode? secondaryNode? lambdaNode? matTag? width? gapTol? "
              "forceTol? <cSwitch>?\n";
    return 0;
  }

  int iData[6];
  double dData[3];
  int icSwitch = 0;

  int numData = 6;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element BeamContact2D " << iData[0]
           << endln;
    return 0;
  }

  numData = 3;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element BeamContact2D " << dData[0]
           << endln;
    return 0;
  }

  int matID = iData[5];
  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == 0) {
    opserr << "WARNING element BeamContact2D " << iData[0] << endln;
    opserr << " Material: " << matID << "not found\n";
    return 0;
  }

  int i = argc - 9;
  while (i >= 1) {
    numData = 1;
    if (OPS_GetIntInput(&numData, &icSwitch) != 0) {
      opserr << "WARNING invalid initial contact flag: element BeamContact2D "
             << iData[0] << endln;
      return 0;
    }
    i -= 1;
  }

  // Parsing was successful, allocate the element
  theElement =
      new BeamContact2D(iData[0], iData[1], iData[2], iData[3], iData[4],
                        *theMaterial, dData[0], dData[1], dData[2], icSwitch);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type BeamContact2DElement\n";
    return 0;
  }

  return theElement;
}


#include <element/community/UWelements/BeamContact2Dp.h>
static Element*
TclDispatch_BeamContact2Dp(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  if (num_BeamContact2Dp == 0) {
    num_BeamContact2Dp++;
    opslog << "BeamContact2Dp element - Written: C.McGann, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to an element that will be returned
  Element *theElement = 0;

  int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingInputArgs < 7) {
    opserr << "Invalid #args, want: element BeamContact2Dp eleTag? iNode? "
              "jNode? secondaryNode? matTag? width? penalty? <cSwitch>?\n";
    return 0;
  }

  int iData[5];
  double dData[2];
  int icSwitch = 0;

  int numData = 5;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element BeamContact2Dp "
           << iData[0] << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element BeamContact2Dp " << iData[0]
           << endln;
    return 0;
  }

  int matID = iData[4];
  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == 0) {
    opserr << "WARNING element BeamContact2Dp " << iData[0] << endln;
    opserr << " Material: " << matID << "not found\n";
    return 0;
  }

  numRemainingInputArgs -= 7;
  while (numRemainingInputArgs >= 1) {
    numData = 1;
    if (OPS_GetIntInput(&numData, &icSwitch) != 0) {
      opserr << "WARNING invalid initial contact flag: element BeamContact2Dp "
             << iData[0] << endln;
      return 0;
    }
    numRemainingInputArgs -= 1;
  }

  // Parsing was successful, allocate the element
  theElement = new BeamContact2Dp(iData[0], iData[1], iData[2], iData[3],
                                  *theMaterial, dData[0], dData[1], icSwitch);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type BeamContact2Dp\n";
    return 0;
  }

  return theElement;
}




#include <element/community/UWelements/BeamContact3D.h>
static Element*
TclDispatch_BeamContact3D(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  if (num_BeamContact3D == 0) {
    num_BeamContact3D++;
    opslog << "BeamContact3D element - Written: K.Petek, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to a uniaxial material that will be returned
  Element *theElement = 0;

  int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingInputArgs < 10) {
    opserr << "Invalid #args,  want: element BeamContact3D eleTag?  iNode? "
              "jNode? secondaryNode? lambdaNode? radius? crdTransf? matTag? "
              "tolGap? tolF? <cSwitch>?\n";
    return 0;
  }

  int iData[7];
  double dData[3];
  int icSwitch = 0;

  int numData = 5;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element BeamContact3DElement"
           << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element BeamContact3D " << iData[0]
           << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetIntInput(&numData, &iData[5]) != 0) {
    opserr << "WARNING invalid integer data: element BeamContact3DElement"
           << iData[0] << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
    opserr << "WARNING invalid data: element BeamContact3D " << iData[0]
           << endln;
    return 0;
  }

  int transfTag = iData[5];
  CrdTransf *theTransf = G3_getSafeBuilder(rt)->getTypedObject<CrdTransf>(transfTag);
  if (theTransf == 0) {
    opserr << "WARNING element BeamContact3D " << iData[0] << endln;
    opserr << " coordTransf: " << transfTag << "not found\n";
    return 0;
  }

  int matID = iData[6];
  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == 0) {
    opserr << "WARNING element BeamContact3D " << iData[0] << endln;
    opserr << " Material: " << matID << "not found\n";
    return 0;
  }

  numRemainingInputArgs -= 10;
  while (numRemainingInputArgs >= 1) {
    numData = 1;
    if (OPS_GetIntInput(&numData, &icSwitch) != 0) {
      opserr << "WARNING invalid initial contact flag: element BeamContact3D "
             << iData[0] << endln;
      return 0;
    }
    numRemainingInputArgs -= 1;
  }

  // Parsing was successful, allocate the element
  theElement = new BeamContact3D(iData[0], iData[1], iData[2], iData[3],
                                 iData[4], dData[0], *theTransf, *theMaterial,
                                 dData[1], dData[2], icSwitch);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type BeamContact3DElement\n";
    return 0;
  }

  return theElement;
}




#include <element/community/UWelements/BeamContact3Dp.h>
static Element*
TclDispatch_BeamContact3Dp(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  if (num_BeamContact3Dp == 0) {
    num_BeamContact3Dp++;
    opslog << "BeamContact3Dp element - Written: K.Petek, C.McGann, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to a uniaxial material that will be returned
  Element *theElement = 0;

  int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingInputArgs < 8) {
    opserr << "Invalid #args,  want: element BeamContact3Dp eleTag?  iNode? "
              "jNode? secondaryNode? radius? crdTransf? matTag? penalty? "
              "<cSwitch>?\n";
    return 0;
  }

  int iData[6];
  double dData[2];
  int icSwitch = 0;

  int numData = 4;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element BeamContact3DpElement"
           << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element BeamContact3Dp " << iData[0]
           << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetIntInput(&numData, &iData[4]) != 0) {
    opserr << "WARNING invalid integer data: element BeamContact3DpElement"
           << iData[0] << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
    opserr << "WARNING invalid data: element BeamContact3Dp " << iData[0]
           << endln;
    return 0;
  }

  int transfTag = iData[4];
  CrdTransf *theTransf = G3_getSafeBuilder(rt)->getTypedObject<CrdTransf>(transfTag);
  if (theTransf == 0) {
    opserr << "WARNING element BeamContact3Dp " << iData[0] << endln;
    opserr << " coordTransf: " << transfTag << "not found\n";
    return 0;
  }

  int matID = iData[5];
  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == 0) {
    opserr << "WARNING element BeamContact3Dp " << iData[0] << endln;
    opserr << " Material: " << matID << "not found\n";
    return 0;
  }

  numRemainingInputArgs -= 8;
  while (numRemainingInputArgs >= 1) {
    numData = 1;
    if (OPS_GetIntInput(&numData, &icSwitch) != 0) {
      opserr << "WARNING invalid initial contact flag: element BeamContact3Dp "
             << iData[0] << endln;
      return 0;
    }
    numRemainingInputArgs -= 1;
  }

  // Parsing was successful, allocate the element
  theElement =
      new BeamContact3Dp(iData[0], iData[1], iData[2], iData[3], dData[0],
                         *theTransf, *theMaterial, dData[1], icSwitch);

  if (theElement == 0) {
    opserr
        << "WARNING could not create element of type BeamContact3DpElement\n";
    return 0;
  }

  return theElement;
}




#include <element/community/UWelements/BeamEndContact3D.h>
static Element*
TclDispatch_BeamEndContact3D(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  if (num_BeamEndContact3D == 0) {
    num_BeamEndContact3D++;
    opslog << "BeamEndContact3D element - Written: C.McGann, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to an element that will be returned
  Element *theElement = 0;

  int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingInputArgs < 8) {
    opserr << "Invalid #args, want: element BeamEndContact3D eleTag? iNode? "
              "jNode? secondaryNode? lambdaNode? radius? gapTol? forceTol "
              "<cFlag>?\n";
    return 0;
  }

  int iData[5];
  double dData[3];
  int icSwitch = 0;

  int numData = 5;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element BeamEndContact3D "
           << iData[0] << endln;
    return 0;
  }

  numData = 3;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid double data: element BeamEndContact3D "
           << iData[0] << endln;
    return 0;
  }

  numRemainingInputArgs -= 8;
  while (numRemainingInputArgs >= 1) {
    numData = 1;
    if (OPS_GetIntInput(&numData, &icSwitch) != 0) {
      opserr
          << "WARNING invalid initial contact flag: element BeamEndContact3D "
          << iData[0] << endln;
      return 0;
    }
    numRemainingInputArgs -= 1;
  }

  // Parsing was successful, allocate the element
  theElement =
      new BeamEndContact3D(iData[0], iData[1], iData[2], iData[3], iData[4],
                           dData[0], dData[1], dData[2], icSwitch);

  if (theElement == 0) {
    opserr
        << "WARNING could not create element of type BeamEndContact3DElement\n";
    return 0;
  }

  return theElement;
}




#include <element/community/UWelements/BeamEndContact3Dp.h>
static Element*
TclDispatch_BeamEndContact3Dp(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  if (num_BeamEndContact3Dp == 0) {
    num_BeamEndContact3Dp++;
    opslog << "BeamEndContact3Dp element - Written: C.McGann, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to an element that will be returned
  Element *theElement = 0;

  int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingInputArgs < 6) {
    opserr << "Invalid #args, want: element BeamEndContact3Dp eleTag? iNode? "
              "jNode? sNode? radius? penalty? <cFlag>?\n";
    return 0;
  }

  int iData[4];
  double dData[2];
  int icSwitch = 0;

  int numData = 4;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element BeamEndContact3Dp "
           << iData[0] << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid double data: element BeamEndContact3Dp "
           << iData[0] << endln;
    return 0;
  }

  numRemainingInputArgs -= 6;
  while (numRemainingInputArgs >= 1) {
    numData = 1;
    if (OPS_GetIntInput(&numData, &icSwitch) != 0) {
      opserr
          << "WARNING invalid initial contact flag: element BeamEndContact3Dp "
          << iData[0] << endln;
      return 0;
    }
    numRemainingInputArgs -= 1;
  }

  // Parsing was successful, allocate the element
  theElement = new BeamEndContact3Dp(iData[0], iData[1], iData[2], iData[3],
                                     dData[0], dData[1], icSwitch);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type "
              "BeamEndContact3DpElement\n";
    return 0;
  }

  return theElement;
}




#include <element/community/UWelements/Brick8FiberOverlay.h>
static Element *OPS_Brick8FiberOverlay(void)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  if (num_Brick8FiberOverlay == 0) {
    num_Brick8FiberOverlay++;
    opslog << "Brick8FiberOverlay element - Written: M.Chiaramonte, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  Element *theElement = 0;

  if (OPS_GetNumRemainingInputArgs() != 15) {
    opserr << "Want: Brick8FiberOverlay tag? nd1? nd2? nd3? nd4? nd5? nd6? "
              "nd7? nd8? matTag? AreaFiber? B1? B2? B3? B4?\n";
    return 0;
  }

  int iData[9];
  double dData[5];
  int matTag = 0;
  int eleTag = 0;
  int numData = 9;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element Brick8FiberOverlay"
           << endln;
    return 0;
  }
  eleTag = iData[0];

  numData = 1;
  if (OPS_GetIntInput(&numData, &matTag) != 0) {
    opserr << "WARNING element Brick8FiberOverlay: invalid matTag for element: "
           << eleTag << "\n";
    return 0;
  }
  numData = 5;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element Brick8FiberOverlay " << eleTag
           << endln;
    return 0;
  }

  UniaxialMaterial *theMaterial = OPS_GetUniaxialMaterial(matTag);

  if (theMaterial == 0) {
    opserr << "WARNING material with tag " << matTag << "not found for element "
           << eleTag << endln;
    return 0;
  }

  // Parsing was successful, allocate the material
  theElement = new Brick8FiberOverlay(iData[0], iData[1], iData[2], iData[3],
                                      iData[4], iData[5], iData[6], iData[7],
                                      iData[8], *theMaterial, dData[0],
                                      dData[1], dData[2], dData[3], dData[4]);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type Brick8FiberOverlay\n";
    return 0;
  }

  return theElement;
}




#include <element/community/UWelements/EmbeddedBeamInterfaceL.h>
static Element *OPS_EmbeddedBeamInterfaceL(void)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  if (num_EmbeddedBeamInterfaceL == 0) {
    num_EmbeddedBeamInterfaceL++;
    opslog << "EmbeddedBeamInterfaceL element - Written: A.Ghofrani, "
              "D.Turello, P.Arduino, U.Washington\n";
  }

  Element *theElement = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 1) {
    opserr << "Want: EmbeddedBeamInterfaceL tag? \n";
    return 0;
  }

  int iData[1];
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element EmbeddedBeamInterfaceL"
           << endln;
    return 0;
  }

  int eleTag = iData[0];

  theElement = new EmbeddedBeamInterfaceL(eleTag);

  if (theElement == 0) {
    opserr
        << "WARNING could not create element of type EmbeddedBeamInterfaceL\n";
    return 0;
  }

  return theElement;
}




#include <element/community/UWelements/EmbeddedBeamInterfaceP.h>
static Element *OPS_EmbeddedBeamInterfaceP(void)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  if (num_EmbeddedBeamInterfaceP == 0) {
    num_EmbeddedBeamInterfaceP++;
    opslog << "EmbeddedBeamInterfaceP element - Written: A.Ghofrani, "
              "D.Turello, P.Arduino, U.Washington\n";
  }

  Element *theElement = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 1) {
    opserr << "Want: EmbeddedBeamInterfaceP tag? \n";
    return 0;
  }

  int iData[1];
  int eleTag = 0;
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element EmbeddedBeamInterfaceP"
           << endln;
    return 0;
  }

  eleTag = iData[0];

  theElement = new EmbeddedBeamInterfaceP(eleTag);

  if (theElement == 0) {
    opserr
        << "WARNING could not create element of type EmbeddedBeamInterfaceP\n";
    return 0;
  }

  return theElement;
}




#include <element/community/UWelements/EmbeddedEPBeamInterface.h>
static Element *OPS_EmbeddedEPBeamInterface(void)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  if (num_EmbeddedEPBeamInterface == 0) {
    num_EmbeddedEPBeamInterface++;
    opslog << "EmbeddedEPBeamInterface element - Written: A.Ghofrani, "
              "D.Turello, P.Arduino, U.Washington\n";
  }

  Element *theElement = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 1) {
    opserr << "Want: EmbeddedEPBeamInterface tag? \n";
    return 0;
  }

  int iData[1];
  int eleTag = 0;
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element EmbeddedEPBeamInterface"
           << endln;
    return 0;
  }

  eleTag = iData[0];

  theElement = new EmbeddedEPBeamInterface(eleTag);

  if (theElement == 0) {
    opserr
        << "WARNING could not create element of type EmbeddedEPBeamInterface\n";
    return 0;
  }

  return theElement;
}




#include <element/community/UWelements/PileToe3D.h>
static Element*
TclDispatch_PileToe3D(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  if (num_PileToe3D == 0) {
    num_PileToe3D++;
    // U.Washington\n", 1);
    opslog << "PileToe3D element - Written: P.Arduino, P.Mackenzie-Helnwein, "
              "U.Washington\n";
  }

  // Pointer to a uniaxial material that will be returned
  Element *theElement = 0;

  int numRemainingInputArgs = OPS_GetNumRemainingInputArgs();

  if (numRemainingInputArgs < 7) {
    opserr << "Invalid #args,  want: element PileToe3D eleTag?  iNode? BiNode? "
              "BjNode? radius? k? crdTransf?\n";
    return 0;
  }

  int iData[5];
  double dData[2];

  int numData = 4;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element PileToe3D" << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid radius data: element PileToe3D " << iData[0]
           << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetDoubleInput(&numData, &dData[1]) != 0) {
    opserr << "WARNING invalid  k data: element PileToe3D " << iData[0]
           << endln;
    return 0;
  }

  numData = 1;
  if (OPS_GetIntInput(&numData, &iData[4]) != 0) {
    opserr << "WARNING invalid integer crdTransf data: element PileToe3D"
           << iData[0] << endln;
    return 0;
  }

  int transfTag = iData[4];
  CrdTransf *theTransf = G3_getSafeBuilder(rt)->getTypedObject<CrdTransf>(transfTag);
  if (theTransf == 0) {
    opserr << "WARNING element PileToe3D " << iData[0] << endln;
    opserr << " coordTransf: " << transfTag << "not found\n";
    return 0;
  }

  // Parsing was successful, allocate the element
  theElement = new PileToe3D(iData[0], iData[1], iData[2], iData[3], dData[0],
                             dData[1], *theTransf);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type PileToe3D\n";
    return 0;
  }

  return theElement;
}




#include <element/community/UWelements/Quad4FiberOverlay.h>
static Element *OPS_Quad4FiberOverlay(void)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  if (num_Quad4FiberOverlay == 0) {
    num_Quad4FiberOverlay++;
    opslog << "Quad4FiberOverlay element - Written: M.Chiaramonte, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  Element *theElement = 0;

  if (OPS_GetNumRemainingInputArgs() != 9) {
    opserr << "Want: Quad4FiberOverlay tag? nd1? nd2? nd3? nd4? matTag? "
              "CrossSectionArea? B1?  B2? \n";
    return 0;
  }

  int iData[6];
  double dData[3];
  int matTag = 0;
  int eleTag = 0;
  int numData = 0;
  numData = 5;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element Quad4FiberOverlay"
           << endln;
    return 0;
  }

  eleTag = iData[0];

  numData = 1;
  if (OPS_GetIntInput(&numData, &matTag) != 0) {
    opserr << "WARNING element Quad4FiberOverlay: invalid matTag for element: "
           << eleTag << "\n";
    return 0;
  }
  numData = 3;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element Quad4FiberOverlay " << eleTag
           << endln;
    return 0;
  }

  UniaxialMaterial *theMaterial = OPS_GetUniaxialMaterial(matTag);

  if (theMaterial == 0) {
    opserr << "WARNING material with tag " << matTag << "not found for element "
           << eleTag << endln;
    return 0;
  }

  // Parsing was successful, allocate the material
  theElement =
      new Quad4FiberOverlay(iData[0], iData[1], iData[2], iData[3], iData[4],
                            *theMaterial, dData[0], dData[1], dData[2]);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type Quad4FiberOverlay\n";
    return 0;
  }

  return theElement;
}




#include <element/community/UWelements/QuadBeamEmbedContact.h>
static Element *OPS_QuadBeamEmbedContact(void)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  if (num_QuadBeamEmbedContact == 0) {
    num_QuadBeamEmbedContact++;
    opslog << "QuadBeamEmbedContact element - Written: A.Ghofrani, P.Arduino, "
              "U.Washington\n";
  }

  Element *theElement = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs < 10) {
    opserr << "Want: QuadBeamEmbedContact tag? Qnd1? Qnd2? Qnd3? Qnd4? Bnd1? "
              "Bnd2? radius? fricCoeff? normalPenalty? <tangentialPenalty?> \n";
    return 0;
  }

  int iData[7];
  int eleTag = 0;
  int numData = 7;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element QuadBeamEmbedContact"
           << endln;
    return 0;
  }

  numData = 3;
  double dData[3];
  double oData[1];

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element QuadBeamEmbedContact" << endln;
    return 0;
  }

  oData[0] = dData[2];
  numData = numArgs - 10;
  if (numData != 0)
    if (OPS_GetDouble(&numData, oData) != 0) {
      opserr << "WARNING invalid data: element QuadBeamEmbedContact" << endln;
      return 0;
    }

  eleTag = iData[0];

  theElement = new QuadBeamEmbedContact(iData[0], iData[1], iData[2], iData[3],
                                        iData[4], iData[5], iData[6], dData[0],
                                        dData[1], dData[2], oData[0]);

  if (theElement == 0) {
    opserr << "WARNING could not create element of type QuadBeamEmbedContact\n";
    return 0;
  }

  return theElement;
}




#include <element/community/UWelements/SimpleContact2D.h>
static Element*
TclDispatch_SimpleContact2D(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  if (num_SimpleContact2D == 0) {
    num_SimpleContact2D++;
    // OPS_Error("SimpleContact2D element - Written: K.Petek, P.Arduino,
    // P.Mackenzie-Helnwein, U.Washington\n", 1);
    opslog << "SimpleContact2D element - Written: K.Petek, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to a uniaxial material that will be returned
  Element *theElement = 0;

  if (OPS_GetNumRemainingInputArgs() != 8) {
    opserr << "Invalid #args,  want: element SimpleContact2D eleTag? iNode? "
              "jNode? secondaryNode? lambdaNode? matTag? tolGap? tolForce?\n";
    return 0;
  }

  int iData[6];
  double dData[2];

  int numData = 6;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element SimpleContact2DElement"
           << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element SimpleContact2D " << iData[0]
           << endln;
    return 0;
  }

  int matID = iData[5];
  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == 0) {
    opserr << "WARNING element SimpleContact2D " << iData[0] << endln;
    opserr << " Material: " << matID << "not found\n";
    return 0;
  }

  // Parsing was successful, allocate the material
  theElement = new SimpleContact2D(iData[0], iData[1], iData[2], iData[3],
                                   iData[4], *theMaterial, dData[0], dData[1]);

  if (theElement == 0) {
    opserr
        << "WARNING could not create element of type SimpleContact2DElement\n";
    return 0;
  }

  return theElement;
}




#include <element/community/UWelements/SimpleContact3D.h>
static Element*
TclDispatch_SimpleContact3D(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  static int num_SimpleContact3D;
  if (num_SimpleContact3D == 0) {
    num_SimpleContact3D++;
    opslog << "SimpleContact3D element - Written: K.Petek, P.Arduino, "
              "P.Mackenzie-Helnwein, U.Washington\n";
  }

  // Pointer to a uniaxial material that will be returned
  Element *theElement = 0;

  if (OPS_GetNumRemainingInputArgs() != 10) {
    opserr << "Invalid #args,  want: element SimpleContact3D eleTag? iNode? "
              "jNode? kNode? lNode? secondaryNode? lambdaNode? matTag? tolGap? "
              "tolForce?\n";
    return 0;
  }

  int iData[8];
  double dData[2];

  int numData = 8;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer data: element SimpleContact3DElement"
           << endln;
    return 0;
  }

  numData = 2;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid data: element SimpleContact3D " << iData[0]
           << endln;
    return 0;
  }

  int matID = iData[7];
  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == 0) {
    opserr << "WARNING element SimpleContact3D " << iData[0] << endln;
    opserr << " Material: " << matID << "not found\n";
    return 0;
  }

  // Parsing was successful, allocate the material
  theElement =
      new SimpleContact3D(iData[0], iData[1], iData[2], iData[3], iData[4],
                          iData[5], iData[6], *theMaterial, dData[0], dData[1]);

  if (theElement == 0) {
    opserr
        << "WARNING could not create element of type SimpleContact3DElement\n";
    return 0;
  }

  return theElement;
}
#endif
