#include <tcl.h>
#include <assert.h>
#include <Logging.h>
#include <Parsing.h>
#include <BasicModelBuilder.h>
#include <element/Shell/ASDShellQ4.h>
#include <element/Shell/ShellANDeS.h>
#include <element/Shell/ShellDKGQ.h>
#include <element/Shell/ShellDKGT.h>
#include <element/Shell/ShellMITC4.h>
#include <element/Shell/ShellMITC9.h>
#include <element/Shell/ShellNLDKGQ.h>
#include <element/Shell/ShellMITC4Thermal.h>
#include <element/Shell/ShellNLDKGQThermal.h>
#include <element/Shell/ShellNLDKGT.h>
using namespace OpenSees;

Element*
TclDispatch_newASDShellQ4(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);

  Vector3D local_x{};
  bool corotational = false;
  bool use_eas = true;
  bool use_drill_stab = false;
  double drilling_stab = 0.01;
  ASDShellQ4::DrillingDOFMode drill_mode = ASDShellQ4::DrillingDOF_Elastic;

  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  if (argc < 6) {
    opserr << "Want: element ASDShellQ4 $tag $iNode $jNode $kNode $lNode "
              "$secTag <-corotational>";
    return nullptr;
  }

  int tag;
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid element tag\n";
    return nullptr;
  }
  int nodes[4];
  for (int i=3; i<7; i++) {
    if (Tcl_GetInt(interp, argv[i], &nodes[i-3]) != TCL_OK) {
      opserr << "WARNING invalid node tag \"" << argv[i] << "\" \n";
      return nullptr;
    }
  }
  int stag;
  if (Tcl_GetInt(interp, argv[7], &stag) != TCL_OK) {
    opserr << "WARNING invalid section tag\n";
    return nullptr;
  }
  for (int i=6; i<argc; i++) {
      if ((strcmp(argv[i], "-corotational") == 0) ||
          (strcmp(argv[i], "-Corotational") == 0))
        corotational = true;
    
      else if (strcmp(argv[i], "-noeas") == 0) {
          use_eas = false;
      }
      else if (strcmp(argv[i], "-drillingStab") == 0) {
          if (drill_mode != ASDShellQ4::DrillingDOF_Elastic) {
              opserr << "Error: element ASDShellQ4: -drillingStab and -drillingNL options are mutually exclusive\n";
              return 0;
          }
          if (argc < i + 2) {
              opserr << "Error: drilling stabilization parameter not provided with -drillingStab option\n";
              return nullptr;
          }
          if (Tcl_GetDouble(interp, argv[i+1], &drilling_stab) != TCL_OK) {
              opserr << "Error: cannot get drilling stabilization parameter with -drillingStab option\n";
              return nullptr;
          }
          drilling_stab = std::max(0.0, std::min(1.0, drilling_stab));
          drill_mode = ASDShellQ4::DrillingDOF_Elastic;
          use_drill_stab = true;
          i++;
      }
      else if (strcmp(argv[i], "-drillingNL") == 0) {
          if (use_drill_stab) {
              opserr << "Error: element ASDShellQ4: -drillingStab and -drillingNL options are mutually exclusive\n";
              return 0;
          }
          drill_mode = ASDShellQ4::DrillingDOF_NonLinear;
          drilling_stab = 1.0;
      }
      // else if (strcmp(argv[i], "-local") == 0) {
      //     if (OPS_GetNumRemainingInputArgs() < 3) {
      //         opserr << "Error: element ASDShellQ4: not enough arguments for -local options (3 components are required)\n";
      //         return 0;
      //     }
      //     for (int i = 0; i < 3; ++i) {
      //         double local_x_com;
      //         if (OPS_GetDoubleInput(&numData, &local_x_com) == 0) {
      //             local_x(i) = local_x_com;
      //         }
      //         else {
      //             opserr << "Error: element ASDShellQ4: cannot get the component " << i + 1 << " for the local X axis\n";
      //             return 0;
      //         }
      //     }
      // }
  }

  SectionForceDeformation *section = builder->getTypedObject<SectionForceDeformation>(stag);
  if (section == nullptr)
    return nullptr;

  return new ASDShellQ4(tag, nodes[0], nodes[1], nodes[2], nodes[3],
                        section, local_x, 
                        corotational, use_eas, drill_mode, drilling_stab);
}


Element*
TclDispatch_newShellANDeS(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{

  if (argc < 6) {
    opserr << "Want: element ShellANDeS $tag $iNode $jNode $kNode $thick $E "
              "$nu $rho";
    return nullptr;
  }

  int numArgs = OPS_GetNumRemainingInputArgs();

  int iData[4];
  int numData = 4;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag\n";
    return nullptr;
  }

  double dData[11];
  numArgs = OPS_GetNumRemainingInputArgs();
  if (OPS_GetDoubleInput(&numArgs, dData) != 0) {
    opserr << "WARNING invalid double thickness: element ShellANDeS \n";
    return nullptr;
  }

  Element *theElement = nullptr;

  if (numArgs == 4) {
    theElement = new ShellANDeS(iData[0], iData[1], iData[2], iData[3],
                                dData[0], dData[1], dData[2], dData[3]);
  } else if (numArgs == 11) {
    theElement =
        new ShellANDeS(iData[0], iData[1], iData[2], iData[3], dData[0],
                       dData[1], dData[2], dData[3], dData[4], dData[5],
                       dData[6], dData[7], dData[8], dData[9], dData[10]);
  }

  return theElement;
}

Element*
TclDispatch_newShellDKGQ(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  if (argc < 6) {
    opserr << "Want: element ShellDKGQ $tag $iNode $jNoe $kNode $lNode $secTag";
    return nullptr;
  }

  int iData[6];
  int numData = 6;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag\n";
    return nullptr;
  }

  SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(iData[5]);
  if (theSection == nullptr)
    return nullptr;

  return new ShellDKGQ(iData[0], iData[1], iData[2], iData[3], iData[4],
                             *theSection);
}

Element*
TclDispatch_newShellDKGT(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 5) {
    opserr << "Want: element ShellDKGT $tag $iNode $jNoe $kNode $secTag";
    return nullptr;
  }

  int iData[5];
  int numData = 5;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag\n";
    return nullptr;
  }

  SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(iData[4]);
  if (theSection == nullptr)
    return nullptr;


  double b_data[3] = {0, 0, 0};

  int num_remaining_args = OPS_GetNumRemainingInputArgs();

  if (num_remaining_args > 3) {
    num_remaining_args = 3;
  }
  if (num_remaining_args > 0) {
    if (OPS_GetDoubleInput(&num_remaining_args, b_data) < 0) {
      opserr << "WARNING: invalid double b_data\n";
      return nullptr;
    }
  }

  return new ShellDKGT(iData[0], iData[1], iData[2], iData[3],
                       *theSection, b_data[0], b_data[1], b_data[2]);
}



int
TclDispatch_newShellMITC4(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  bool updateBasis = false;
  Element *theElement = nullptr;

  if (argc < 6) {
    opserr << "Want: element ShellMITC4 $tag $iNode $jNode $kNode $lNode "
              "$secTag<-updateBasis>";
    return TCL_ERROR;
  }

  int iData[6];
  int numData = 6;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag\n";
    return TCL_ERROR;
  }

  if (argc == 7) {
    const char *type = OPS_GetString();
    if (strcmp(type, "-updateBasis") == 0)
      updateBasis = true;
  }

  SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(iData[5]);
  if (theSection == nullptr)
    return TCL_ERROR;

  theElement = new ShellMITC4(iData[0], iData[1], iData[2], iData[3], iData[4],
                              *theSection, updateBasis);

  if (builder->getDomain()->addElement(theElement) == false)
    return TCL_ERROR;
  return TCL_OK;
}



Element*
TclDispatch_newShellMITC4Thermal(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 6) {
    opserr << "Want: element ShellMITC4Thermal $tag $iNode $jNoe $kNode $lNode "
              "$secTag";
    return nullptr;
  }

  int iData[6];
  int numData = 6;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag\n";
    return nullptr;
  }

  SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(iData[5]);
  if (theSection == nullptr)
    return nullptr;

  return new ShellMITC4Thermal(iData[0], iData[1], iData[2], iData[3], iData[4], *theSection);
}


Element*
TclDispatch_newShellMITC9(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 11) {
    opserr << "Want: element ShellMITC9 $tag $node1 $node2 .... $node9 $secTag";
    return nullptr;
  }

  int iData[11];
  int numData = 11;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag\n";
    return nullptr;
  }

  SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(iData[10]);
  if (theSection == nullptr)
    return nullptr;

  return
      new ShellMITC9(iData[0], iData[1], iData[2], iData[3], iData[4], iData[5],
                     iData[6], iData[7], iData[8], iData[9], *theSection);
}



Element*
TclDispatch_newShellNLDKGQ(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 6) {
    opserr
        << "Want: element ShellNLDKGQ $tag $iNode $jNoe $kNode $lNode $secTag";
    return nullptr;
  }

  int iData[6];
  int numData = 6;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag\n";
    return nullptr;
  }

  SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(iData[5]);
  if (theSection == nullptr)
    return nullptr;

  return new ShellNLDKGQ(iData[0], iData[1], iData[2], iData[3], iData[4],
                               *theSection);

}


Element*
TclDispatch_newShellNLDKGQThermal(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 6) {
    opserr << "Want: element ShellNLDKGQThermal $tag $iNode $jNoe $kNode "
              "$lNode $secTag";
    return nullptr;
  }

  int iData[6];
  int numData = 6;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag\n";
    return nullptr;
  }

  SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(iData[5]);
  if (theSection == nullptr)
    return nullptr;

  return new ShellNLDKGQThermal(iData[0], iData[1], iData[2], iData[3],
                                      iData[4], *theSection);

}


Element*
TclDispatch_newShellNLDKGT(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 5) {
    opserr << "Want: element ShellNLDKGT $tag $iNode $jNoe $kNode $secTag";
    return nullptr;
  }

  int iData[5];
  int numData = 5;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag\n";
    return nullptr;
  }

  SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(iData[4]);
  if (theSection == nullptr)
    return nullptr;

  return
      new ShellNLDKGT(iData[0], iData[1], iData[2], iData[3], *theSection);

}
