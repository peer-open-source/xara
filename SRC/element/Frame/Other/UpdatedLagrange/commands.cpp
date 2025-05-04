//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
#include <tcl.h>
#include <stdlib.h>
#include <string.h>
#include <OPS_Stream.h>
#include <Logging.h>
#include <Parsing.h>

#include <Domain.h>
#include <Node.h>
#include <Matrix.h>

#include <CyclicModel.h>
#include <LinearCyclic.h>
#include <BilinearCyclic.h>
#include <QuadraticCyclic.h>
#include <BasicModelBuilder.h>

int
TclBasicBuilder_addLinearCylic(ClientData clientData, Tcl_Interp *interp, int argc,
                               TCL_Char ** const argv)
{
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  int tag;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid CyclicModel tag" << "\n";
    return TCL_ERROR;
  }

  CyclicModel *cModel = new LinearCyclic(tag);
  if (!cModel) {
    opserr << "TclBasicBuilder_addLinearCycylic - could not allocate memory\n";
    return TCL_ERROR;
  }
  if (builder->addTaggedObject<CyclicModel>(*cModel) < 0) {
    opserr << "WARNING TclElmtBuilder - could not add cycModel to domain ";
    opserr << tag << "\n";
    opserr << "\a";
    return TCL_ERROR;
  }

  return TCL_OK;
}

int
TclBasicBuilder_addBilinearCyclic(ClientData clientData, Tcl_Interp *interp,
                                  int argc, TCL_Char ** const argv)
{
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  int tag;
  double wt;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid CyclicModel tag" << "\n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &wt) != TCL_OK) {
    opserr << "WARNING invalid arg[3]" << "\n";
    return TCL_ERROR;
  }

  CyclicModel *cModel = new BilinearCyclic(tag, wt);
  if (builder->addTaggedObject<CyclicModel>(*cModel) < 0) {
    opserr << "WARNING TclElmtBuilder - could not add cycModel to domain ";
    opserr << tag << "\n";
    opserr << "\a";
    return TCL_ERROR;
  }

  return TCL_OK;
}

int
TclBasicBuilder_addQuadraticCyclic(ClientData clientData, Tcl_Interp *interp,
                                   int argc, TCL_Char ** const argv)
{
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  int tag;
  double wt, qy;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid CyclicModel tag" << "\n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &wt) != TCL_OK) {
    opserr << "WARNING invalid arg[3]" << "\n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[4], &qy) != TCL_OK) {
    opserr << "WARNING invalid arg[4]" << "\n";
    return TCL_ERROR;
  }

  CyclicModel *cModel = new QuadraticCyclic(tag, wt, qy);
  if (builder->addTaggedObject<CyclicModel>(*cModel) < 0) {
    opserr << "WARNING TclElmtBuilder - could not add cycModel to domain ";
    opserr << tag << "\n";
    opserr << "\a";
    return TCL_ERROR;
  }
  return TCL_OK;
}

/*******************************************************************************************/
int
TclCommand_addCyclicModel(ClientData clientData, Tcl_Interp *interp,
                                  int argc, TCL_Char ** const argv)
{

  if (strcmp(argv[1], "linear") == 0) {
    int result = TclBasicBuilder_addLinearCylic(clientData, interp, argc, argv);
    return result;

  } else if (strcmp(argv[1], "bilinear") == 0) {
    int result = TclBasicBuilder_addBilinearCyclic(clientData, interp, argc, argv);
    return result;
  }

  else if (strcmp(argv[1], "quadratic") == 0) {
    int result = TclBasicBuilder_addQuadraticCyclic(clientData, interp, argc, argv);
    return result;
  }

  else

    return TCL_ERROR;
}

#include <tcl.h>
#include <stdlib.h>
#include <string.h>
#include <OPS_Stream.h>

#include <Domain.h>
#include <Node.h>
#include <Matrix.h>

#include <Elastic2DGNL.h>

#include <BasicModelBuilder.h>

#define tcl_debug 0

// Elastic2DGNL(int tag, double A, double E, double I, int Nd1, int Nd2,
//             double rho = 0.0, bool islinear = false);

int
TclBasicBuilder_addElastic2dGNL(ClientData clientData, Tcl_Interp *interp, int argc,
                                TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);

  if (tcl_debug)
    opserr << " TclBasicBuilder_addElastic2dGNL \n";

  if (argc < 8) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "element element2dGNL int tag, int Nd1, int Nd2, double A, "
              "double E, double Iz, <int linear>\n";

    return TCL_ERROR;
  }

  int tag, ndI, ndJ;
  double E, A, I;
//double massDens = 0.0;
  bool linear = false;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid Elastic2dGNL tag" << "\n";
    return TCL_ERROR;
  }
  if (tcl_debug)
    opserr << "\tElement tag = " << tag << "\n";

  if (Tcl_GetInt(interp, argv[3], &ndI) != TCL_OK) {
    opserr << "WARNING invalid node I\n";
    opserr << "Elastic2dGNL: " << tag << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[4], &ndJ) != TCL_OK) {
    opserr << "WARNING invalid node J\n";
    opserr << "Elastic2dGNL: " << tag << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[5], &A) != TCL_OK) {
    opserr << "WARNING invalid A\n";
    opserr << "Elastic2dGNL: " << tag << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[6], &E) != TCL_OK) {
    opserr << "WARNING invalid E\n";
    opserr << "Elastic2dGNL: " << tag << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[7], &I) != TCL_OK) {
    opserr << "WARNING invalid I\n";
    opserr << "Elastic2dGNL: " << tag << "\n";
    return TCL_ERROR;
  }

  if (argc == 9) {
    int lin = 0;
    if (Tcl_GetInt(interp, argv[8], &lin) != TCL_OK) {
      opserr << "WARNING invalid Linear Flag\n";
      opserr << "Elastic2dGNL: " << tag << "\n";
      return TCL_ERROR;
    }

    if (lin == 1)
      linear = true;

    if (tcl_debug)
      opserr << " 9 arguments - " << lin << "\n";
  }

  // if(tcl_debug) opserr << "\tAdded upto mass - input parameters\n";

  Element *theElement =
      new Elastic2dGNL(tag, A, E, I, ndI, ndJ, linear); //, false, massDens);

  if (tcl_debug)
    opserr << "\tElement created\n";

  // Ensure we have created the element, out of memory if got here and no
  // element
  if (theElement == 0) {
    opserr << "WARNING ran out of memory creating element\n";
    opserr << "Elastic2dGNL: " << tag << "\n";
    opserr << "\a";
    return TCL_ERROR;
  }

  if (builder->getDomain()->addElement(theElement) == false) {
    opserr << "WARNING TclElmtBuilder - addElastic2dGNL - could not add "
              "element to domain ";
    opserr << tag << "\n";
    opserr << "\a";
    return TCL_ERROR;
  }

  if (tcl_debug)
    opserr << "\tElement number " << tag << " added to domain - returning\n";

  return TCL_OK;
}

#include <string.h>
#include <OPS_Stream.h>

#include <Domain.h>
#include <Node.h>
#include <Matrix.h>

#include <CyclicModel.h>
#include <Inelastic2DYS01.h>
#include <Inelastic2DYS02.h>
#include <Inelastic2DYS03.h>
//#include <Inelastic2DYS04.h>
//#include <Inelastic2DYS05.h>

#include <YieldSurface_BC.h>
#include <BasicModelBuilder.h>

#define tcl_debug 0

// Element2dGNL(int tag, double A, double E, double I, int Nd1, int Nd2,
//             double rho = 0.0, bool islinear = false);

static int
TclBasicBuilder_addElement2dYS01(ClientData clientData, Tcl_Interp *interp,
                                 int argc, TCL_Char ** const argv)
{
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  if (tcl_debug)
    opserr << " TclBasicBuilder_addElement2dGNL \n";

  if (argc < 11) {
    opserr << "WARNING insufficient arguments\n";
    opserr
        << "element element2dYS tag? Nd1? Nd2? A? E? Iz? ysID1? ysID2? algo?";

    return TCL_ERROR;
  }

  int tag, ndI, ndJ;
  double E, A, I;
  //	double massDens = 0.0;
  int ysID1, ysID2;
  int rf_algo;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid element2dYS tag" << "\n";
    return TCL_ERROR;
  }
  if (tcl_debug)
    opserr << "\tElement tag = " << tag << "\n";

  if (Tcl_GetInt(interp, argv[3], &ndI) != TCL_OK) {
    opserr << "WARNING invalid node I\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[4], &ndJ) != TCL_OK) {
    opserr << "WARNING invalid node J\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[5], &A) != TCL_OK) {
    opserr << "WARNING invalid A\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[6], &E) != TCL_OK) {
    opserr << "WARNING invalid E\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[7], &I) != TCL_OK) {
    opserr << "WARNING invalid I\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[8], &ysID1) != TCL_OK) {
    opserr << "WARNING invalid ysID1\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[9], &ysID2) != TCL_OK) {
    opserr << "WARNING invalid ysID2\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[10], &rf_algo) != TCL_OK) {
    opserr << "WARNING invalid ysID1\n";
    return TCL_ERROR;
  }

  YieldSurface_BC *theYS1 = builder->getTypedObject<YieldSurface_BC>(ysID1);
  if (theYS1 == 0) {
    opserr << "WARNING element2dYS: " << tag << "\n";
    opserr << " no yield surface exists with tag: " << ysID1 << "\n";
    return TCL_ERROR;
  }

  YieldSurface_BC *theYS2 = builder->getTypedObject<YieldSurface_BC>(ysID2);
  if (theYS2 == 0) {
    return TCL_ERROR;
  }

  // 		Inelastic2DYS(	int tag, double A, double E, double I, int Nd1,
  // int Nd2, 						YieldSurface_BC *ysEnd1,
  // YieldSurface_BC *ysEnd2, 						int rf_algo = -1,
  // // updated
  Element *theElement =
      new Inelastic2DYS01(tag, A, E, I, ndI, ndJ, theYS1, theYS2, rf_algo);

  if (tcl_debug)
    opserr << "\tElement created\n";

  // Ensure we have created the element, out of memory if got here and no
  // element
  if (theElement == 0) {
    opserr << "WARNING ran out of memory creating element\n";
    opserr << "\a";
    return TCL_ERROR;
  }

  if (builder->getDomain()->addElement(theElement) == false) {
    opserr << "WARNING TclElmtBuilder - addelement2dYS - could not add element "
              "to domain ";
    opserr << tag << "\n";
    opserr << "\a";
    return TCL_ERROR;
  }

  if (tcl_debug)
    opserr << "\tElement number " << tag << " added to domain - returning\n";

  return TCL_OK;
}

int
TclBasicBuilder_addElement2dYS02(ClientData clientData, Tcl_Interp *interp,
                                 int argc, TCL_Char ** const argv)
{
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  if (tcl_debug)
    opserr << " TclBasicBuilder_addElement2dGNL \n";

  if (argc < 14) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "element element2dYS tag? Nd1? Nd2? A? E? Iz? ysID1? ysID2? "
              "cycType? wt? power? algo?";

    return TCL_ERROR;
  }

  int tag, ndI, ndJ;
  double E, A, I;
  //	double massDens = 0.0;
  int ysID1, ysID2;
  int cyc_type;
  //	double wt;

  int rf_algo = -1;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid element2dYS tag" << "\n";
    return TCL_ERROR;
  }
  if (tcl_debug)
    opserr << "\tElement tag = " << tag << "\n";

  if (Tcl_GetInt(interp, argv[3], &ndI) != TCL_OK) {
    opserr << "WARNING invalid node I\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[4], &ndJ) != TCL_OK) {
    opserr << "WARNING invalid node J\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[5], &A) != TCL_OK) {
    opserr << "WARNING invalid A\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[6], &E) != TCL_OK) {
    opserr << "WARNING invalid E\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[7], &I) != TCL_OK) {
    opserr << "WARNING invalid I\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[8], &ysID1) != TCL_OK) {
    opserr << "WARNING invalid ysID1\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[9], &ysID2) != TCL_OK) {
    opserr << "WARNING invalid ysID2\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[10], &cyc_type) != TCL_OK) {
    opserr << "WARNING invalid cyc_type\n";
    return TCL_ERROR;
  }

  /*if (Tcl_GetDouble(interp, argv[11], &wt) != TCL_OK)
  {
          opserr << "WARNING invalid wt\n";
          opserr << "element2dYS: " << tag << "\n";
          return TCL_ERROR;
  }*/

  double delpmax, alfa, beta;
  if (Tcl_GetDouble(interp, argv[11], &delpmax) != TCL_OK) {
    opserr << "WARNING invalid power\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[12], &alfa) != TCL_OK) {
    opserr << "WARNING invalid power\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[13], &beta) != TCL_OK) {
    opserr << "WARNING invalid rfalgo\n";
    return TCL_ERROR;
  }

  YieldSurface_BC *theYS1 = builder->getTypedObject<YieldSurface_BC>(ysID1);
  if (theYS1 == 0) {
    opserr << "WARNING element2dYS: " << tag << "\n";
    opserr << " no yield surface exists with tag: " << ysID1 << "\n";
    return TCL_ERROR;
  }

  YieldSurface_BC *theYS2 = builder->getTypedObject<YieldSurface_BC>(ysID2);
  if (theYS2 == 0) {
    opserr << "WARNING element2dYS: " << tag << "\n";
    opserr << " no yield surface exists with tag: " << ysID2 << "\n";
    return TCL_ERROR;
  }

  // Inelastic2DYS02(int tag, double a, double e, double i, int Nd1, int Nd2,
  //				YieldSurface_BC *ysEnd1,  YieldSurface_BC
  //*ysEnd2, 				int rf_algo, bool islinear, double rho)

  CyclicModel *theModel = builder->getTypedObject<CyclicModel>(cyc_type);
  // Element *theElement = new Inelastic2DYS02(tag, A, E, I, ndI, ndJ, theYS1,
  // theYS2, cyc_type, wt, delpmax, alfa, beta, rf_algo);
  Element *theElement =
      new Inelastic2DYS02(tag, A, E, I, ndI, ndJ, theYS1, theYS2, theModel,
                          delpmax, alfa, beta, rf_algo);
  opserr << "Inelastic2DYS02 created\n";

  if (tcl_debug)
    opserr << "\tElement created\n";

  // Ensure we have created the element, out of memory if got here and no
  // element
  if (theElement == 0) {
    opserr << "WARNING ran out of memory creating element\n";
    opserr << "\a";
    return TCL_ERROR;
  }

  opserr << "Inelastic2DYS02 adding to domain\n";

  if (builder->getDomain()->addElement(theElement) == false) {
    opserr << "WARNING TclElmtBuilder - addelement2dYS - could not add element "
              "to domain ";
    opserr << tag << "\n";
    opserr << "\a";
    return TCL_ERROR;
  }

  // if(tcl_debug)
  opserr << "Inelastic2DYS02 #" << tag << " added to domain - returning\n";

  return TCL_OK;
}

int
TclBasicBuilder_addElement2dYS03(ClientData clientData, Tcl_Interp *interp,
                                 int argc, TCL_Char ** const argv)
{
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  // cerr << "Press key to continue...\n";
  // cin.get();

  if (tcl_debug)
    opserr << " TclBasicBuilder_addElement2dGNL \n";

  if (argc < 11) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "element element2dYS03 tag? Nd1? Nd2? A_ten? A_com? E? IzPos? "
              "IzNeg? ysID1? ysID2? algo?";

    return TCL_ERROR;
  }

  int tag, ndI, ndJ;
  double E, aTens, aComp, Ipos, Ineg;
  //	double massDens = 0.0;
  int ysID1, ysID2;

  int rf_algo;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid element2dYS tag" << "\n";
    return TCL_ERROR;
  }
  if (tcl_debug)
    opserr << "\tElement tag = " << tag << "\n";

  if (Tcl_GetInt(interp, argv[3], &ndI) != TCL_OK) {
    opserr << "WARNING invalid node I\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[4], &ndJ) != TCL_OK) {
    opserr << "WARNING invalid node J\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[5], &aTens) != TCL_OK) {
    opserr << "WARNING invalid A\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[6], &aComp) != TCL_OK) {
    opserr << "WARNING invalid A\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[7], &E) != TCL_OK) {
    opserr << "WARNING invalid E\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[8], &Ipos) != TCL_OK) {
    opserr << "WARNING invalid I\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[9], &Ineg) != TCL_OK) {
    opserr << "WARNING invalid I\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[10], &ysID1) != TCL_OK) {
    opserr << "WARNING invalid ysID1\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[11], &ysID2) != TCL_OK) {
    opserr << "WARNING invalid ysID2\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[12], &rf_algo) != TCL_OK) {
    opserr << "WARNING invalid ysID1\n";
    return TCL_ERROR;
  }

  YieldSurface_BC *theYS1 = builder->getTypedObject<YieldSurface_BC>(ysID1);
  if (theYS1 == 0) {
    opserr << "WARNING element2dYS: " << tag << "\n";
    opserr << " no yield surface exists with tag: " << ysID1 << "\n";
    return TCL_ERROR;
  }

  YieldSurface_BC *theYS2 = builder->getTypedObject<YieldSurface_BC>(ysID2);
  if (theYS2 == 0) {
    opserr << "WARNING element2dYS: " << tag << "\n";
    opserr << " no yield surface exists with tag: " << ysID2 << "\n";
    return TCL_ERROR;
  }

  //	Inelastic2DYS03(int tag, double a_ten, double a_com, double e,
  //	                double iz_pos, double iz_neg, int Nd1, int Nd2,
  //                    YieldSurface_BC *ysEnd1,  YieldSurface_BC *ysEnd2,
  //                    int rf_algo, bool islinear, double rho);

  Element *theElement = new Inelastic2DYS03(tag, aTens, aComp, E, Ipos, Ineg,
                                            ndI, ndJ, theYS1, theYS2, rf_algo);

  opserr << "Inelastic2DYS03 created\n";

  if (tcl_debug)
    opserr << "\tElement created\n";

  // Ensure we have created the element, out of memory if got here and no
  // element
  if (theElement == 0) {
    opserr << "WARNING ran out of memory creating element\n";
    opserr << "\a";
    return TCL_ERROR;
  }

  opserr << "Inelastic2DYS03 adding to domain\n";

  if (builder->getDomain()->addElement(theElement) == false) {
    opserr << "WARNING TclElmtBuilder - addelement2dYS - could not add element "
              "to domain ";
    opserr << tag << "\n";
    opserr << "\a";
    return TCL_ERROR;
  }

  if (tcl_debug)
    opserr << "Inelastic2DYS03 #" << tag << " added to domain - returning\n";

  return TCL_OK;
}

#if 0
int
TclBasicBuilder_addElement2dYS04 (ClientData clientData, Tcl_Interp *interp,
                                                                   int argc,
                                  char **argv)
{
    BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

        if (argc < 11)
        {
                opserr << "WARNING insufficient arguments\n";
                opserr << "element element2dYS04 tag? Nd1? Nd2? A? E? Iz? ysID1? ysID2? algo?";
                return TCL_ERROR;
        }

        int tag, ndI, ndJ;
        double E, A, I;
//	double massDens = 0.0;
        int ysID1, ysID2;
        int rf_algo;

        if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
        {
                opserr << "WARNING invalid element2dYS04 tag" << "\n";
                return TCL_ERROR;
        }
    if(tcl_debug) opserr << "\tElement tag = " << tag << "\n";

    if (Tcl_GetInt(interp, argv[3], &ndI) != TCL_OK)
        {
                opserr << "WARNING invalid node I\n";
                opserr << "element2dYS04: " << tag << "\n";
                return TCL_ERROR;
        }

        if (Tcl_GetInt(interp, argv[4], &ndJ) != TCL_OK)
        {
                opserr << "WARNING invalid node J\n";
                opserr << "element2dYS04: " << tag << "\n";
                return TCL_ERROR;
        }


        if (Tcl_GetDouble(interp, argv[5], &A) != TCL_OK)
        {
                opserr << "WARNING invalid A\n";
                opserr << "element2dYS04: " << tag << "\n";
                return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[6], &E) != TCL_OK)
        {
                opserr << "WARNING invalid E\n";
                opserr << "element2dYS04: " << tag << "\n";
                return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[7], &I) != TCL_OK)
        {
                opserr << "WARNING invalid I\n";
                opserr << "element2dYS04: " << tag << "\n";
                return TCL_ERROR;
        }

        if (Tcl_GetInt(interp, argv[8], &ysID1) != TCL_OK)
        {
                opserr << "WARNING invalid ysID1\n";
                opserr << "element2dYS04: " << tag << "\n";
                return TCL_ERROR;
        }

        if (Tcl_GetInt(interp, argv[9], &ysID2) != TCL_OK)
        {
                opserr << "WARNING invalid ysID2\n";
                opserr << "element2dYS04: " << tag << "\n";
                return TCL_ERROR;
        }

        if (Tcl_GetInt(interp, argv[10], &rf_algo) != TCL_OK)
        {
                opserr << "WARNING invalid ysID1\n";
                opserr << "element2dYS04: " << tag << "\n";
                return TCL_ERROR;
        }

                YieldSurface_BC *theYS1 = builder->getTypedObject<YieldSurface_BC>(ysID1);
        if(theYS1 == 0)
        {
                opserr << "WARNING element2dYS: " << tag << "\n";
                opserr <<  " no yield surface exists with tag: " << ysID1 <<
"\n"; return TCL_ERROR;
        }

        YieldSurface_BC *theYS2 = builder->getTypedObject<YieldSurface_BC>(ysID2);
        if(theYS2 == 0)
        {
                opserr << "WARNING element2dYS: " << tag << "\n";
                opserr <<  " no yield surface exists with tag: " << ysID2 << "\n"; 
                return TCL_ERROR;
        }

//	Inelastic2DYS(	int tag, double A, double E, double I, int Nd1, int Nd2,
//					YieldSurface_BC *ysEnd1, YieldSurface_BC *ysEnd2,
//					int rf_algo = -1, // updated

        Element *theElement = new Inelastic2DYS04(tag, A, E, I, ndI, ndJ, theYS1, theYS2, rf_algo);

        if(tcl_debug) 
          opserr << "\tElement created\n";

        // Ensure we have created the element, out of memory if got here and no element 

        if (theElement == 0)
        {
                opserr << "WARNING ran out of memory creating element\n";
                opserr << "element2dYS04: " << tag << "\n";
                opserr << "\a";
                return TCL_ERROR;
        }

        if (theDomain->addElement(theElement) == false)
        {
                opserr << "WARNING TclElmtBuilder - addelement2dYS - could not add element to domain ";
                opserr << tag << "\n"; 
                opserr << "\a";
                return TCL_ERROR;
        }

        if(tcl_debug)
          opserr << "\tElement number " << tag << " added to domain - returning\n";
        return TCL_OK;
}
#endif

#if 0
int
TclBasicBuilder_addElement2dYS05 (ClientData clientData, Tcl_Interp *interp,
                                  int argc, char **argv)
{
        //cerr << "Press key to continue...\n";
        //cin.get();

        if (argc < 11)
        {
                opserr << "WARNING insufficient arguments\n";
                opserr << "element element2dYS04 tag? Nd1? Nd2? A? E? Iz? ysID1? ysID2? algo?";
                return TCL_ERROR;
        }

        int tag, ndI, ndJ;
        double E, A, I;
//	double massDens = 0.0;
        int ysID1, ysID2;
        int rf_algo;

        if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)
        {
                opserr << "WARNING invalid element2dYS05 tag" << "\n";
                return TCL_ERROR;
        }

        if(tcl_debug) opserr << "\tElement tag = " << tag << "\n";

        if (Tcl_GetInt(interp, argv[3], &ndI) != TCL_OK)
        {
                opserr << "WARNING invalid node I\n";
                opserr << "element2dYS05: " << tag << "\n";
                return TCL_ERROR;
        }

        if (Tcl_GetInt(interp, argv[4], &ndJ) != TCL_OK)
        {
                opserr << "WARNING invalid node J\n";
                opserr << "element2dYS05: " << tag << "\n";
                return TCL_ERROR;
        }


        if (Tcl_GetDouble(interp, argv[5], &A) != TCL_OK)
        {
                opserr << "WARNING invalid A\n";
                opserr << "element2dYS05: " << tag << "\n";
                return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[6], &E) != TCL_OK)
        {
                opserr << "WARNING invalid E\n";
                opserr << "element2dYS05: " << tag << "\n";
                return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[7], &I) != TCL_OK)
        {
                opserr << "WARNING invalid I\n";
                opserr << "element2dYS05: " << tag << "\n";
                return TCL_ERROR;
        }

        if (Tcl_GetInt(interp, argv[8], &ysID1) != TCL_OK)
        {
                opserr << "WARNING invalid ysID1\n";
                opserr << "element2dYS05: " << tag << "\n";
                return TCL_ERROR;
        }

        if (Tcl_GetInt(interp, argv[9], &ysID2) != TCL_OK)
        {
                opserr << "WARNING invalid ysID2\n";
                opserr << "element2dYS05: " << tag << "\n";
                return TCL_ERROR;
        }

        if (Tcl_GetInt(interp, argv[10], &rf_algo) != TCL_OK)
        {
                opserr << "WARNING invalid ysID1\n";
                opserr << "element2dYS05: " << tag << "\n";
                return TCL_ERROR;
        }

        YieldSurface_BC *theYS1 = builder->getTypedObject<YieldSurface_BC>(ysID1);
        if (theYS1 == nullptr)
            return TCL_ERROR;

        YieldSurface_BC *theYS2 = builder->getTypedObject<YieldSurface_BC>(ysID2);

        if (theYS2 == nullptr)
            return TCL_ERROR;


        Element *theElement = new Inelastic2DYS05(tag, A, E, I, ndI, ndJ, theYS1, theYS2, rf_algo);

        if(tcl_debug) opserr << "\tElement created\n";

        // Ensure we have created the element, out of memory if got here and no element 
        if (theElement == 0)
        {
                opserr << "WARNING ran out of memory creating element\n";
                opserr << "element2dYS05: " << tag << "\n";
                opserr << "\a";
                return TCL_ERROR;
        }

        if (theDomain->addElement(theElement) == false) {
          opserr << "WARNING TclElmtBuilder - addelement2dYS - could not add element to domain "; 
          opserr << tag << "\n"; opserr << "\a"; 
          return TCL_ERROR;
        }

        if(tcl_debug) 
          opserr << "\tElement number " << tag << " added to domain - returning\n";

        return TCL_OK;
}
#endif

/*******************************************************************************************/
int
TclBasicBuilder_addElement2dYS(ClientData clientData, Tcl_Interp *interp, int argc,
                               TCL_Char ** const argv)
{

  if (strcmp(argv[1], "inelastic2dYS01") == 0) {
    return TclBasicBuilder_addElement2dYS01(clientData, interp, argc, argv);

  } else if (strcmp(argv[1], "inelastic2dYS02") == 0) {
    return TclBasicBuilder_addElement2dYS02(clientData, interp, argc, argv);
  }

  else if (strcmp(argv[1], "inelastic2dYS03") == 0) {
    return TclBasicBuilder_addElement2dYS03(clientData, interp, argc, argv);
  }

  /*	else if (strcmp(argv[1],"inelastic2dYS04") == 0) {
            int result = TclBasicBuilder_addElement2dYS04
            (clientData, interp, argc, argv,
                                                   theTclDomain, theTclBuilder);
      return result;
    }*/
  /*else if (strcmp(argv[1],"inelastic2dYS05") == 0) {
          int result = TclBasicBuilder_addElement2dYS05
          (clientData, interp, argc, argv,
                                                 theTclDomain, theTclBuilder);
    return result;
  }*/
  else

    return TCL_ERROR;
}
