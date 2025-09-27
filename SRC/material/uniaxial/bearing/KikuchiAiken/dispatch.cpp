
#include <Logging.h>
#include <Parsing.h>
#include <ModelRegistry.h>
#include <string.h>
#include "KikuchiAikenHDR.h"
#include "KikuchiAikenLRB.h"

int
TclCommand_KikuchiAikenHDR(ClientData cd, Tcl_Interp* interp, int argc, G3_Char ** const argv)
{
  // uniaxialMaterial KikuchiAikenHDR matTag? tp? ar? hr? 
  //           <-coGHU cg? ch? cu?> <-coMSS rs? rf?>

  // arguments (necessary)
  int tag;
  int tp = -1;
  double ar;
  double hr;

  // arguments (optional)
  double cg = 1.0;
  double ch = 1.0;
  double cu = 1.0;
  double rs = 1.0;
  double rf = 1.0;

  //
  UniaxialMaterial *theMaterial = nullptr;

  // error flag
  bool ifNoError = true;

  if (argc < 6) { // uniaxialMaterial KikuchiAikenHDR matTag? tp? ar? hr?

    opserr << "WARNING invalid number of arguments\n";
    ifNoError = false;

  } else {

    // argv[2~5]
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid KikuchiAikenHDR tag" << OpenSees::SignalMessageEnd;
      ifNoError = false;
    }

    if ((strcmp(argv[3], "X0.6") == 0) || (strcmp(argv[3], "1") == 0)) {
      tp = 1;
    } else if ((strcmp(argv[3], "X0.6-0MPa") == 0) ||
               (strcmp(argv[3], "2") == 0)) {
      tp = 2;
    } else if ((strcmp(argv[3], "X0.4") == 0) || (strcmp(argv[3], "3") == 0)) {
      tp = 3;
    } else if ((strcmp(argv[3], "X0.4-0MPa") == 0) ||
               (strcmp(argv[3], "4") == 0)) {
      tp = 4;
    } else if ((strcmp(argv[3], "X0.3") == 0) || (strcmp(argv[3], "5") == 0)) {
      tp = 5;
    } else if ((strcmp(argv[3], "X0.3-0MPa") == 0) ||
               (strcmp(argv[3], "6") == 0)) {
      tp = 6;

    } else {
      opserr << "WARNING invalid KikuchiAikenHDR tp" << OpenSees::SignalMessageEnd;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[4], &ar) != TCL_OK || ar <= 0.0) {
      opserr << "WARNING invalid ar\n";
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[5], &hr) != TCL_OK || hr <= 0.0) {
      opserr << "WARNING invalid hr\n";
      ifNoError = false;
    }

    // argv[6~]
    for (int i = 6; i <= (argc - 1); i++) {

      if (strcmp(argv[i], "-coGHU") == 0 &&
          (i + 3) <= (argc - 1)) { // <-coGHU cg? ch? cu?>

        if (Tcl_GetDouble(interp, argv[i + 1], &cg) != TCL_OK || cg < 0.0) {
          opserr << "WARNING invalid cg\n";
          ifNoError = false;
        }

        if (Tcl_GetDouble(interp, argv[i + 2], &ch) != TCL_OK || ch < 0.0) {
          opserr << "WARNING invalid ch\n";
          ifNoError = false;
        }

        if (Tcl_GetDouble(interp, argv[i + 3], &cu) != TCL_OK || cu < 0.0) {
          opserr << "WARNING invalid cu" << OpenSees::SignalMessageEnd;
          ifNoError = false;
        }

        i += 3;

      } else if (strcmp(argv[i], "-coMSS") == 0 &&
                 (i + 2) <= (argc - 1)) { // <-coMSS rs? rf?>

        if (Tcl_GetDouble(interp, argv[i + 1], &rs) != TCL_OK || rs < 0.0) {
          opserr << "WARNING invalid rs" << OpenSees::SignalMessageEnd;
          ifNoError = false;
        }

        if (Tcl_GetDouble(interp, argv[i + 2], &rf) != TCL_OK || rf < 0.0) {
          opserr << "WARNING invalid rf" << OpenSees::SignalMessageEnd;
          ifNoError = false;
        }

        i += 2;

      } else { // invalid option
        opserr << "WARNING invalid optional arguments" << OpenSees::SignalMessageEnd;
        ifNoError = false;
        break;
      }
    }
  }

  // if error detected
  if (!ifNoError) {
    return TCL_ERROR;
  }

  // regard 0.0 input as mistake (substitute 1.0 for 0.0)
  if (cg == 0.0)
    cg = 1.0;
  if (ch == 0.0)
    ch = 1.0;
  if (cu == 0.0)
    cu = 1.0;
  if (rs == 0.0)
    rs = 1.0;
  if (rf == 0.0)
    rf = 1.0;

  // Parsing was successful, allocate the material
  theMaterial = new KikuchiAikenHDR(tag, tp, ar, hr, cg, ch, cu, rs, rf);

  // succeeded
  return ((ModelRegistry*)cd)->addTaggedObject<UniaxialMaterial>(*theMaterial);
}

int
TclCommand_KikuchiAikenLRB(ClientData cd, Tcl_Interp* interp, int argc, G3_Char ** const argv)
{

  // uniaxialMaterial KikuchiAikenLRB matTag? type? ar? hr? gr? "
  //        "ap? tp? alph? beta? <-T temp? > <-coKQ rk? rq?> <-coMSS rs? rf?>"

  // arguments (necessary)
  int tag;
  int type = 1;
  double ar = 0.0;
  double hr = 0.0;
  double gr = 0.392e6;
  double ap = 0.0;
  double tp = 8.33e6;
  double alph = 0.588e6;
  double beta = 13.0;

  // arguments (optional)
  double temp = 15.0;
  double rk = 1.0;
  double rq = 1.0;
  double rs = 1.0;
  double rf = 1.0;

  //
  UniaxialMaterial *theMaterial = 0;

  // error flag
  bool ifNoError = true;

  if (argc < 11) { // uniaxialMaterial KikuchiAikenLRB matTag? type? ar? hr? gr?
                   // ap? tp? alph? beta?

    opserr << "WARNING KikuchiAikenLRB invalid number of arguments\n";
    ifNoError = false;

  } else {

    // argv[2~10]
    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING KikuchiAikenLRB invalid tag" << OpenSees::SignalMessageEnd;
      ifNoError = false;
    }

    if (Tcl_GetInt(interp, argv[3], &type) != TCL_OK) {
      opserr << "WARNING KikuchiAikenLRB invalid type" << OpenSees::SignalMessageEnd;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[4], &ar) != TCL_OK || ar <= 0.0) {
      opserr << "WARNING KikuchiAikenLRB invalid ar" << OpenSees::SignalMessageEnd;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[5], &hr) != TCL_OK || hr <= 0.0) {
      opserr << "WARNING KikuchiAikenLRB invalid hr" << OpenSees::SignalMessageEnd;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[6], &gr) != TCL_OK || gr <= 0.0) {
      opserr << "WARNING KikuchiAikenLRB invalid gr" << OpenSees::SignalMessageEnd;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[7], &ap) != TCL_OK || ap <= 0.0) {
      opserr << "WARNING KikuchiAikenLRB invalid ap" << OpenSees::SignalMessageEnd;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[8], &tp) != TCL_OK || tp <= 0.0) {
      opserr << "WARNING KikuchiAikenLRB invalid tp" << OpenSees::SignalMessageEnd;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[9], &alph) != TCL_OK || alph <= 0.0) {
      opserr << "WARNING KikuchiAikenLRB invalid alph" << OpenSees::SignalMessageEnd;
      ifNoError = false;
    }

    if (Tcl_GetDouble(interp, argv[10], &beta) != TCL_OK || beta <= 0.0) {
      opserr << "WARNING KikuchiAikenLRB invalid beta" << OpenSees::SignalMessageEnd;
      ifNoError = false;
    }

    // argv[11~]
    for (int i = 11; i <= (argc - 1); i++) {

      if (strcmp(argv[i], "-T") == 0 && (i + 1) <= (argc - 1)) { // <-T temp?>

        if (Tcl_GetDouble(interp, argv[i + 1], &temp) != TCL_OK) {
          opserr << "WARNING KikuchiAikenLRB invalid temp" << OpenSees::SignalMessageEnd;
          ifNoError = false;
        }

        i += 1;

      } else if (strcmp(argv[i], "-coKQ") == 0 &&
                 (i + 2) <= (argc - 1)) { // <-coKQ rk? rq?>

        if (Tcl_GetDouble(interp, argv[i + 1], &rk) != TCL_OK || rk < 0.0) {
          opserr << "WARNING KikuchiAikenLRB invalid rk" << OpenSees::SignalMessageEnd;
          ifNoError = false;
        }

        if (Tcl_GetDouble(interp, argv[i + 2], &rq) != TCL_OK || rq < 0.0) {
          opserr << "WARNING KikuchiAikenLRB invalid rq" << OpenSees::SignalMessageEnd;
          ifNoError = false;
        }

        i += 2;
      } else if (strcmp(argv[i], "-coMSS") == 0 &&
                 (i + 2) <= (argc - 1)) { // <-coMSS rs? rf?>

        if (Tcl_GetDouble(interp, argv[i + 1], &rs) != TCL_OK || rs < 0.0) {
          opserr << "WARNING KikuchiAikenLRB invalid rs" << OpenSees::SignalMessageEnd;
          ifNoError = false;
        }

        if (Tcl_GetDouble(interp, argv[i + 2], &rf) != TCL_OK || rf < 0.0) {
          opserr << "WARNING KikuchiAikenLRB invalid rf" << OpenSees::SignalMessageEnd;
          ifNoError = false;
        }

        i += 2;

      } else { // invalid option
        opserr << "WARNING KikuchiAikenLRB invalid optional arguments" << OpenSees::SignalMessageEnd;
        ifNoError = false;
        break;
      }
    }
  }

  // if error detected
  if (!ifNoError) {
    return TCL_ERROR;
  }

  // regard 0.0 input as misteke (substitute 1.0 for 0.0)
  if (rk == 0.0)
    rk = 1.0;
  if (rq == 0.0)
    rq = 1.0;
  if (rs == 0.0)
    rs = 1.0;
  if (rf == 0.0)
    rf = 1.0;

  // Parsing was successful, allocate the material
  theMaterial = new KikuchiAikenLRB(tag, type, ar, hr, gr, ap, tp, alph, beta,
                                    temp, rk, rq, rs, rf);

  return ((ModelRegistry*)cd)->addTaggedObject<UniaxialMaterial>(*theMaterial);
}

