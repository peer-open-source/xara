//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
#include <tcl.h>
#include <assert.h>
#include <Logging.h>
#include <Parsing.h>
#include <elementAPI.h>
#include <ArgumentTracker.h>


#include <ID.h>
#include <FrameSection.h>
#include <BasicModelBuilder.h>
#include <BeamIntegration.h>
#include <FrameQuadrature.hpp>
#include <quadrature/GaussLobatto1D.hpp>
#include <LobattoBeamIntegration.h>
#include <LegendreBeamIntegration.h>
#include <RadauBeamIntegration.h>
#include <NewtonCotesBeamIntegration.h>
#include <UserDefinedBeamIntegration.h>

#include <HingeMidpointBeamIntegration.h>
#include <HingeEndpointBeamIntegration.h>
#include <HingeRadauBeamIntegration.h>
#include <HingeRadauTwoBeamIntegration.h>
#include <UserDefinedHingeIntegration.h>
#include <DistHingeIntegration.h>
#include <RegularizedHingeIntegration.h>

#include <TrapezoidalBeamIntegration.h>
#include <CompositeSimpsonBeamIntegration.h>
#include <FixedLocationBeamIntegration.h>
#include <LowOrderBeamIntegration.h>
#include <MidDistanceBeamIntegration.h>
//#include <GaussQBeamIntegration.h>

extern "C" int OPS_ResetInputNoBuilder(ClientData clientData, Tcl_Interp *interp, int cArg,
                          int mArg, TCL_Char ** const argv, Domain *domain);

extern void *OPS_LobattoBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_LegendreBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_NewtonCotesBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_RadauBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_TrapezoidalBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_CompositeSimpsonBeamIntegration(int &integrationTag,
                                                 ID &secTags);
extern void *OPS_UserDefinedBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_FixedLocationBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_LowOrderBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_MidDistanceBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_UserHingeBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_HingeMidpointBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_HingeRadauBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_HingeRadauTwoBeamIntegration(int &integrationTag, ID &secTags);
extern void *OPS_HingeEndpointBeamIntegration(int &integrationTag, ID &secTags);
extern void* OPS_ConcentratedPlasticityBeamIntegration(int&, ID&);
extern void* OPS_ConcentratedCurvatureBeamIntegration(int&, ID&);


BeamIntegration*
GetBeamIntegration(TCL_Char* type)
{
  if (strcmp(type, "Lobatto") == 0)
    return new FrameQuadrature<GaussLobatto<1,20>>; // LobattoBeamIntegration();

  else if (strcmp(type, "Legendre") == 0)
    return new LegendreBeamIntegration();

  else if (strcmp(type, "Radau") == 0)
    return new RadauBeamIntegration();

  else if (strcmp(type, "NewtonCotes") == 0)
    return new NewtonCotesBeamIntegration();

  else if (strcmp(type, "Trapezoidal") == 0)
    return new TrapezoidalBeamIntegration();

  else if (strcmp(type, "CompositeSimpson") == 0)
    return new CompositeSimpsonBeamIntegration();
  else
    return nullptr;
}

extern int
TclCommand_addBeamIntegration(ClientData clientData, Tcl_Interp *interp,
                              int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  if (argc < 2) {
    opserr << G3_ERROR_PROMPT << "want beamIntegration type tag...\n";
    return TCL_ERROR;
  }

  OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, nullptr);

  int iTag;
  ID secTags;
  BeamIntegration *bi = nullptr;
  if (strcmp(argv[1], "Lobatto") == 0) {
    bi = (BeamIntegration *)OPS_LobattoBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "Legendre") == 0) {
    bi = (BeamIntegration *)OPS_LegendreBeamIntegration(iTag, secTags);

  } else if (strcmp(argv[1], "NewtonCotes") == 0) {
    bi = (BeamIntegration *)OPS_NewtonCotesBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "Radau") == 0) {
    bi = (BeamIntegration *)OPS_RadauBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "Trapezoidal") == 0) {
    bi = (BeamIntegration *)OPS_TrapezoidalBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "CompositeSimpson") == 0) {
    bi = (BeamIntegration *)OPS_CompositeSimpsonBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "UserDefined") == 0) {
    bi = (BeamIntegration *)OPS_UserDefinedBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "FixedLocation") == 0) {
    bi = (BeamIntegration *)OPS_FixedLocationBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "LowOrder") == 0) {
    bi = (BeamIntegration *)OPS_LowOrderBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "MidDistance") == 0) {
    bi = (BeamIntegration *)OPS_MidDistanceBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "UserHinge") == 0) {
    bi = (BeamIntegration *)OPS_UserHingeBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "HingeMidpoint") == 0) {
    bi = (BeamIntegration *)OPS_HingeMidpointBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "HingeRadau") == 0) {
    bi = (BeamIntegration *)OPS_HingeRadauBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "HingeRadauTwo") == 0) {
    bi = (BeamIntegration *)OPS_HingeRadauTwoBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1], "HingeEndpoint") == 0) {
    bi = (BeamIntegration *)OPS_HingeEndpointBeamIntegration(iTag, secTags);
  } else if (strcmp(argv[1],"ConcentratedPlasticity") == 0) {
    bi = (BeamIntegration*)OPS_ConcentratedPlasticityBeamIntegration(iTag,secTags);
  } else if (strcmp(argv[1],"ConcentratedCurvature") == 0) {
    bi = (BeamIntegration*)OPS_ConcentratedCurvatureBeamIntegration(iTag,secTags);
  } else {
    opserr << "WARNING: integration type " << argv[1] << " is unknown\n";
    return TCL_ERROR;
  }

  assert(bi);
  BeamIntegrationRule *rule = new BeamIntegrationRule(iTag, bi, secTags);


  if (builder->addTypedObject<BeamIntegrationRule>(iTag, rule) < 0) {
    opserr << G3_ERROR_PROMPT << "could not add BeamIntegrationRule.";
    delete rule;
    return TCL_ERROR;
  }

  return TCL_OK;
}


int
TclCommand_CreateHingeStencil(ClientData clientData, Tcl_Interp *interp,
                              int argc, TCL_Char **const argv)
{

  bool parse_tag = false;
  const char* const type = argv[6];

  int nIP;
  std::vector<int> sections;
  BeamIntegration *beamIntegr = nullptr;

  if (strcmp(type, "Lobatto") == 0 || 
      strcmp(type, "Legendre") == 0 ||
      strcmp(type, "Radau") == 0 ||
      strcmp(type, "NewtonCotes") == 0 ||
      strcmp(type, "Trapezoidal") == 0 ||
      strcmp(type, "CompositeSimpson") == 0) {
    // 0      1     2
    // $Type tag $secTag $nIP
    enum class Argument {Tag, Section, nIP, End};
    ArgumentTracker<Argument> tracker;

    if (!parse_tag)
      tracker.consume(Argument::Tag);

    if (argc < 3) { // 3 = 9 - 6
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      opserr << "Want: "
             //           0       1
             << " Lobatto secTag? nIP?\n";
      return TCL_ERROR;
    }

    int secTag;
    if (Tcl_GetInt(interp, argv[0], &secTag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid secTag\n";
      return TCL_ERROR;
    }
    int nIP;
    if (Tcl_GetInt(interp, argv[1], &nIP) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid nIP\n";
      return TCL_ERROR;
    }

    for (int i = 0; i < nIP; i++)
      sections.push_back(secTag);

    if (strcmp(type, "Lobatto") == 0)
      beamIntegr = new LobattoBeamIntegration();
    else if (strcmp(type, "Legendre") == 0)
      beamIntegr = new LegendreBeamIntegration();
    else if (strcmp(type, "Radau") == 0)
      beamIntegr = new RadauBeamIntegration();
    else if (strcmp(type, "NewtonCotes") == 0)
      beamIntegr = new NewtonCotesBeamIntegration();
    else if (strcmp(type, "Trapezoidal") == 0)
      beamIntegr = new TrapezoidalBeamIntegration();
    else if (strcmp(type, "CompositeSimpson") == 0)
      beamIntegr = new CompositeSimpsonBeamIntegration();
    else {
      opserr << "ERROR: invalid integration type: " << type << "\n";
      return TCL_ERROR;
    }
  }

  else if (strcmp(type, "UserDefined") == 0) {

    if (argc < 3) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      opserr << "Want: "
             //              0    1
             << "UserDefined nIP? secTag1? "
                "... pt1? ... wt1? ...\n";
      return TCL_ERROR;
    }

    int nIP;
    if (Tcl_GetInt(interp, argv[0], &nIP) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid nIP\n";
      return TCL_ERROR;
    }

    sections.resize(nIP, -1);
    Vector pts(nIP);
    Vector wts(nIP);

    int i, j;
    for (i = 0, j = 1; i < nIP; i++, j++) {
      int sec;
      double pt, wt;
      if (Tcl_GetInt(interp, argv[j], &sec) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid sec\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j + nIP], &pt) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid pt\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j + 2 * nIP], &wt) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid wt\n";
        return TCL_ERROR;
      }
      sections[i] = sec;
      pts[i] = pt;
      wts[i] = wt;
    }

    beamIntegr = new UserDefinedBeamIntegration(nIP, pts, wts);
  }

  //
  //
  //
  else if (strcmp(type, "HingeMidpoint") == 0 ||
           strcmp(type, "HingeRadau") == 0 ||
           strcmp(type, "HingeRadauTwo") == 0 ||
           strcmp(type, "HingeEndpoint") == 0) {

    if (argc < 6) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      opserr << "Want: "
             //        0        1    2        3    4
             << " type secTagI? lpI? secTagJ? lpJ? secTagE?\n";
      return TCL_ERROR;
    }

    double lpI, lpJ;
    int secTagI, secTagJ, secTagE;

    if (Tcl_GetInt(interp, argv[0], &secTagI) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid secTagI\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[1], &lpI) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid lpI\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[2], &secTagJ) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid secTagJ\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[3], &lpJ) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid lpJ\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[4], &secTagE) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid secTagE\n";
      return TCL_ERROR;
    }

    int nIP;

    if (strcmp(type, "HingeMidpoint") == 0) {
      beamIntegr = new HingeMidpointBeamIntegration(lpI, lpJ);
      nIP = 4;
      sections.resize(nIP, -1);
      sections[0] = secTagI;
      sections[1] = secTagE;
      sections[2] = secTagE;
      sections[3] = secTagJ;

    } else if (strcmp(type, "HingeRadau") == 0) {
      beamIntegr = new HingeRadauBeamIntegration(lpI, lpJ);
      nIP = 6;
      sections.resize(nIP, -1);
      sections[0] = secTagI;
      sections[1] = secTagE;
      sections[2] = secTagE;
      sections[3] = secTagE;
      sections[4] = secTagE;
      sections[5] = secTagJ;

    } else if (strcmp(type, "HingeRadauTwo") == 0) {
      beamIntegr = new HingeRadauTwoBeamIntegration(lpI, lpJ);
      nIP = 6;
      sections.resize(nIP, -1);
      sections[0] = secTagI;
      sections[1] = secTagI;
      sections[2] = secTagE;
      sections[3] = secTagE;
      sections[4] = secTagJ;
      sections[5] = secTagJ;

    } else {
      beamIntegr = new HingeEndpointBeamIntegration(lpI, lpJ);
      nIP = 4;
      sections.resize(nIP, -1);
      sections[0] = secTagI;
      sections[1] = secTagE;
      sections[2] = secTagE;
      sections[3] = secTagJ;
    }
  }

  else if (strcmp(type, "UserHinge") == 0) {

    if (argc < 3) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      opserr << "Want: "
             //             0        1
             << " UserHinge secTagE? npL? "
                "secTagL1? ... ptL1? ... wtL1? ... npR? secTagR1? ... ptR1? "
                "... wtR1? ...\n";
      return TCL_ERROR;
    }

    int secTagE;
    if (Tcl_GetInt(interp, argv[0], &secTagE) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid secTagE\n";
      return TCL_ERROR;
    }

    int argStart = 1;

    int npL, npR;

    if (Tcl_GetInt(interp, argv[argStart], &npL) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid npL\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[argStart + 3 * npL + 1], &npR) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid npR\n";
      return TCL_ERROR;
    }

    int nIP = npL + npR;

    sections.resize(nIP + 2, -1);
    Vector ptsL(npL);
    Vector wtsL(npL);
    Vector ptsR(npR);
    Vector wtsR(npR);

    int i, j;
    for (i = 0, j = argStart + 1; i < npL; i++, j++) {
      int sec;
      double pt, wt;
      if (Tcl_GetInt(interp, argv[j], &sec) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid sec\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j + npL], &pt) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid pt\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j + 2 * npL], &wt) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid wt\n";
        return TCL_ERROR;
      }
      sections[i+1] = sec;
      ptsL[i] = pt;
      wtsL[i] = wt;
    }

    for (i = 0, j = 1 + (argStart + 1) + 3 * npL; i < npR; i++, j++) {
      int sec;
      double pt, wt;
      if (Tcl_GetInt(interp, argv[j], &sec) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid sec\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j + npR], &pt) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid pt\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j + 2 * npR], &wt) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid wt\n";
        return TCL_ERROR;
      }
      sections[i + npL] = sec;
      ptsR[i] = pt;
      wtsR[i] = wt;
    }

    sections[nIP] = secTagE;
    sections[nIP + 1] = secTagE;

    beamIntegr =
        new UserDefinedHingeIntegration(npL, ptsL, wtsL, npR, ptsR, wtsR);

    nIP += 2;
  }


  else if (strcmp(type, "DistHinge") == 0) {

    if (argc < 8) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      opserr << "Want: "
             //        0        1    2
             << " type distType nIP? secTagI? "
                "lpI? secTagJ? lpJ? secTagE?\n";
      return TCL_ERROR;
    }

    BeamIntegration *otherBeamInt = nullptr;
    if (strcmp(argv[0], "Lobatto") == 0)
      otherBeamInt = new LobattoBeamIntegration();
    else if (strcmp(argv[0], "Legendre") == 0)
      otherBeamInt = new LegendreBeamIntegration();
    else if (strcmp(argv[0], "Radau") == 0)
      otherBeamInt = new RadauBeamIntegration();
    else if (strcmp(argv[0], "NewtonCotes") == 0)
      otherBeamInt = new NewtonCotesBeamIntegration();
    else if (strcmp(argv[0], "Trapezoidal") == 0)
      otherBeamInt = new TrapezoidalBeamIntegration();
    else if (strcmp(argv[0], "CompositeSimpson") == 0)
      otherBeamInt = new CompositeSimpsonBeamIntegration();
    else {
      opserr << "ERROR: invalid integration type: " << argv[0] << "\n";
      return TCL_ERROR;
    }


    int nIP;
    if (Tcl_GetInt(interp, argv[1], &nIP) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid nIP\n";
      return TCL_ERROR;
    }
    int secTagI, secTagJ, secTagE;
    if (Tcl_GetInt(interp, argv[2], &secTagI) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid secTagI\n";
      return TCL_ERROR;
    }
    double lpI, lpJ;
    if (Tcl_GetDouble(interp, argv[3], &lpI) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid lpI\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[4], &secTagJ) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid secTagJ\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[5], &lpJ) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid lpJ\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[6], &secTagE) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid secTagE\n";
      return TCL_ERROR;
    }


    nIP = 2 * nIP;
//  sections = new FrameSection *[nIP + 2];
    for (int i = 0; i < nIP; i++) {
      sections[i] = secTagI;
      sections[i + nIP] = secTagJ;
    }

    sections[nIP] = secTagE;
    sections[nIP + 1] = secTagE;

    beamIntegr = new DistHingeIntegration(lpI, lpJ, *otherBeamInt);

    nIP += 2;

    if (otherBeamInt != 0)
      delete otherBeamInt;
  }

  //
  //
  //
  else if (strcmp(type, "RegularizedHinge") == 0) {

    if (argc < 10) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      opserr << "Want: "
             //        0        1    2 
             << " type distType nIP? secTagI? "
                "lpI? zetaI? secTagJ? lpJ? zetaJ? secTagE?\n";
      return TCL_ERROR;
    }

    int secTagI, secTagJ, secTagE;
    double lpI, lpJ;
    double zetaI, zetaJ;
    int nIP;

    BeamIntegration *otherBeamInt = nullptr;
    if (strcmp(argv[0], "Lobatto") == 0)
      otherBeamInt = new LobattoBeamIntegration();
    else if (strcmp(argv[0], "Legendre") == 0)
      otherBeamInt = new LegendreBeamIntegration();
    else if (strcmp(argv[0], "Radau") == 0)
      otherBeamInt = new RadauBeamIntegration();
    else if (strcmp(argv[0], "NewtonCotes") == 0)
      otherBeamInt = new NewtonCotesBeamIntegration();
    else if (strcmp(argv[0], "Trapezoidal") == 0)
      otherBeamInt = new TrapezoidalBeamIntegration();
    else if (strcmp(argv[0], "CompositeSimpson") == 0)
      otherBeamInt = new CompositeSimpsonBeamIntegration();
    else {
      opserr << "ERROR: invalid integration type: " << argv[0] << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[1], &nIP) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid nIP\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[2], &secTagI) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid secTagI\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[3], &lpI) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid lpI\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[4], &zetaI) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid zetaI\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[5], &secTagJ) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid secTagJ\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[6], &lpJ) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid lpJ\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[7], &zetaJ) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid zetaI\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[8], &secTagE) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid secTagE\n";
      return TCL_ERROR;
    }

    nIP = nIP;
    sections.resize(nIP+2, -1);
    for (int i = 1; i < nIP - 1; i++) {
      sections[i] = secTagE;
    }

    sections[0]       = secTagI;
    sections[nIP]     = secTagI;
    sections[nIP - 1] = secTagJ;
    sections[nIP + 1] = secTagJ;

    beamIntegr =
        new RegularizedHingeIntegration(*otherBeamInt, lpI, lpJ, zetaI, zetaJ);

    nIP += 2;

    if (otherBeamInt != nullptr)
      delete otherBeamInt;
  }

  //
  //
  //
  else if (strcmp(type, "FixedLocation") == 0) {

    if (argc < 3) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      opserr << "Want: "
             //                 0    1
             << " FixedLocation nIP? secTag1? "
                "... pt1? ... \n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[0], &nIP) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid nIP\n";
      return TCL_ERROR;
    }

    sections.resize(nIP, -1);
    Vector pts(nIP);
    Vector wts(nIP);

    int i, j;
    for (i = 0, j = 1; i < nIP; i++, j++) {
      int sec;
      double pt;
      if (Tcl_GetInt(interp, argv[j], &sec) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid sec\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j + nIP], &pt) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid pt\n";
        return TCL_ERROR;
      }
      sections[i] = sec;
      pts[i]      = pt;
    }

    beamIntegr = new FixedLocationBeamIntegration(nIP, pts);
  }

  //
  //
  //
  else if (strcmp(type, "LowOrder") == 0) {

    if (argc < 3) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      opserr << "Want: "
             //            0    1
             << " LowOrder nIP? secTag1? ... "
                "pt1? ... wt1? ...\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[0], &nIP) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid nIP\n";
      return TCL_ERROR;
    }

    sections.resize(nIP, -1);
    Vector pts(nIP);
    Vector wts(nIP);

    int i, j;
    int nc = 0;
    for (i = 0, j = 1; i < nIP; i++, j++) {
      int sec;
      double pt, wt;
      if (Tcl_GetInt(interp, argv[j], &sec) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid sec\n";
        return TCL_ERROR;
      }
      sections[i] = sec;

      if (Tcl_GetDouble(interp, argv[j + nIP], &pt) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid pt\n";
        return TCL_ERROR;
      }
      pts[i] = pt;

      if (j + 2 * nIP < argc) {
        if (Tcl_GetDouble(interp, argv[j + 2 * nIP], &wt) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid wt\n";
          return TCL_ERROR;
        } else {
          wts(i) = wt;
          nc++;
        }
      }
    }

    beamIntegr = new LowOrderBeamIntegration(nIP, pts, nc, wts);
  }

  //
  //
  //
  else if (strcmp(type, "MidDistance") == 0) {

    if (argc < 3) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      opserr << "Want: "
             //               0    1
             << " MidDistance nIP? secTag1? "
                "... pt1? ... \n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[0], &nIP) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid nIP\n";
      return TCL_ERROR;
    }

    sections.resize(nIP, -1);
    Vector pts(nIP);
    Vector wts(nIP);

    int i, j;
    for (i = 0, j = 1; i < nIP; i++, j++) {
      int sec;
      double pt;
      if (Tcl_GetInt(interp, argv[j], &sec) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid sec\n";
        return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[j + nIP], &pt) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid pt\n";
        return TCL_ERROR;
      }
      sections[i] = sec;
      pts[i] = pt;
    }


    // Sort locations before calling constructor
    for (int i = 0; i < nIP; i++) {
      int key = i;
      for (int j = i + 1; j < nIP; j++) {
        if (pts(j) < pts(key)) {
          key = j;
        }
      }
      if (key != i) {
        // Swap locs
        double temp;
        temp = pts[i];
        pts[i] = pts(key);
        pts(key) = temp;

        // Swap sections
        int tempTag = sections[i];
        sections[i] = sections[key];
        sections[key] = tempTag;
      }
    }

    beamIntegr = new MidDistanceBeamIntegration(nIP, pts);
  }

// TODO
#if 0
  else if (strcmp(type,"GaussQ") == 0) {

    int type, secTag;

    if (argc < 10) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
      opserr << "Want: " << " GaussQ type? secTag? nIP?\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[0], &type) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid type\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[8], &secTag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid secTag\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[9], &nIP) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid nIP\n";
      return TCL_ERROR;
    }

    sections.resize(nIP, -1);
    for (int i = 0; i < nIP; i++)
      sections[i] = theSection;

    beamIntegr = new GaussQBeamIntegration(type);
  }
#endif

  else {
    opserr << "Unknown integration type: " << type << "\n";
    return TCL_ERROR;
  }
}
