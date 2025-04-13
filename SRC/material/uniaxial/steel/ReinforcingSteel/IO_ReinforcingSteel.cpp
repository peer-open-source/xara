

#include <g3_api.h>


#include <SRC/material/uniaxial/ReinforcingSteel.h>
void *OPS_ReinforcingSteel(G3_Runtime *rt)
{
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 7) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "uniaxialMaterial ReinforcingSteel ";
    opserr << "tag? fy? fu? Es? Esh? esh? eult? ";
    opserr << "<-GABuck?> <-DMBuck?> <-CMFatigue?> <-MPCurveParams?> "
              "<-IsoHard?>\n";
    return 0;
  }

  int tag;
  numdata = 1;
  if (OPS_GetIntInput(&numdata, &tag) < 0) {
    opserr << "WARNING invalid tag\n";
    return 0;
  }

  double data[6];
  numdata = 6;
  if (OPS_GetDoubleInput(&numdata, data) < 0) {
    opserr << "WARNING invalid double data\n";
    return 0;
  }

  int buckModel = 0;
  double gabuckdata[4] = {0.0, 1.0, 1.0, 0.5};
  double dmbuckdata[2] = {0.0, 1.0};
  double fatiguedata[3] = {0.0, -4.46, 0.0};
  double mpdata[3] = {1.0 / 3.0, 18.0, 4.0};
  double isodata[2] = {0.0, 0.01};
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *type = OPS_GetString();

    if (strcmp(type, "-GABuck") == 0) {
      numdata = OPS_GetNumRemainingInputArgs();
      if (numdata < 4) {
        opserr << "WARNING insufficient optional arguments for -GABuck\n";
        opserr << "Want: <-GABuck lsr? beta? r? gama?>\n";
        return 0;
      }
      buckModel = 1;
      numdata = 4;
      if (OPS_GetDoubleInput(&numdata, gabuckdata) < 0) {
        opserr << "WARNING invalid double data\n";
        return 0;
      }

    } else if (strcmp(type, "-DMBuck") == 0) {
      numdata = OPS_GetNumRemainingInputArgs();
      if (numdata < 2) {
        opserr << "WARNING insufficient optional arguments for -DMBuck\n";
        opserr << "Want: <-DMBuck lsr? alpha?>\n";
        return 0;
      }
      buckModel = 2;
      numdata = 2;

      if (OPS_GetDoubleInput(&numdata, dmbuckdata) < 0) {
        opserr << "WARNING invalid double data\n";
        return 0;
      }
      if (dmbuckdata[1] < 0.75 || dmbuckdata[1] > 1.0) {
        opserr << "WARNING alpha usually is between 0.75 and 1.0\n";
        return 0;
      }

    } else if (strcmp(type, "-CMFatigue") == 0) {
      numdata = OPS_GetNumRemainingInputArgs();
      if (numdata < 3) {
        opserr << "WARNING insufficient optional arguments for -CMFatigue\n";
        opserr << "Want: <-CMFatigue Cf? alpha? Cd?>\n";
        return 0;
      }
      numdata = 3;
      if (OPS_GetDoubleInput(&numdata, fatiguedata) < 0) {
        opserr << "WARNING invalid double data\n";
        return 0;
      }

    } else if (strcmp(type, "-MPCurveParams") == 0) {
      numdata = OPS_GetNumRemainingInputArgs();
      if (numdata < 3) {
        opserr
            << "WARNING insufficient optional arguments for -MPCurveParams\n";
        opserr << "Want: <-CMFatigue R1? R2? R3?>\n";
        return 0;
      }
      numdata = 3;
      if (OPS_GetDoubleInput(&numdata, mpdata)) {
        opserr << "WARNING invalid double data\n";
        return 0;
      }
    } else if (strcmp(type, "-IsoHard") == 0) {
      numdata = OPS_GetNumRemainingInputArgs();
      if (numdata < 2) {
        opserr << "WARNING insufficient optional arguments for -IsoHard\n";
        opserr << "Want: <-IsoHard a1 limit>\n";
        return 0;
      }
      numdata = 2;
      if (OPS_GetDoubleInput(&numdata, isodata)) {
        opserr << "WARNING invalid double data\n";
        return 0;
      }
    } else {
      opserr << "WARNING did not recognize optional flag\n";
      opserr << "Possible Optional Flags: ";
      opserr << "<-GABuck?> <-DMBuck?> <-CMFatigue?> <-MPCurveParams?> "
                "<-IsoHard?>\n";
      return 0;
    }
  }

  // Parsing was successful, allocate the material
  double slen = 0.0, beta = 1.0;
  if (buckModel == 1) {
    slen = gabuckdata[0];
    beta = gabuckdata[1];
  } else if (buckModel == 2) {
    slen = dmbuckdata[0];
    beta = dmbuckdata[1];
  }
  UniaxialMaterial *theMaterial = 0;
  theMaterial = new ReinforcingSteel(
      tag, data[0], data[1], data[2], data[3], data[4], data[5], buckModel,
      slen, beta, gabuckdata[2], gabuckdata[3], fatiguedata[0], fatiguedata[1],
      fatiguedata[2], mpdata[0], mpdata[1], mpdata[2], isodata[0], isodata[1]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type "
              "ReinforcingSteel\n";
    return 0;
  }

  return theMaterial;
}
