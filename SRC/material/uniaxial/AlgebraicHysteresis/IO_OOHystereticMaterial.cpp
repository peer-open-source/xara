

#include <g3_api.h>


#include <SRC/material/uniaxial/OOHystereticMaterial.h>
void *OPS_OOHystereticMaterial(void)
{
  UniaxialMaterial *theMaterial = 0;

  if (OPS_GetNumRemainingInputArgs() < 5) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial OOHysteretic tag? bTag+? unlRulTag+? "
              "stfDegTag+? strDegTag+? "
           << "<bTag-? unlRulTag-? stfDegTag-? strDegTag-?> <pinchX? pinchY?>"
           << endln;
    return 0;
  }

  int tag;
  int bTagPos, bTagNeg;
  int unlTagPos, unlTagNeg;
  int stfTagPos, stfTagNeg;
  int strTagPos, strTagNeg;
  double pinchX = 0.0;
  double pinchY = 1.0;

  int argc = OPS_GetNumRemainingInputArgs();

  int numData = 1;
  if (OPS_GetIntInput(&numData, &tag) != 0) {
    opserr << "WARNING invalid tag\n";
    opserr << "OOHysteretic material: " << tag << endln;
    return 0;
  }
  if (OPS_GetIntInput(&numData, &bTagPos) != 0) {
    opserr << "WARNING invalid bTag+\n";
    opserr << "OOHysteretic material: " << tag << endln;
    return 0;
  }
  if (OPS_GetIntInput(&numData, &unlTagPos) != 0) {
    opserr << "WARNING invalid unlRulTag+\n";
    opserr << "OOHysteretic material: " << tag << endln;
    return 0;
  }
  if (OPS_GetIntInput(&numData, &stfTagPos) != 0) {
    opserr << "WARNING invalid stfDegTag+\n";
    opserr << "OOHysteretic material: " << tag << endln;
    return 0;
  }
  if (OPS_GetIntInput(&numData, &strTagPos) != 0) {
    opserr << "WARNING invalid strDegTag+\n";
    opserr << "OOHysteretic material: " << tag << endln;
    return 0;
  }

  if (argc == 7) {
    if (OPS_GetDoubleInput(&numData, &pinchX) != 0) {
      opserr << "WARNING invalid pinchX\n";
      opserr << "OOHysteretic material: " << tag << endln;
      return 0;
    }
    if (OPS_GetDoubleInput(&numData, &pinchY) != 0) {
      opserr << "WARNING invalid pinchY\n";
      opserr << "OOHysteretic material: " << tag << endln;
      return 0;
    }
  }

  if (argc > 8) {
    if (OPS_GetIntInput(&numData, &bTagNeg) != 0) {
      opserr << "WARNING invalid bTag-\n";
      opserr << "OOHysteretic material: " << tag << endln;
      return 0;
    }
    if (OPS_GetIntInput(&numData, &unlTagNeg) != 0) {
      opserr << "WARNING invalid unlRulTag-\n";
      opserr << "OOHysteretic material: " << tag << endln;
      return 0;
    }
    if (OPS_GetIntInput(&numData, &stfTagNeg) != 0) {
      opserr << "WARNING invalid stfDegTag-\n";
      opserr << "OOHysteretic material: " << tag << endln;
      return 0;
    }
    if (OPS_GetIntInput(&numData, &strTagNeg) != 0) {
      opserr << "WARNING invalid strDegTag-\n";
      opserr << "OOHysteretic material: " << tag << endln;
      return 0;
    }
  }
  if (argc == 11) {
    if (OPS_GetDoubleInput(&numData, &pinchX) != 0) {
      opserr << "WARNING invalid pinchX\n";
      opserr << "OOHysteretic material: " << tag << endln;
      return 0;
    }
    if (OPS_GetDoubleInput(&numData, &pinchY) != 0) {
      opserr << "WARNING invalid pinchY\n";
      opserr << "OOHysteretic material: " << tag << endln;
      return 0;
    }
  }

  HystereticBackbone *posBB = OPS_getHystereticBackbone(bTagPos);

  if (posBB == 0) {
    opserr << "WARNING backbone does not exist\n";
    opserr << "backbone: " << bTagPos;
    opserr << "\nuniaxialMaterial OOHystereitc: " << tag << endln;
    return 0;
  }
  UnloadingRule *posUnl = OPS_getUnloadingRule(unlTagPos);

  if (posUnl == 0) {
    opserr << "WARNING unloadingRule does not exist\n";
    opserr << "unloadingRule: " << unlTagPos;
    opserr << "\nuniaxialMaterial OOHystereitc: " << tag << endln;
    return 0;
  }

  StiffnessDegradation *posStf = OPS_getStiffnessDegradation(stfTagPos);

  if (posStf == 0) {
    opserr << "WARNING stiffnessDegradation does not exist\n";
    opserr << "stiffnessDegradation: " << stfTagPos;
    opserr << "\nuniaxialMaterial OOHystereitc: " << tag << endln;
    return 0;
  }

  StrengthDegradation *posStr = OPS_getStrengthDegradation(strTagPos);

  if (posStr == 0) {
    opserr << "WARNING strengthDegradation does not exist\n";
    opserr << "strengthDegradation: " << strTagPos;
    opserr << "\nuniaxialMaterial OOHystereitc: " << tag << endln;
    return 0;
  }
  if (argc > 8) {
    HystereticBackbone *negBB = OPS_getHystereticBackbone(bTagNeg);

    if (negBB == 0) {
      opserr << "WARNING backbone does not exist\n";
      opserr << "backbone: " << bTagNeg;
      opserr << "\nuniaxialMaterial OOHystereitc: " << tag << endln;
      return 0;
    }

    UnloadingRule *negUnl = OPS_getUnloadingRule(unlTagNeg);

    if (negUnl == 0) {
      opserr << "WARNING unloadingRule does not exist\n";
      opserr << "unloadingRule: " << unlTagNeg;
      opserr << "\nuniaxialMaterial OOHystereitc: " << tag << endln;
      return 0;
    }

    StiffnessDegradation *negStf = OPS_getStiffnessDegradation(stfTagNeg);

    if (negStf == 0) {
      opserr << "WARNING stiffnessDegradation does not exist\n";
      opserr << "stiffnessDegradation: " << stfTagNeg;
      opserr << "\nuniaxialMaterial OOHystereitc: " << tag << endln;
      return 0;
    }

    StrengthDegradation *negStr = OPS_getStrengthDegradation(strTagNeg);

    if (negStr == 0) {
      opserr << "WARNING strengthDegradation does not exist\n";
      opserr << "strengthDegradation: " << strTagNeg;
      opserr << "\nuniaxialMaterial OOHystereitc: " << tag << endln;
      return 0;
    }

    theMaterial =
        new OOHystereticMaterial(tag, *posBB, *negBB, *posUnl, *negUnl, *posStf,
                                 *negStf, *posStr, *negStr, pinchX, pinchY);
  } else {
    theMaterial = new OOHystereticMaterial(tag, *posBB, *posUnl, *posStf,
                                           *posStr, pinchX, pinchY);
  }

  // opserr << "OOHysteretic " << bTagPos << ' ' << unlTagPos << ' ' <<
  // stfTagPos << ' ' << strTagPos << ' ' << pinchX << ' ' << pinchY << endln;
  // opserr << "\t" << bTagNeg << ' ' << unlTagNeg << ' ' << stfTagNeg << ' ' <<
  // strTagNeg << ' ' << pinchX << ' ' << pinchY << endln;

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type "
              "OOHystereticMaterial\n";
    return 0;
  }

  return theMaterial;
}
