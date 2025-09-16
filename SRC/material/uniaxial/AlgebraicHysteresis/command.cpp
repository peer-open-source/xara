#include <OPS_Globals.h>
#include <Logging.h>
#include <Parsing.h>

#include <ModelRegistry.h>
#include <UnloadingRule.h>
#include <HystereticBackbone.h>
#include <OOHystereticMaterial.h>
#include <StiffnessDegradation.h>
#include <StrengthDegradation.h>

#include <elementAPI.h>

int 
Create_OOHystereticMaterial(ClientData clientData,
                                   Tcl_Interp *interp, 
                                   Tcl_Size argc,
                                   TCL_Char ** const argv)
{
  ModelRegistry *builder = static_cast<ModelRegistry *>(clientData);
  UniaxialMaterial *theMaterial = nullptr;

  if (argc - 2 < 5) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial OOHysteretic tag? bTag+? unlRulTag+? stfDegTag+? strDegTag+? "
	   << "<bTag-? unlRulTag-? stfDegTag-? strDegTag-?> <pinchX? pinchY?>" << "\n";
    return TCL_ERROR;
  }
  
  int tag;
  int bTagPos, bTagNeg;
  int unlTagPos, unlTagNeg;
  int stfTagPos, stfTagNeg;
  int strTagPos, strTagNeg;
  double pinchX = 0.0;
  double pinchY = 1.0;

  argc -= 1; // remove the uniaxialMaterial keyword
  argc -= 1; // remove the OOHysteretic keyword
  
  int i = 3;
  int numData = 1;
  if (Tcl_GetInt(interp, argv[i++], &tag) != TCL_OK) {
    opserr << "WARNING invalid tag\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[i++], &bTagPos) != TCL_OK) {
    opserr << "WARNING invalid bTag+\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[i++], &unlTagPos) != TCL_OK) {
    opserr << "WARNING invalid unlRulTag+\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[i++], &stfTagPos) != TCL_OK) {
    opserr << "WARNING invalid stfDegTag+\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[i++], &strTagPos) != TCL_OK) {
    opserr << "WARNING invalid strDegTag+\n";
    return TCL_ERROR;
  }

  if (argc == 7) {
    if (OPS_GetDoubleInput(&numData, &pinchX) != 0) {
      opserr << "WARNING invalid pinchX\n";
      return TCL_ERROR;
    }
    if (OPS_GetDoubleInput(&numData, &pinchY) != 0) {
      opserr << "WARNING invalid pinchY\n";
      return TCL_ERROR;
    }
  }

  if (argc > 8) {
    if (Tcl_GetInt(interp, argv[i++], &bTagNeg) != TCL_OK) {
      opserr << "WARNING invalid bTag-\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[i++], &unlTagNeg) != TCL_OK) {
      opserr << "WARNING invalid unlRulTag-\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[i++], &stfTagNeg) != TCL_OK) {
      opserr << "WARNING invalid stfDegTag-\n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[i++], &strTagNeg) != TCL_OK) {
      opserr << "WARNING invalid strDegTag-\n";
      return TCL_ERROR;
    }
  }
  if (argc == 11) {
    if (OPS_GetDoubleInput(&numData, &pinchX) != 0) {
      opserr << "WARNING invalid pinchX\n";
      return TCL_ERROR;
    }
    if (OPS_GetDoubleInput(&numData, &pinchY) != 0) {
      opserr << "WARNING invalid pinchY\n";
      return TCL_ERROR;
    }
  }
  
  HystereticBackbone *posBB = builder->getTypedObject<HystereticBackbone>(bTagPos);
  
  if (posBB == 0) {
    opserr << "WARNING backbone does not exist\n";
    opserr << "backbone: " << bTagPos; 
    return TCL_ERROR;
  }
  UnloadingRule *posUnl = builder->getTypedObject<UnloadingRule>(unlTagPos);
  
  if (posUnl == 0) {
    opserr << "WARNING unloadingRule does not exist\n";
    opserr << "unloadingRule: " << unlTagPos; 
    return TCL_ERROR;
  }
  
  StiffnessDegradation *posStf = builder->getTypedObject<StiffnessDegradation>(stfTagPos);
  
  if (posStf == nullptr) {
    opserr << "WARNING stiffnessDegradation does not exist\n";
    opserr << "stiffnessDegradation: " << stfTagPos; 
    return TCL_ERROR;
  }

  StrengthDegradation *posStr = builder->getTypedObject<StrengthDegradation>(strTagPos);

  if (posStr == 0) {
    opserr << "WARNING strengthDegradation does not exist\n";
    opserr << "strengthDegradation: " << strTagPos; 
    return TCL_ERROR;
  }

  if (argc > 8) {
    HystereticBackbone *negBB = OPS_getHystereticBackbone(bTagNeg);
    
    if (negBB == 0) {
      opserr << "WARNING backbone does not exist\n";
      opserr << "backbone: " << bTagNeg; 
      opserr << "\nuniaxialMaterial OOHystereitc: " << tag << endln;
      return TCL_ERROR;
    }
    
    UnloadingRule *negUnl = OPS_getUnloadingRule(unlTagNeg);
    
    if (negUnl == 0) {
      opserr << "WARNING unloadingRule does not exist\n";
      opserr << "unloadingRule: " << unlTagNeg; 
      opserr << "\nuniaxialMaterial OOHystereitc: " << tag << endln;
      return TCL_ERROR;
    }
    
    StiffnessDegradation *negStf = builder->getTypedObject<StiffnessDegradation>(stfTagNeg);
    
    if (negStf == 0) {
      opserr << "WARNING stiffnessDegradation does not exist\n";
      opserr << "stiffnessDegradation: " << stfTagNeg; 
      opserr << "\nuniaxialMaterial OOHystereitc: " << tag << endln;
      return TCL_ERROR;
    }
    
    StrengthDegradation *negStr = builder->getTypedObject<StrengthDegradation>(strTagNeg);
    
    if (negStr == 0) {
      opserr << "WARNING strengthDegradation does not exist\n";
      opserr << "strengthDegradation: " << strTagNeg; 
      opserr << "\nuniaxialMaterial OOHystereitc: " << tag << endln;
      return TCL_ERROR;
    }
    
    theMaterial = 
      new OOHystereticMaterial(tag, *posBB, *negBB, *posUnl, *negUnl, 
			       *posStf, *negStf, *posStr, *negStr,
			       pinchX, pinchY);
  }
  else {
    theMaterial = 
      new OOHystereticMaterial(tag, *posBB, *posUnl, *posStf, *posStr,
			       pinchX, pinchY);
  }

  if (theMaterial == nullptr) {
    opserr << "WARNING could not create uniaxialMaterial of type OOHystereticMaterial\n";
    return TCL_ERROR;
  }

  if (builder->addTaggedObject<UniaxialMaterial>(*theMaterial) != TCL_OK) {
    opserr << "WARNING could not add uniaxialMaterial of type OOHystereticMaterial to the domain\n";
    opserr << *theMaterial << endln;
    delete theMaterial;
    return TCL_ERROR;
  }

  return TCL_OK;
  
}