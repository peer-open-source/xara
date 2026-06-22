//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Written: cmp
//
#include <stdio.h>
#include <string.h>
#include <Parsing.h>
#include <Logging.h>

#include "DegradingUniaxialWrapper.h"
#include <ModelRegistry.h>

#define WRAPPER_CMD "UniaxialDamage"
// #define WRAPPER_CMD "FedeasDamage"


// dmg::load type tag
static int
NewLoading_Cmd([[maybe_unused]] ClientData clientData, Tcl_Interp *interp, int objc,
               Tcl_Obj *const objv[])
{
 
  double E, fy;
  int argi = 0;

  if (objc < 3){
    opserr << OpenSees::PromptValueError << "insufficient parameters.\n";
    return TCL_ERROR;
  }

  DegradingUniaxialWrapper::Data *data = static_cast<DegradingUniaxialWrapper::Data *>(clientData);

  DamageIndex *load = nullptr;
  bool iso = false;
  const char *dir = Tcl_GetString(objv[1]);
  if ((strcmp(dir, "pos") == 0) || (strcmp(dir, "iso") == 0)) {
    load = &data->idx[0];
    if (strcmp(dir, "iso") == 0)
      iso = true;
  } else if (strcmp(dir, "neg") == 0) {
    load = &data->idx[1];
  } else {
    opserr << OpenSees::PromptValueError 
           << "unknown loading direction '" << dir << "'.\n";
    return TCL_ERROR;
  }


  const char *typ = Tcl_GetString(objv[2]);

  if (strcmp(typ, "none") == 0) {
    load->type = DamageIndex::Type::None;
  } else if (strcmp(typ, "mbeta") == 0) {
    load->type = DamageIndex::Type::mBeta;
  } else if (strcmp(typ, "obeta") == 0) {
    load->type = DamageIndex::Type::oBeta;
  } else if (strcmp(typ, "lognormal") == 0) {
    load->type = DamageIndex::Type::Log;
  } else if (strcmp(typ, "weibull") == 0) {
    load->type = DamageIndex::Type::Weibull;
  } else if (strcmp(typ, "bilinear") == 0) {
    load->type = DamageIndex::Type::Bilinear;
  } else if (strcmp(typ, "trilinear") == 0) {
    load->type = DamageIndex::Type::Trilinear;
  }
  else {
    opserr << OpenSees::PromptValueError 
           << "unknown evolution type '" << typ << "'.\n";
    return TCL_ERROR;
  }

  {
    TCL_Char** evol_argv;
    int evol_argc;
    if (Tcl_SplitList(interp, Tcl_GetString(objv[3]), &evol_argc, &evol_argv) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid evolution parameters list.\n";
      return TCL_ERROR;
    }

    for (int i=0; i<evol_argc; i++) {
      if (strcmp(evol_argv[i], "-frac") == 0) {
        // TODO
        return TCL_ERROR;
      } else {
        if (Tcl_GetDouble(interp, evol_argv[i], load->evol.theParam + load->evol.nParam++) != TCL_OK)
          return TCL_ERROR;
      }
    }
  }

  for (; argi < objc; argi++) {
    auto argis = Tcl_GetString(objv[argi]);

    if (strcmp(argis, "-frac") == 0) {
      opserr << OpenSees::PromptValueError << "fracture data not implemented yet.\n";
      return TCL_ERROR;
    }
    else if (strcmp(argis, "-Cd0") == 0) {
      if (Tcl_GetDoubleFromObj(interp, objv[++argi], &load->Cd0) != TCL_OK)
        return TCL_ERROR;
    } 
    else if ( (strcmp(argis, "-ccd") == 0)    ||
              (strcmp(argis, "-Ccd") == 0))   {
      if (Tcl_GetDoubleFromObj(interp, objv[++argi], &load->Ccd) != TCL_OK)
        return TCL_ERROR;
    }
    else if (strcmp(Tcl_GetString(objv[argi]), "-Psi_d0") == 0) {
      if (Tcl_GetDoubleFromObj(interp, objv[++argi], &load->psi_d0) != TCL_OK)
        return TCL_ERROR;

    } else if (strcmp(Tcl_GetString(objv[argi]), "-Cd1") == 0) {
      if (Tcl_GetDoubleFromObj(interp, objv[++argi], &load->Cd1) != TCL_OK)
        return TCL_ERROR;

    } else if (strcmp(Tcl_GetString(objv[argi]), "-Psi_d1") == 0) {
      if (Tcl_GetDoubleFromObj(interp, objv[++argi], &load->psi_d1) != TCL_OK)
        return TCL_ERROR;

    } else if (strcmp(Tcl_GetString(objv[argi]), "-E") == 0) {
      if (Tcl_GetDoubleFromObj(interp, objv[++argi], &E) != TCL_OK)
        return TCL_ERROR;

    } else if (strcmp(Tcl_GetString(objv[argi]), "-fy") == 0) {
      if (Tcl_GetDoubleFromObj(interp, objv[++argi], &fy) != TCL_OK)
        return TCL_ERROR;

    } else if (strcmp(Tcl_GetString(objv[argi]), "-Cwc") == 0) {
      if (Tcl_GetDoubleFromObj(interp, objv[++argi], &load->Cwc) != TCL_OK)
        return TCL_ERROR;
    }
  }
  if (load->psi_d0 == 0.0 && load->psi_d1 == 0.0
      && E > 0.0 && fy > 0.0){
    double yield_energy = fy*fy * 0.5 / E;
    load->psi_d0 = load->Cd0 * yield_energy;
    load->psi_d1 = load->Cd1 * yield_energy;
  }

  // TODO
  load->frac = DmgFrac {false, 1.0, 1.0};

  load->valid = true;
  
  if (iso) {
    data->idx[1] = *load;
  }

  return TCL_OK;
}



int
TclCommand_newFedeasUniaxialDamage(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  ModelRegistry* builder = static_cast<ModelRegistry*>(clientData);

  if (argc < 2) {
    opserr << "WARNING invalid uniaxialMaterial " WRAPPER_CMD " $tag "
              "$wrapTag <-damage $damageTag>"
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;;
  }

  

  // Get wrapper tag
  int tags[2];
  if (Tcl_GetInt(interp, argv[2], &tags[0]) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial tag" << OpenSees::SignalMessageEnd;
    return TCL_ERROR;;
  }
  // Get base tag
  if (Tcl_GetInt(interp, argv[3], &tags[1]) != TCL_OK) {
    opserr << "WARNING invalid uniaxialMaterial tag" << OpenSees::SignalMessageEnd;
    return TCL_ERROR;;
  }

  // Get base material
  UniaxialMaterial *theWrappedMaterial = builder->getTypedObject<UniaxialMaterial>(tags[1]);
  if (theWrappedMaterial == nullptr) {
    opserr << "WARNING unable to retrieve uniaxialMaterial with tag" WRAPPER_CMD " tag: "
           << tags[1] << OpenSees::SignalMessageEnd;
    return TCL_ERROR;;
  }

  DegradingUniaxialWrapper::Data data;

  int argn = 4;
  double Ccd = 0.5;
  while (argn < argc) {
    const char *param = argv[argn];

    if ((strcmp(param, "-damage") == 0) || 
        (strcmp(param, "-evol") == 0)   ||
        (strcmp(param, "-evolution") == 0))   {

      Tcl_CreateObjCommand(interp, "dmg::evol", NewLoading_Cmd,   (void*)&data, NULL);

      Tcl_Eval(interp, argv[++argn]);

      Tcl_DeleteCommand(interp, "dmg::evol");

    } else {
      break;
    }
    argn++;
  }

  if (data.idx[0].valid == false || data.idx[1].valid == false) {
    opserr << "WARNING no damage data provided for uniaxialMaterial "
           << tags[0] 
           << OpenSees::SignalMessageEnd;
    return TCL_ERROR;;
  }

  // Parsing was successful, allocate the material
  DegradingUniaxialWrapper *theMaterial =
    new DegradingUniaxialWrapper(tags[0], *theWrappedMaterial, data);


  if (builder->addTaggedObject<UniaxialMaterial>(*theMaterial) != TCL_OK) {
    delete theMaterial;
    return TCL_ERROR;;
  }
  return TCL_OK;
}