//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
#include <tcl.h>
#include <string>
#include <assert.h>
#include <Parsing.h>
#include <unordered_map>
#include <elementAPI.h>
#include <NDMaterial.h>

static Tcl_CmdProc TclCommand_addPlaneWrapper;
extern Tcl_CmdProc TclCommand_newJ2Material;
extern Tcl_CmdProc TclCommand_newJ2Simplified;
extern Tcl_CmdProc TclCommand_newPlasticMaterial;
extern Tcl_CmdProc TclCommand_newElasticMaterial;

extern OPS_Routine OPS_ElasticOrthotropicPlaneStress;
extern OPS_Routine OPS_OrthotropicMaterial;
extern OPS_Routine OPS_Series3DMaterial;
extern OPS_Routine OPS_Parallel3DMaterial;
extern OPS_Routine OPS_J2Plasticity;
extern OPS_Routine OPS_J2PlateFibreMaterial;
extern OPS_Routine OPS_J2CyclicBoundingSurfaceMaterial;
extern OPS_Routine OPS_ASDConcrete3DMaterial;
extern OPS_Routine OPS_ReinforcedConcretePlaneStressMaterial;
extern OPS_Routine OPS_FAReinforcedConcretePlaneStressMaterial;
extern OPS_Routine OPS_FAFourSteelRCPlaneStressMaterial;
extern OPS_Routine OPS_RAFourSteelRCPlaneStressMaterial;
extern OPS_Routine OPS_PrestressedConcretePlaneStressMaterial;
extern OPS_Routine OPS_FAPrestressedConcretePlaneStressMaterial;
extern OPS_Routine OPS_FAFourSteelPCPlaneStressMaterial;
extern OPS_Routine OPS_RAFourSteelPCPlaneStressMaterial;
// extern OPS_Routine OPS_MaterialCMM;
// extern OPS_Routine OPS_NewMaterialCMM;
extern OPS_Routine OPS_NewPlasticDamageConcrete3d;
extern OPS_Routine OPS_NewPlasticDamageConcretePlaneStress;
extern OPS_Routine OPS_ElasticIsotropicMaterial;
extern OPS_Routine OPS_ElasticIsotropic3D;
extern OPS_Routine OPS_IncrementalElasticIsotropicThreeDimensional;
extern OPS_Routine OPS_ElasticOrthotropicMaterial;
extern OPS_Routine OPS_DruckerPragerMaterial;
extern OPS_Routine OPS_BoundingCamClayMaterial;
extern OPS_Routine OPS_ContactMaterial2DMaterial;
extern OPS_Routine OPS_ContactMaterial3DMaterial;
extern OPS_Routine OPS_InitialStateAnalysisWrapperMaterial;
extern OPS_Routine OPS_ManzariDafaliasMaterial;
extern OPS_Routine OPS_ManzariDafaliasMaterialRO;
extern OPS_Routine OPS_PM4SandMaterial;
extern OPS_Routine OPS_PM4SiltMaterial;
extern OPS_Routine OPS_CycLiqCPMaterial;
extern OPS_Routine OPS_CycLiqCPSPMaterial;
extern OPS_Routine OPS_InitStressNDMaterial;
extern OPS_Routine OPS_InitStrainNDMaterial;
extern OPS_Routine OPS_StressDensityMaterial;
extern OPS_Routine OPS_PlaneStressLayeredMaterial;
extern OPS_Routine OPS_PlaneStressRebarMaterial;
extern OPS_Routine OPS_PlateFiberMaterial;
extern OPS_Routine OPS_BeamFiberMaterial;
extern OPS_Routine OPS_BeamFiberMaterial2d;
extern OPS_Routine OPS_BeamFiberMaterial2dPS;
extern OPS_Routine OPS_LinearCap;
extern OPS_Routine OPS_AcousticMedium;
extern OPS_Routine OPS_UVCmultiaxial;
extern OPS_Routine OPS_UVCplanestress;
extern OPS_Routine OPS_SAniSandMSMaterial;
extern OPS_Routine OPS_OrthotropicRotatingAngleConcreteT2DMaterial01;	// M. J. Nunez - UChile
extern OPS_Routine OPS_SmearedSteelDoubleLayerT2DMaterial01;		// M. J. Nunez - UChile

extern OPS_Routine OPS_ElasticIsotropicMaterialThermal;           // L.Jiang [SIF]
extern OPS_Routine OPS_DruckerPragerMaterialThermal;              // L.Jiang [SIF]
extern OPS_Routine OPS_PlasticDamageConcretePlaneStressThermal;   // L.Jiang [SIF]


extern OPS_Routine OPS_AllASDPlasticMaterials;

#ifdef _HAVE_Faria1998
extern OPS_Routine OPS_NewFaria1998Material;
extern OPS_Routine OPS_NewConcreteMaterial;
#endif

extern OPS_Routine OPS_FSAMMaterial; // K Kolozvari

#ifdef _HAVE_Damage2p
extern OPS_Routine OPS_Damage2p;
#endif


extern "C" int OPS_ResetInputNoBuilder(ClientData clientData, Tcl_Interp *interp, int cArg,
  int mArg, TCL_Char ** const argv, Domain *domain);

template <OPS_Routine fn> static int
dispatch(ClientData clientData, Tcl_Interp* interp, int argc, G3_Char** const argv)
{
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

  OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, 0);

  G3_Runtime *rt = G3_getRuntime(interp);
  NDMaterial* theMaterial = (NDMaterial*)fn( rt, argc, argv );
  if (theMaterial == nullptr) {
    return TCL_ERROR;
  }

  if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "Failed to add material to the model builder.\n";
    delete theMaterial;
    return TCL_ERROR;
  }
  return TCL_OK;
}

template <int (*fn)(ClientData clientData, Tcl_Interp* interp, int, G3_Char** const)> 
static int
dispatch(ClientData clientData, Tcl_Interp* interp, int argc, G3_Char** const argv)
{
  assert(clientData != nullptr);
  return fn( clientData, interp, argc, argv );
}

static std::unordered_map<std::string, Tcl_CmdProc*> material_dispatch2 = {
//
// Elastic 
//
// Isotropic
  {"ElasticIsotropic3D",               dispatch<OPS_ElasticIsotropic3D>},
  {"ElasticIsotropic",                 dispatch<TclCommand_newElasticMaterial>},
  {"ElasticIsotropic3DThermal",        dispatch<OPS_ElasticIsotropicMaterialThermal>},
// Orthotropic
  {"ElasticOrthotropic",               dispatch<OPS_ElasticOrthotropicMaterial>},
  {"ElasticOrthotropicPlaneStress",    dispatch<OPS_ElasticOrthotropicPlaneStress>},

//
// Plasticity
//
  {"J2",                               dispatch<OPS_J2Plasticity>},
  {"J2Plasticity",                     dispatch<OPS_J2Plasticity>},
  {"SimplifiedJ2",                     dispatch<TclCommand_newPlasticMaterial>},
  {"J2Simplified",                     dispatch<TclCommand_newPlasticMaterial>},
  {"Simplified3DJ2",                   dispatch<TclCommand_newJ2Simplified>},
  {"3DJ2",                             dispatch<TclCommand_newJ2Simplified>},
  {"PlaneStressSimplifiedJ2",          dispatch<TclCommand_newJ2Simplified>},

  {"UVCplanestress",                   dispatch<OPS_UVCplanestress       > },
  {"UVCmultiaxial",                    dispatch<OPS_UVCmultiaxial        > },
  {"J2PlateFibre",                     dispatch<OPS_J2PlateFibreMaterial>}, 
  {"PlateFiber",                       dispatch<OPS_PlateFiberMaterial>},
//
  {"ManzariDafalias",                  dispatch<OPS_ManzariDafaliasMaterial>},
  {"ManzariDafaliasRO",                dispatch<OPS_ManzariDafaliasMaterialRO>},

  // Beam fiber
  {"BeamFiber",                        dispatch<OPS_BeamFiberMaterial> },
  {"BeamFiber2d",                      dispatch<OPS_BeamFiberMaterial2d> },
  {"BeamFiber2dPS",                    dispatch<OPS_BeamFiberMaterial2dPS> },
  {"DruckerPragerThermal",             dispatch<OPS_DruckerPragerMaterialThermal> },
  {"DruckerPrager",                    dispatch<OPS_DruckerPragerMaterial> },
  {"TruncatedDP",                      dispatch<OPS_LinearCap            > },
  {"FSAM",                             dispatch<OPS_FSAMMaterial         > },
  {"AcousticMedium",                   dispatch<OPS_AcousticMedium       > },
  {"CycLiqCP",                         dispatch<OPS_CycLiqCPMaterial>},
  {"CycLiqCPSP",                       dispatch<OPS_CycLiqCPSPMaterial>},
  {"BoundingCamClay",                  dispatch<OPS_BoundingCamClayMaterial>},
//
// Wrapper
//
  {"PlaneStressMaterial",              dispatch<TclCommand_addPlaneWrapper>},
  {"PlaneStress",                      dispatch<TclCommand_addPlaneWrapper>},
  {"PlaneStrainMaterial",              dispatch<TclCommand_addPlaneWrapper>},
  {"PlaneStrain",                      dispatch<TclCommand_addPlaneWrapper>},
//
// Other

  {"InitStressMaterial",               dispatch<OPS_InitStressNDMaterial>},
  {"InitStrainMaterial",               dispatch<OPS_InitStrainNDMaterial>},
  {"InitStrain",                       dispatch<OPS_InitStrainNDMaterial>},
  {"ReinforcedConcretePlaneStress",    dispatch<OPS_ReinforcedConcretePlaneStressMaterial>},
  {"PlaneStressLayeredMaterial",       dispatch<OPS_PlaneStressLayeredMaterial>},
  {"PlaneStressRebarMaterial",         dispatch<OPS_PlaneStressRebarMaterial>},
  {"ASDConcrete3D",                    dispatch<OPS_ASDConcrete3DMaterial>},
  {"PlasticDamageConcrete",            dispatch<OPS_NewPlasticDamageConcrete3d>},
  {"PlasticDamageConcretePlaneStress", dispatch<OPS_NewPlasticDamageConcretePlaneStress>},
};

static std::unordered_map<std::string, OPS_Routine*> material_dispatch = {

#ifdef OPS_USE_ASDPlasticMaterials
  {"ASDPlasticMaterial",            OPS_AllASDPlasticMaterials},
#endif

#if 0
  {"CDPPlaneStressThermal", OPS_PlasticDamageConcretePlaneStressThermal},
#endif

#ifdef _HAVE_Faria1998
  {"Faria1998", OPS_NewFaria1998Material},  
  {"Concrete", OPS_NewConcreteMaterial},
#endif

  {"FAReinforcedConcretePlaneStress", OPS_FAReinforcedConcretePlaneStressMaterial},
  {"RAFourSteelRCPlaneStress",        OPS_RAFourSteelRCPlaneStressMaterial},
  {"FAFourSteelRCPlaneStress",        OPS_FAFourSteelRCPlaneStressMaterial},

#ifdef _HAVE_Damage2p
  {"Damage2p",                        OPS_Damage2p},
#endif

  {"PrestressedConcretePlaneStress",   OPS_PrestressedConcretePlaneStressMaterial},
  {"FAPrestressedConcretePlaneStress", OPS_FAPrestressedConcretePlaneStressMaterial},
  {"RAFourSteetPCPlaneStress",         OPS_RAFourSteelPCPlaneStressMaterial},

  {"FAFourSteelPCPlaneStress",         OPS_FAFourSteelPCPlaneStressMaterial},


//{"MaterialCMM",    OPS_MaterialCMM},

  {"PM4Sand", OPS_PM4SandMaterial},

  {"J2CyclicBoundingSurface", OPS_J2CyclicBoundingSurfaceMaterial},

  {"PM4Silt", OPS_PM4SiltMaterial},

  {"ContactMaterial2D", OPS_ContactMaterial2DMaterial},

  {"ContactMaterial3D", OPS_ContactMaterial3DMaterial},

  {"InitialStateAnalysisWrapper", OPS_InitialStateAnalysisWrapperMaterial},

  {"stressDensity", OPS_StressDensityMaterial},

  {"IncrementalElasticIsotropic3D", OPS_IncrementalElasticIsotropicThreeDimensional},

  {"OrthotropicRAConcrete",         OPS_OrthotropicRotatingAngleConcreteT2DMaterial01},
  {"SmearedSteelDoubleLayer",       OPS_SmearedSteelDoubleLayerT2DMaterial01},

  {"SAniSandMS",                    OPS_SAniSandMSMaterial},
//
// Wrapper
//
  {"OrthotropicMaterial",           OPS_OrthotropicMaterial},
  {"Series3DMaterial",              OPS_Series3DMaterial},
  {"Parallel3DMaterial",            OPS_Parallel3DMaterial},
  {"Parallel3D",                    OPS_Parallel3DMaterial},
};

