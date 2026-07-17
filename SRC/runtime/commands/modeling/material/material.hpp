//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
// |       Stress            |          Strain         |      Reduced Ordering     |
// | xx  yy   zz xy  yz  xz  | xx  yy   zz xy  yz  xz  |   0   1   2   3   4   5   |
// |  0   1   2   3  [.] [.] |  0  1   2   3  (.) (.)  |  xx  yy  zz  xy ......... | Axisymmetric      |
// |  0   1  (.)  2  (.) (.) |  0  1  [.]  2  [.] [.]  |  xx  yy  xy ............. | PlaneStress/Membrane |
// |  0   1  [.]  2  [.] [.] |  0  1  (.)  2  (.) (.)  |  xx  yy  xy ............. | PlaneStrain       |
// |  0   1  (.)  2  [.] [.] |  0  1  [.]  2  (.) (.)  |  xx  yy  xy ............. | ThinShell         |
// |  0  (.) (.)  1  (.)  2  |  0 [.] [.]  1  [.]  2   |  xx  xy  xz ............. | BeamFiber         |
// |  0  (.) (.)  1  [.] [.] |  0 [.] [.]  1  (.) (.)  |  xx  xy  ................ | BeamFiber2D02     |
// |  0  (.) (.)  1  (.) (.) |  0 [.] [.]  1  [.] [.]  |  xx  xy  ................ | BeamFiber2D       |
// |  0  (.) (.) [.] (.) [.] |  0 [.] [.] (.) [.] (.)  |  xx  .................... | Uniaxial          |
// |  0   1  (.)  2   3   4  |  0  1  [.]  2   3   4   |  xx  yy  xy  yz  xz ..... | ThickShell        |
// |  0   1   2   3   4   5  |  0  1   2   3   4   5   |  xx  yy  zz  xy  yz  xz   | ThreeDimensional  |
//
// (.): zero
// [.]: implicitly determined
//
// |        Explicit          |   Implicit Stress        |  Implicit Strain         |
// |   0   1   2   3   4   5  |   0   1   2   3   4   5  |   0   1   2   3   4   5  |
// |  xx  yy  zz  xy ........ |  yz  xz  ..............  |  ......................  | Generalized Plane Strain
// |  xx  yy  zz  xy ........ |  yz  xz  ..............  |  ......................  | Axisymmetric
// |  xx  yy  xy ............ |  ......................  |  zz  yz  xz ...........  | PlaneStress
// |  xx  yy  xy ............ |  zz  yz  xz ...........  |  ......................  | PlaneStrain
// |  xx  yy  xy ............ |  yz  xz  ..............  |  zz  ..................  | ThinShell
// |  xx  xy  xz ............ |  ......................  |  yy  zz  yz ...........  | 3D Shear Beam
// |  xx  xy  ............... |  ......................  |  yy  zz  yz  xz .......  | BeamFiber2D
// |  xx  xy  ............... |  yz  xz  ..............  |  yy  zz  ..............  | BeamFiber2D02
// |  xx  ................... |  xy  xz  ..............  |  yy  zz  yz ...........  | 2D/3D Euler Beam
// |  xx  yy  xy  yz  xz .... |  ......................  |  zz  ..................  | Mindlin Shell
// |  xx  ................... |  ......................  |  yy  zz  xy  yz  xz ...  | Uniaxial Stress
// |  xx  ................... |  yy  zz  xy  yz  xz ...  |  ......................  | Uniaxial Strain
// |  yz  xz  ............... |  xx  yy  zz  xy .......  |  ......................  | Anti-plane Strain
// |  xx  yy  zz  xy  yz  xz  |  ......................  |  ......................  | Continuum Solid
//
//
//
//
// strains ordered : eps11, eps22, eps33, 2*eps12, 2*eps23, 2*eps31 
// NDmaterial  strain order       = 11, 22, 33, 12, 23, 31 
// PlaneStress strain order       = 11, 22, 12, 33, 23, 31
// BeamFiber3D strain order       = 11, 12, 31, 22, 33, 23
// PlateFiber strain order        = 11, 22, 12, 23, 31, 33
// 
//                                   0   1   2   3   4   5

// Platefiber: 22, 33, 13, and 23 are condensed out.


//      | Stress     |  Strain      | 
//      0  1  2  3  4  5
//   ND : 11  22  33   12   23   31
//   PS : 11  22  12   33   23   31 
//   PF : 11  22  12   23   31   33 | 
//   F3 : 11  12  13 | 22   33   23 
//   F2 : 11  12     | 11   12      | 
//   AS : 11  22  33 12
//  
// ---------------------------------------
//
#include <tcl.h>
#include <string>
#include <assert.h>
#include <Parsing.h>
#include <unordered_map>
#include <elementAPI.h>
#include <runtimeAPI.h>
#include <NDMaterial.h>

extern Tcl_CmdProc TclCommand_newElasticMaterial;
extern Tcl_CmdProc TclCommand_newElasticOrthotropic;
// 
extern Tcl_CmdProc TclCommand_newJ2Material;
extern Tcl_CmdProc TclCommand_newPlasticMaterial;
extern Tcl_CmdProc TclCommand_newConcreteMaterial;
// concrete_asd.cpp
extern Tcl_CmdProc TclCommand_addASDConcrete3D;
// wrapper.cpp
extern Tcl_CmdProc TclCommand_addWrappingMaterial;
extern Tcl_CmdProc TclCommand_addParallel3DMaterial;
extern Tcl_CmdProc TclCommand_newPlateRebar;
extern Tcl_CmdProc TclCommand_newPlateFiber;
extern Tcl_CmdProc TclCommand_addPlaneWrapper;
extern Tcl_CmdProc TclCommand_addOrthotropicWrapper;

extern OPS_Routine OPS_ElasticOrthotropicPlaneStress;
extern OPS_Routine OPS_Series3DMaterial;
extern OPS_Routine OPS_J2CyclicBoundingSurfaceMaterial;
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
extern OPS_Routine OPS_NewPlasticDamageConcretePlaneStress;
extern OPS_Routine OPS_ElasticIsotropicMaterial;
extern OPS_Routine OPS_IncrementalElasticIsotropicThreeDimensional;
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
extern OPS_Routine OPS_StressDensityMaterial;
extern OPS_Routine OPS_PlaneStressLayeredMaterial;
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


extern "C"
int OPS_ResetInputNoBuilder(ClientData, Tcl_Interp *, int cArg,
                            int mArg, TCL_Char ** const argv, Domain *);


namespace {
template <OPS_Routine fn> static int
dispatch(ClientData clientData, Tcl_Interp* interp, int argc, G3_Char** const argv)
{
  ModelRegistry *builder = static_cast<ModelRegistry*>(clientData);

  OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, 0);

  G3_Runtime *rt = G3_getRuntime(interp);
  NDMaterial* theMaterial = (NDMaterial*)fn( rt, argc, argv );
  if (theMaterial == nullptr) {
    return TCL_ERROR;
  }

  if (builder->addTaggedObject<NDMaterial>(*theMaterial) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "Failed to add material to the model builder.\n";
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
}

namespace OpenSees {
namespace Library {

static std::unordered_map<std::string, Tcl_CmdProc*> MaterialLibrary = {
//
// Elastic 
//
// Isotropic
  {"ElasticIsotropic3D",               dispatch<TclCommand_newElasticMaterial>},
  {"ElasticIsotropic",                 dispatch<TclCommand_newElasticMaterial>},
  {"ElasticIsotropic3DThermal",        dispatch<OPS_ElasticIsotropicMaterialThermal>},
// Orthotropic
  {"ElasticOrthotropic",               dispatch<TclCommand_newElasticOrthotropic>},
  {"ElasticOrthotropicPlaneStress",    dispatch<OPS_ElasticOrthotropicPlaneStress>},
//
// Plasticity
//
  {"J2",                               dispatch<TclCommand_newPlasticMaterial>},
  {"J2Plasticity",                     dispatch<TclCommand_newPlasticMaterial>},
  {"GeneralizedJ2",                    dispatch<TclCommand_newPlasticMaterial>},

  {"PlasticJ2",                        dispatch<TclCommand_newPlasticMaterial>},
  {"NonlinearJ2",                      dispatch<TclCommand_newPlasticMaterial>},
  {"NonlinearJ2-UVC",                  dispatch<TclCommand_newPlasticMaterial>},

  {"J2N",                              dispatch<TclCommand_newPlasticMaterial>},
  {"J2L",                              dispatch<TclCommand_newPlasticMaterial>},
  {"J2Thermal",                        dispatch<TclCommand_newPlasticMaterial>},
  {"J2PlasticityThermal",              dispatch<TclCommand_newPlasticMaterial>},
  {"J2BeamFiber",                      dispatch<TclCommand_newPlasticMaterial>},
  {"J2BeamThread",                     dispatch<TclCommand_newPlasticMaterial>},

  {"SimplifiedJ2",                     dispatch<TclCommand_newPlasticMaterial>},
  {"J2Simplified",                     dispatch<TclCommand_newPlasticMaterial>},
  {"Simplified3DJ2",                   dispatch<TclCommand_newPlasticMaterial>},
  {"3DJ2",                             dispatch<TclCommand_newPlasticMaterial>},
  {"PlaneStressSimplifiedJ2",          dispatch<TclCommand_newPlasticMaterial>},
  {"DruckerPrager",                    dispatch<TclCommand_newPlasticMaterial>},

  {"UVCplanestress",                   dispatch<OPS_UVCplanestress       > },
  {"UVCmultiaxial",                    dispatch<OPS_UVCmultiaxial        > },
  {"J2PlateFibre",                     dispatch<TclCommand_newPlasticMaterial >},
//
  {"ManzariDafalias",                  dispatch<OPS_ManzariDafaliasMaterial>},
  {"ManzariDafaliasRO",                dispatch<OPS_ManzariDafaliasMaterialRO>},

//
  {"DruckerPragerThermal",             dispatch<OPS_DruckerPragerMaterialThermal> },
  {"TruncatedDP",                      dispatch<OPS_LinearCap     > },
  {"FSAM",                             dispatch<OPS_FSAMMaterial  > },
  {"AcousticMedium",                   dispatch<OPS_AcousticMedium> },
  {"CycLiqCP",                         dispatch<OPS_CycLiqCPMaterial>},
  {"CycLiqCPSP",                       dispatch<OPS_CycLiqCPSPMaterial>},
  {"BoundingCamClay",                  dispatch<OPS_BoundingCamClayMaterial>},
//
// Wrapper
//
  {"InitStrainMaterial",               dispatch<TclCommand_addWrappingMaterial>},
  {"InitStrain",                       dispatch<TclCommand_addWrappingMaterial>},
  {"InitialStrain",                    dispatch<TclCommand_addWrappingMaterial>},
  {"InitStressMaterial",               dispatch<OPS_InitStressNDMaterial>},
  {"Orthotropic",                      dispatch<TclCommand_addOrthotropicWrapper>},
  {"Series3DMaterial",                 dispatch<OPS_Series3DMaterial>},
  {"Parallel3DMaterial",               dispatch<TclCommand_addParallel3DMaterial>},
  {"Parallel3D",                       dispatch<TclCommand_addParallel3DMaterial>},
// Beam fiber (             22, 33, and 23 == 0)
  {"BeamFiber",                        dispatch<TclCommand_newPlateFiber>},
  {"BeamFiber2d",                      dispatch<TclCommand_newPlateFiber>},
  {"BeamFiber2dPS",                    dispatch<TclCommand_newPlateFiber>},
// Plane 
  {"PlaneStressMaterial",              dispatch<TclCommand_addPlaneWrapper>},
  {"PlaneStress",                      dispatch<TclCommand_addPlaneWrapper>},
  {"PlaneStrainMaterial",              dispatch<TclCommand_addPlaneWrapper>},
  {"PlaneStrain",                      dispatch<TclCommand_addPlaneWrapper>},
  {"PlaneStressRebarMaterial",         dispatch<TclCommand_newPlateRebar>},
// Plate  (constrain stress 33 == 13 == 23 == 0) 
  {"PlateRebarMaterial",               dispatch<TclCommand_newPlateRebar>},
  {"PlateRebar",                       dispatch<TclCommand_newPlateRebar>},
  {"PlateFiberMaterial",               dispatch<TclCommand_newPlateFiber>},
  {"PlateFiber",                       dispatch<TclCommand_newPlateFiber>},
//
// Other
//
  {"ReinforcedConcretePlaneStress",    dispatch<OPS_ReinforcedConcretePlaneStressMaterial>},
  {"PlaneStressLayeredMaterial",       dispatch<OPS_PlaneStressLayeredMaterial>},
  {"ASDConcrete3D",                    dispatch<TclCommand_addASDConcrete3D>},
  {"PlasticDamageConcrete",            dispatch<TclCommand_newConcreteMaterial>},
  {"FariaPlasticDamage",               dispatch<TclCommand_newConcreteMaterial>},
  {"PlasticDamageConcretePlaneStress", dispatch<OPS_NewPlasticDamageConcretePlaneStress>},
};

} // namespace Library

static std::unordered_map<std::string, OPS_Routine*> OldMaterialCommands = {
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
#ifdef _HAVE_Damage2p
  {"Damage2p",                        OPS_Damage2p},
#endif

  {"FAReinforcedConcretePlaneStress", OPS_FAReinforcedConcretePlaneStressMaterial},
  {"RAFourSteelRCPlaneStress",        OPS_RAFourSteelRCPlaneStressMaterial},
  {"FAFourSteelRCPlaneStress",        OPS_FAFourSteelRCPlaneStressMaterial},


  {"PrestressedConcretePlaneStress",   OPS_PrestressedConcretePlaneStressMaterial},
  {"FAPrestressedConcretePlaneStress", OPS_FAPrestressedConcretePlaneStressMaterial},
  {"RAFourSteetPCPlaneStress",         OPS_RAFourSteelPCPlaneStressMaterial},
  {"FAFourSteelPCPlaneStress",         OPS_FAFourSteelPCPlaneStressMaterial},

//{"MaterialCMM",    OPS_MaterialCMM},

  {"PM4Sand",                       OPS_PM4SandMaterial},
  {"J2CyclicBoundingSurface",       OPS_J2CyclicBoundingSurfaceMaterial},
  {"PM4Silt",                       OPS_PM4SiltMaterial},
  {"ContactMaterial2D",             OPS_ContactMaterial2DMaterial},
  {"ContactMaterial3D",             OPS_ContactMaterial3DMaterial},
  {"InitialStateAnalysisWrapper",   OPS_InitialStateAnalysisWrapperMaterial},
  {"stressDensity",                 OPS_StressDensityMaterial},
  {"IncrementalElasticIsotropic3D", OPS_IncrementalElasticIsotropicThreeDimensional},
  {"OrthotropicRAConcrete",         OPS_OrthotropicRotatingAngleConcreteT2DMaterial01},
  {"SmearedSteelDoubleLayer",       OPS_SmearedSteelDoubleLayerT2DMaterial01},
  {"SAniSandMS",                    OPS_SAniSandMSMaterial},
};
} // namespace OpenSees
