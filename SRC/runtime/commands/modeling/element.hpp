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
#include <CommandLibrary.h>
#include <string>
#include <unordered_map>
#include <algorithm>

class G3_Runtime;
typedef void *OPS_Routine(G3_Runtime* , int, const char** const);

extern OPS_Routine OPS_ElasticTubularJoint;
extern OPS_Routine OPS_ZeroLengthContactNTS2D;
extern OPS_Routine OPS_ZeroLengthVG_HG;
extern OPS_Routine OPS_ZeroLengthInterface2D;
extern OPS_Routine OPS_ZeroLengthImpact3D;
extern OPS_Routine OPS_ZeroLengthContactASDimplex; 
extern "C" OPS_Routine OPS_PY_Macro2D;
extern OPS_Routine OPS_SimpleContact2D;
extern OPS_Routine OPS_SimpleContact3D;

extern OPS_Routine OPS_SurfaceLoad;
extern OPS_Routine OPS_TriSurfaceLoad;

extern OPS_Routine OPS_ModElasticBeam2d;
extern OPS_Routine OPS_ModElasticBeam3d;


extern OPS_Routine OPS_TPB1D;
extern OPS_Routine OPS_TFP_Bearing;
extern OPS_Routine OPS_FPBearingPTV;
extern OPS_Routine OPS_MultiFP2d;
extern OPS_Routine OPS_CoupledZeroLength;
extern OPS_Routine OPS_FourNodeQuad3d;
extern OPS_Routine OPS_Quad4FiberOverlay;
extern OPS_Routine OPS_QuadBeamEmbedContact;
extern OPS_Routine OPS_ASID8QuadWithSensitivity;
extern OPS_Routine OPS_AV3D4QuadWithSensitivity;

extern OPS_Routine OPS_Brick8FiberOverlay;
extern OPS_Routine OPS_PML3D;
extern OPS_Routine OPS_PML2D;
extern OPS_Routine OPS_CorotTruss2;
extern OPS_Routine OPS_N4BiaxialTruss;
extern OPS_Routine OPS_AC3D8HexWithSensitivity;
extern OPS_Routine OPS_VS3D4QuadWithSensitivity;
extern OPS_Routine OPS_MVLEM;        // Kristijan Kolozvari
extern OPS_Routine OPS_SFI_MVLEM;    // Kristijan Kolozvari
extern OPS_Routine OPS_MVLEM_3D;     // Kristijan Kolozvari
extern OPS_Routine OPS_SFI_MVLEM_3D; // Kristijan Kolozvari
extern OPS_Routine OPS_E_SFI_MVLEM_3D;
extern OPS_Routine OPS_E_SFI;
extern OPS_Routine OPS_MEFI;
extern OPS_Routine OPS_InertiaTrussElement; // Added by Xiaodong Ji, Yuhao Cheng, Yue Yu
extern OPS_Routine OPS_CatenaryCableElement;
extern OPS_Routine OPS_LysmerTriangle;
extern OPS_Routine OPS_ASDEmbeddedNodeElement;     // Massimo Petracca (ASDEA)
extern OPS_Routine OPS_ASDAbsorbingBoundary2D;     // Massimo Petracca (ASDEA)
extern OPS_Routine OPS_ASDAbsorbingBoundary3D;     // Massimo Petracca (ASDEA)
// Fluid
extern OPS_Routine OPS_FSIInterfaceElement2D;      // Massimo Petracca (ASDEA)
extern OPS_Routine OPS_FSIFluidBoundaryElement2D;  // Massimo Petracca (ASDEA)
extern OPS_Routine OPS_FSIFluidElement2D;          // Massimo Petracca (ASDEA)
// Other
extern OPS_Routine OPS_LinearElasticSpring;
extern OPS_Routine OPS_Inerter;
extern OPS_Routine OPS_Inno3DPnPJoint;
extern OPS_Routine OPS_Adapter;
extern OPS_Routine OPS_Actuator;
extern OPS_Routine OPS_ActuatorCorot;

extern OPS_Routine OPS_RJWatsonEQS2d;
extern OPS_Routine OPS_RJWatsonEQS3d;
extern OPS_Routine OPS_RockingBC;
extern OPS_Routine OPS_LehighJoint2d;
extern OPS_Routine OPS_MasonPan12;
extern OPS_Routine OPS_MasonPan3D;


namespace {
static
std::string toLower( const std::string & s )
{
  std::string copy = s;
  transform( copy.begin( ), copy.end( ), copy.begin( ), 
      [](unsigned char c) { return std::tolower(c); });
  return copy;
}

static bool 
equalsIgnoreCase( const std::string & lhs, const std::string & rhs )
{
  return toLower( lhs ) == toLower( rhs );
}

class CaseInsensitive
{
public:
  size_t operator( ) ( const std::string & s ) const
  {  
      static std::hash<std::string> hf;
      return hf( toLower( s ) );
  }
  
  bool operator( ) ( const std::string & lhs, const std::string & rhs ) const
  {
    return equalsIgnoreCase( lhs, rhs );
  }
};
}

// Zero-length
Tcl_CmdProc TclCommand_addZeroLength;
Tcl_CmdProc TclCommand_addZeroLengthSection;
Tcl_CmdProc TclCommand_addZeroLengthContact2D;
Tcl_CmdProc TclCommand_addZeroLengthContact3D;
Tcl_CmdProc TclCommand_addZeroLengthRocking;
Tcl_CmdProc TclCommand_addZeroLengthND;
// Truss
extern OPS_Routine OPS_Truss2;
Tcl_CmdProc TclCommand_addTruss;
Tcl_CmdProc TclCommand_addTwoNodeLink;
Tcl_CmdProc TclCommand_addTwoNodeLinkSection;
// Plane
Tcl_CmdProc TclBasicBuilder_addFourNodeQuad;
Tcl_CmdProc TclBasicBuilder_addFourNodeQuadWithSensitivity;
Tcl_CmdProc TclBasicBuilder_addConstantPressureVolumeQuad;
Tcl_CmdProc TclBasicBuilder_addNineNodeMixedQuad;
Tcl_CmdProc TclCommand_SSPquadUP;
// Tcl_CmdProc TclBasicBuilder_addSixNodeTri;
Tcl_CmdProc TclBasicBuilder_addFourNodeQuadUP;
Tcl_CmdProc TclBasicBuilder_addNineFourNodeQuadUP;
Tcl_CmdProc TclBasicBuilder_addBBarFourNodeQuadUP;
// Frame
Tcl_CmdProc TclBasicBuilder_addElasticBeam;
Tcl_CmdProc TclBasicBuilder_addGradientInelasticBeamColumn;
Tcl_CmdProc TclBasicBuilder_addForceBeamColumn;
Tcl_CmdProc TclBasicBuilder_addBeamWithHinges;
Tcl_CmdProc TclBasicBuilder_addDispBeamColumnInt;
extern OPS_Routine OPS_ElasticTimoshenkoBeam2d;
extern OPS_Routine OPS_ElasticTimoshenkoBeam3d;
extern OPS_Routine OPS_AxEqDispBeamColumn2d;
extern OPS_Routine OPS_BeamGT;
extern OPS_Routine OPS_ComponentElement2d;
extern OPS_Routine OPS_ComponentElement3d;
#if defined(_HAVE_LHNMYS) || defined(OPSDEF_ELEMENT_LHNMYS)
  extern void *OPS_BeamColumn2DwLHNMYS(G3_Runtime*);
  extern void *OPS_Beam2dDamage(G3_Runtime*);
  extern void *OPS_BeamColumn2DwLHNMYS_Damage(G3_Runtime*);
  extern void *OPS_BeamColumn3DwLHNMYS(G3_Runtime*);
#endif
// Shell
Tcl_CmdProc TclBasicBuilder_addShell;
Tcl_CmdProc TclDispatch_newShellANDeS;
extern OPS_Routine OPS_ASDShellT3;
// Solid
extern OPS_Routine OPS_FourNodeTetrahedron;
extern OPS_Routine OPS_TenNodeTetrahedron;
// Brick
Tcl_CmdProc TclBasicBuilder_addBrickUP;
Tcl_CmdProc TclBasicBuilder_addBBarBrickUP;
Tcl_CmdProc TclBasicBuilder_addTwentyEightNodeBrickUP;
Tcl_CmdProc TclBasicBuilder_addTwentyNodeBrick;
Tcl_CmdProc TclBasicBuilder_addBrick;
Tcl_CmdProc TclCommand_SSP_Element;
//
Tcl_CmdProc TclCommand_addActuator;
Tcl_CmdProc TclCommand_addActuatorCorot;
Tcl_CmdProc TclCommand_addAdapter;
// Bearing
Tcl_CmdProc TclCommand_addFlatSliderBearing;
Tcl_CmdProc TclCommand_addSingleFPBearing;
Tcl_CmdProc TclBasicBuilder_addRJWatsonEqsBearing;
Tcl_CmdProc TclBasicBuilder_addYamamotoBiaxialHDR;
extern OPS_Routine OPS_ElastomericBearingPlasticity2d;
extern OPS_Routine OPS_ElastomericBearingPlasticity3d;
extern OPS_Routine OPS_ElastomericBearingBoucWen2d;
extern OPS_Routine OPS_ElastomericBearingBoucWen3d;
extern OPS_Routine OPS_ElastomericBearingUFRP2d;
extern OPS_Routine OPS_ElastomericBearingBoucWenMod3d;
extern OPS_Routine OPS_TripleFrictionPendulum;
extern OPS_Routine OPS_TripleFrictionPendulumX;
extern OPS_Routine OPS_HDR;
extern OPS_Routine OPS_LeadRubberX;
extern OPS_Routine OPS_LeadRubberY;
extern OPS_Routine OPS_ElastomericX;
// Joint
Tcl_CmdProc TclBasicBuilder_addJoint2D;
Tcl_CmdProc TclBasicBuilder_addJoint3D;
Tcl_CmdProc TclBasicBuilder_addBeamColumnJoint;
// Other
Tcl_CmdProc TclBasicBuilder_addWheelRail;
Tcl_CmdProc TclBasicBuilder_addElement2dYS;
Tcl_CmdProc TclBasicBuilder_addElastic2dGNL;
Tcl_CmdProc TclBasicBuilder_addKikuchiBearing;
Tcl_CmdProc TclBasicBuilder_addGenericCopy;
Tcl_CmdProc TclBasicBuilder_addGenericClient;

namespace OpenSees {

namespace Library {

const static
CommandLibrary 
// std::unordered_map<std::string, Tcl_CmdProc *, CaseInsensitive, CaseInsensitive> 
ElementLibrary = {
// Link
  {"twoNodeLink",               TclCommand_addTwoNodeLink},
  {"Link",                      TclCommand_addTwoNodeLink},
  {"twoNodeLinkSection",        TclCommand_addTwoNodeLinkSection},
// Truss
  {"Truss",                     TclCommand_addTruss},
  {"TrussSection",              TclCommand_addTruss},
  {"CorotTruss",                TclCommand_addTruss},
  {"CorotTrussSection",         TclCommand_addTruss},
//
// Plane
//
  {"Q4",                        TclBasicBuilder_addFourNodeQuad},
  {"stdQuad",                   TclBasicBuilder_addFourNodeQuad},
  {"LagrangeQuad",              TclBasicBuilder_addFourNodeQuad},
  {"enhancedQuad",              TclBasicBuilder_addFourNodeQuad},
  {"quad",                      TclBasicBuilder_addFourNodeQuad},
  {"SSPquad",                   TclBasicBuilder_addFourNodeQuad},
  // Q9
  {"Q9",                        TclBasicBuilder_addFourNodeQuad},
  {"quad9n",                    TclBasicBuilder_addFourNodeQuad},
  // Q8
  {"Q8",                        TclBasicBuilder_addFourNodeQuad},
  {"quad8n",                    TclBasicBuilder_addFourNodeQuad},
  // T3
  {"T3",                        TclBasicBuilder_addFourNodeQuad},
  {"tri31",                     TclBasicBuilder_addFourNodeQuad},
  {"CST",                       TclBasicBuilder_addFourNodeQuad},
  // T6
  {"T6",                        TclBasicBuilder_addFourNodeQuad},
  {"tri6n",                     TclBasicBuilder_addFourNodeQuad},
  // {"tri6n",                     TclBasicBuilder_addSixNodeTri},


  {"quadWithSensitivity",       TclBasicBuilder_addFourNodeQuadWithSensitivity},

  {"Q1/P0",                     TclBasicBuilder_addConstantPressureVolumeQuad},
  {"bbarQuad",                  TclBasicBuilder_addConstantPressureVolumeQuad},
  {"mixedQuad",                 TclBasicBuilder_addConstantPressureVolumeQuad},

  {"nineNodeMixedQuad",         TclBasicBuilder_addNineNodeMixedQuad},
  {"nineNodeQuad",              TclBasicBuilder_addNineNodeMixedQuad}, // ??

//
// Frame
//
  {"elasticBeamColumn",            TclBasicBuilder_addElasticBeam},
  {"elasticBeam",                  TclBasicBuilder_addElasticBeam},
  {"PrismFrame",                   TclBasicBuilder_addElasticBeam},
  {"ElasticFrame",                 TclBasicBuilder_addElasticBeam},

  // Nonlinear, nonstandard
  {"BeamWithHinges",               TclBasicBuilder_addBeamWithHinges},
  {"dispBeamColumnInt",            TclBasicBuilder_addDispBeamColumnInt},
  {"gradientInelasticBeamColumn",  TclBasicBuilder_addGradientInelasticBeamColumn},
  // Nonlinear
  {"DisplFrame",                   TclBasicBuilder_addForceBeamColumn},
  {"CubicFrame",                   TclBasicBuilder_addForceBeamColumn},
  {"EulerFrame",                   TclBasicBuilder_addForceBeamColumn},
  {"ForceFrame",                   TclBasicBuilder_addForceBeamColumn},
  {"MixedFrame",                   TclBasicBuilder_addForceBeamColumn},
  {"ExactFrame",                   TclBasicBuilder_addForceBeamColumn},
  {"ExactFrame02",                 TclBasicBuilder_addForceBeamColumn},
  {"CosseratFrame",                TclBasicBuilder_addForceBeamColumn},
  {"CosseratFrame01",              TclBasicBuilder_addForceBeamColumn},
  {"CosseratFrame02",              TclBasicBuilder_addForceBeamColumn},
  {"ShearFrame",                   TclBasicBuilder_addForceBeamColumn},
  {"ForceDeltaFrame",              TclBasicBuilder_addForceBeamColumn},
  {"ForceBeamColumn",              TclBasicBuilder_addForceBeamColumn},
  {"DispBeamColumn",               TclBasicBuilder_addForceBeamColumn},
  {"DispBeamColumnAsym",           TclBasicBuilder_addForceBeamColumn},
  {"TimoshenkoBeamColumn",         TclBasicBuilder_addForceBeamColumn},
  {"ForceBeamColumnCBDI",          TclBasicBuilder_addForceBeamColumn},
  {"ForceBeamColumnCSBDI",         TclBasicBuilder_addForceBeamColumn},
  {"ForceBeamColumnWarping",       TclBasicBuilder_addForceBeamColumn},
  {"ForceBeamColumnThermal",       TclBasicBuilder_addForceBeamColumn},
  {"ElasticForceBeamColumnWarping",TclBasicBuilder_addForceBeamColumn},
  {"DispBeamColumnNL",             TclBasicBuilder_addForceBeamColumn},
  {"DispBeamColumnThermal",        TclBasicBuilder_addForceBeamColumn},
  {"ElasticForceBeamColumn",       TclBasicBuilder_addForceBeamColumn},
  {"nonlinearBeamColumn",          TclBasicBuilder_addForceBeamColumn},
  {"DispBeamColumnWithSensitivity",TclBasicBuilder_addForceBeamColumn},

// Shell
  {"ThickShell01",                 TclBasicBuilder_addShell},
  {"ThickShell02",                 TclBasicBuilder_addShell},
  {"ThickShell03",                 TclBasicBuilder_addShell},
  {"ThickShell04",                 TclBasicBuilder_addShell},
  {"ThickShell05",                 TclBasicBuilder_addShell},
  {"ASDShellQ4",                   TclBasicBuilder_addShell},
  {"ShellMITC4",                   TclBasicBuilder_addShell},
  {"ShellMITC9",                   TclBasicBuilder_addShell},
  {"ShellDKGQ",                    TclBasicBuilder_addShell},
  {"ShellDKGT",                    TclBasicBuilder_addShell},
  {"ShellNLDKGQ",                  TclBasicBuilder_addShell},
  {"ShellNLDKGT",                  TclBasicBuilder_addShell},
  {"ShellMITC4Thermal",            TclBasicBuilder_addShell},
  {"ShellNLDKGQThermal",           TclBasicBuilder_addShell},
  {"ShellANDeS",                   TclDispatch_newShellANDeS},

  {"ShellQ4/F",                    TclBasicBuilder_addShell},
  {"ShellQ4/T01",                  TclBasicBuilder_addShell},
  {"ShellQ4/T02",                  TclBasicBuilder_addShell},
  {"ShellQ4/ASD",                  TclBasicBuilder_addShell},
  {"ShellQ4/E5",                   TclBasicBuilder_addShell},
  // {"ShellQ4/S",                    TclBasicBuilder_addShell},
  // {"ShellQ4/U",                    TclBasicBuilder_addShell},
  // {"ShellQ4/P0",                   TclBasicBuilder_addShell},

  {"PlateQ4/F",                    TclBasicBuilder_addShell},
  {"PlateQ4/S",                    TclBasicBuilder_addShell},
  {"PlateQ4/U",                    TclBasicBuilder_addShell},
  {"PlateQ4/L01",                  TclBasicBuilder_addShell},
  {"PlateQ4/L02",                  TclBasicBuilder_addShell},
  {"PlateQ4/E5",                   TclBasicBuilder_addShell},
  {"PlateQ4/P0",                   TclBasicBuilder_addShell},

  {"HeterosisPlate",               TclBasicBuilder_addShell},

// U-P
  {"quadUP",                    TclBasicBuilder_addFourNodeQuadUP},
  {"SSPquadUP",                 TclCommand_SSPquadUP},
  {"9_4_QuadUP",                TclBasicBuilder_addNineFourNodeQuadUP},
  {"bbarQuadUP",                TclBasicBuilder_addBBarFourNodeQuadUP},
//
// Brick
//
  {"stdBrick",                  TclBasicBuilder_addBrick},
  {"H8E12",                     TclBasicBuilder_addBrick},
  {"bbarBrick",                 TclBasicBuilder_addBrick},
  {"bbarBrickWithSensitivity",  TclBasicBuilder_addBrick},
  {"flBrick",                   TclBasicBuilder_addBrick},
  {"SSPbrick",                  TclCommand_SSP_Element},

  {"BrickUP",                   TclBasicBuilder_addBrickUP},
  {"20_8_BrickUP",              TclBasicBuilder_addTwentyEightNodeBrickUP},
  {"20NodeBrick",               TclBasicBuilder_addTwentyNodeBrick},
  {"bbarBrickUP",               TclBasicBuilder_addBBarBrickUP},

//
// Joint
//
  {"Joint2D",                   TclBasicBuilder_addJoint2D},
  {"Joint3D",                   TclBasicBuilder_addJoint3D},
  {"BeamColumnJoint",           TclBasicBuilder_addBeamColumnJoint},

// Zero-Length
  {"zeroLength",                TclCommand_addZeroLength},
  {"zeroLengthSection",         TclCommand_addZeroLengthSection},
  {"zeroLengthRocking",         TclCommand_addZeroLengthRocking},
  {"zeroLengthContact2D",       TclCommand_addZeroLengthContact2D},
  {"zeroLengthContact3D",       TclCommand_addZeroLengthContact3D},
  {"zeroLengthND",              TclCommand_addZeroLengthND},

// Actuators
  {"actuator",                  TclCommand_addActuator},
  {"corotActuator",             TclCommand_addActuatorCorot},
  {"adapter",                   TclCommand_addAdapter},

// Bearing
  {"RJWatsonEqsBearing",        TclBasicBuilder_addRJWatsonEqsBearing},
  {"RJWatsonBearing",           TclBasicBuilder_addRJWatsonEqsBearing},
  {"EQSBearing",                TclBasicBuilder_addRJWatsonEqsBearing},
  {"KikuchiBearing",            TclBasicBuilder_addKikuchiBearing},
  {"YamamotoBiaxialHDR",        TclBasicBuilder_addYamamotoBiaxialHDR},
  {"FlatSliderBearing",         TclCommand_addFlatSliderBearing},
  {"SingleFPBearing",           TclCommand_addSingleFPBearing},
  {"SinglePFBearing",           TclCommand_addSingleFPBearing},
  {"SFPBearing",                TclCommand_addSingleFPBearing},
  {"SPFBearing",                TclCommand_addSingleFPBearing},

// Other
  {"WheelRail",                 TclBasicBuilder_addWheelRail},
};

} // namespace OpenSees 
} // namespace ModelingCommands


static
std::unordered_map<std::string, OPS_Routine *, CaseInsensitive, CaseInsensitive> 
element_dispatch = {
// Truss
  {"N4BiaxialTruss",               OPS_N4BiaxialTruss},
  {"Truss2",                       OPS_Truss2},
  {"CorotTruss2",                  OPS_CorotTruss2},
  {"InertiaTruss",                 OPS_InertiaTrussElement},


// Shell
  {"ASDShellT3",                   OPS_ASDShellT3},

// Point
  {"zeroLengthContactNTS2D",       OPS_ZeroLengthContactNTS2D},
  {"zeroLengthInterface2D",        OPS_ZeroLengthInterface2D},
  {"zeroLengthImpact3D",           OPS_ZeroLengthImpact3D},

// Frame
  {"componentElement2d",           OPS_ComponentElement2d},
  {"componentElement3d",           OPS_ComponentElement3d},
#if 0
  {"componentElementDamp2d",       OPS_ComponentElementDamp2d},
#endif
  {"ModElasticBeam2d",             OPS_ModElasticBeam2d},
  {"ModElasticBeam3d",             OPS_ModElasticBeam3d},

// Solid
  {"FourNodeTetrahedron",          OPS_FourNodeTetrahedron},
  {"TenNodeTetrahedron",           OPS_TenNodeTetrahedron},

// Bearing
  {"FPBearingPTV",                 OPS_FPBearingPTV},
  {"TripleFrictionPendulum",       OPS_TripleFrictionPendulum},
  {"TripleFrictionPendulumX",      OPS_TripleFrictionPendulumX},
  {"HDR",                          OPS_HDR},
//{"LeadRubberX",                  OPS_LeadRubberX},
  {"LeadRubberX",                  OPS_LeadRubberY},
  {"ElastomericX",                 OPS_ElastomericX},

  {"AxEqDispBeamColumn2d",         OPS_AxEqDispBeamColumn2d},

#ifdef XARA_HAVE_MVLEM
// MVLEM
  {"MVLEM",                        OPS_MVLEM},        // Kristijan Kolozvari
  {"SFI_MVLEM",                    OPS_SFI_MVLEM},    // Kristijan Kolozvari
  {"MVLEM_3D",                     OPS_MVLEM_3D},     // Kristijan Kolozvari
  {"SFI_MVLEM_3D",                 OPS_SFI_MVLEM_3D}, // Kristijan Kolozvari
  {"E_SFI_MVLEM_3D",               OPS_E_SFI_MVLEM_3D},
  {"E_SFI",                        OPS_E_SFI},
  {"MEFI",                         OPS_MEFI},
#endif

// Fluid
  {"FSIFluidElement2D",            OPS_FSIFluidElement2D },
  {"FSIInterfaceElement2D",        OPS_FSIInterfaceElement2D },
  {"FluidInterface",               OPS_FSIInterfaceElement2D },
  {"FSIFluidBoundaryElement2D",    OPS_FSIFluidBoundaryElement2D },

// Joint
  {"ElasticTubularJoint",          OPS_ElasticTubularJoint},
  {"Inno3DPnPJoint",               OPS_Inno3DPnPJoint},

// Other
  {"MasonPan12",                   OPS_MasonPan12},
  {"MasonPan3D",                   OPS_MasonPan3D},
  {"BeamGT",                       OPS_BeamGT},
  {"ZeroLengthVG_HG",              OPS_ZeroLengthVG_HG},
  {"ZeroLengthContactASDimplex",   OPS_ZeroLengthContactASDimplex},
  {"SurfaceLoad",                  OPS_SurfaceLoad},
  {"TriSurfaceLoad",               OPS_TriSurfaceLoad},
  {"TPB1D",                        OPS_TPB1D},
  {"quad3d",                       OPS_FourNodeQuad3d},
  {"AC3D8",                        OPS_AC3D8HexWithSensitivity},
  {"ASI3D8",                       OPS_ASID8QuadWithSensitivity},
  {"AV3D4",                        OPS_AV3D4QuadWithSensitivity},
  {"ElastomericBearingBoucWenMod", OPS_ElastomericBearingBoucWenMod3d},
  {"VS3D4",                        OPS_VS3D4QuadWithSensitivity},
  {"CatenaryCable",                OPS_CatenaryCableElement},
  {"ASDEmbeddedNodeElement",       OPS_ASDEmbeddedNodeElement},
  {"LysmerTriangle",               OPS_LysmerTriangle},
  {"ASDAbsorbingBoundary2D",       OPS_ASDAbsorbingBoundary2D},
  {"ASDAbsorbingBoundary3D",       OPS_ASDAbsorbingBoundary3D},


  {"LinearElasticSpring",          OPS_LinearElasticSpring},
  {"Inerter",                      OPS_Inerter},
  {"Adapter",                      OPS_Adapter},
  {"Actuator",                     OPS_Actuator},
  {"CorotActuator",                OPS_ActuatorCorot},
  {"RockingBC",                    OPS_RockingBC},
  {"LehighJoint2D",                OPS_LehighJoint2d},
};

