#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(XaraCmds_Modeling PRIVATE Tcl_UWElements.cpp)

target_sources(OPS_Element
  PRIVATE
    SSPquad.cpp
    SSPquadUP.cpp
    SSPbrick.cpp
    SSPbrickUP.cpp
    Quad4FiberOverlay.cpp
    Brick8FiberOverlay.cpp

    SimpleContact2D.cpp
    SimpleContact3D.cpp
    BeamContact2D.cpp
    BeamContact2Dp.cpp
    BeamContact3D.cpp
    BeamContact3Dp.cpp
    BeamEndContact3D.cpp
    BeamEndContact3Dp.cpp


    QuadBeamEmbedContact.cpp
    PileToe3D.cpp
    EmbeddedBeamInterfaceL.cpp
    EmbeddedBeamInterfaceP.cpp
    EmbeddedEPBeamInterface.cpp
  PUBLIC
    SSPquad.h
    SSPquadUP.h
    SSPbrick.h
    SSPbrickUP.h
    Quad4FiberOverlay.h
    Brick8FiberOverlay.h
    PileToe3D.h
    EmbeddedEPBeamInterface.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
