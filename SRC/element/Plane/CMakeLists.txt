#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Element
  PRIVATE
    IGAQuad.cpp
    ConstantPressureVolumeQuad.cpp
    EnhancedQuad.cpp
    FourNodeQuad.cpp
    FourNodeQuad3d.cpp
    FourNodeQuadWithSensitivity.cpp
    EightNodeQuad.cpp
    NineNodeMixedQuad.cpp
    NineNodeQuad.cpp
    SixNodeTri.cpp
    BBarFourNodeQuadUP.cpp
    FourNodeQuadUP.cpp
    Nine_Four_Node_QuadUP.cpp
    Tri31.cpp
  PUBLIC
    ConstantPressureVolumeQuad.h
    EnhancedQuad.h
    FourNodeQuad.h
    FourNodeQuad3d.h
    FourNodeQuadWithSensitivity.h
    EightNodeQuad.h
    NineNodeMixedQuad.h
    NineNodeQuad.h
    SixNodeTri.h
    BBarFourNodeQuadUP.h
    FourNodeQuadUP.h
    Nine_Four_Node_QuadUP.h
    Tri31.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
