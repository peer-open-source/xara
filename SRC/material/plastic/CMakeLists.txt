#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Material 
  PRIVATE
    J2Plasticity.cpp
    J2ThreeDimensional.cpp
    J2BeamFiber2d.cpp
    J2BeamFiber3d.cpp
    UVCmultiaxial.cpp
    UVCplanestress.cpp
    SimplifiedJ2.cpp
  PUBLIC
    J2Plasticity.h
    SimplifiedJ2.h
    J2BeamFiber2d.h
    J2BeamFiber3d.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})
