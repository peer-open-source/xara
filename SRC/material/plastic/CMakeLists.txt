#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Material 
  PRIVATE
    J2Plasticity.cpp
    UVCmultiaxial.cpp
    UVCplanestress.cpp
    SimplifiedJ2.cpp
  PUBLIC
    J2Plasticity.h
    SimplifiedJ2.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})
