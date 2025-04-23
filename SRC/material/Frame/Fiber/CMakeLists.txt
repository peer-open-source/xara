#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Material
  PRIVATE
    BeamFiberMaterial.cpp
    BeamFiberMaterial2d.cpp
    BeamFiberMaterial2dPS.cpp
  PUBLIC
    BeamFiberMaterial.h
    BeamFiberMaterial2d.h
    BeamFiberMaterial2dPS.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})
