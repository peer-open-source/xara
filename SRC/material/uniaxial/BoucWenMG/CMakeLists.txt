#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_include_directories(XaraCmds_Modeling PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_sources(XaraCmds_Modeling PRIVATE commands.cpp)

target_sources(OPS_Material
  PRIVATE
    BoucWenMG.cpp
  PUBLIC
    BoucWenMG.h
)

