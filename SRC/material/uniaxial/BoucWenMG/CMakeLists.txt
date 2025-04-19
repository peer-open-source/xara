#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_include_directories(OPS_Runtime PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_sources(OPS_Runtime PRIVATE commands.cpp)

target_sources(OPS_Material
  PRIVATE
    BoucWenMG.cpp
  PUBLIC
    BoucWenMG.h
)

