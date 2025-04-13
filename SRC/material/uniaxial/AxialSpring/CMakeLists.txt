#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Runtime PRIVATE dispatch.cpp)

target_sources(OPS_Material
  PRIVATE
    AxialSp.cpp
    AxialSpHD.cpp
  PUBLIC
    AxialSp.h
    AxialSpHD.h
)
