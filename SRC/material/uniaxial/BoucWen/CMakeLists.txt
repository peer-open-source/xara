#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Runtime PRIVATE dispatch.cpp)

target_sources(OPS_Material PRIVATE
  BoucWenMaterial.cpp
  BWBN.cpp
  DegradingPinchedBW.cpp
  BoucWenOriginal.cpp
  BoucWenInfill.cpp
)
