#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(XaraCmds_Modeling PRIVATE dispatch.cpp)

target_sources(OPS_Material PRIVATE
  BoucWenMaterial.cpp
  BWBN.cpp
  BWBF.cpp
  DegradingPinchedBW.cpp
  BoucWenOriginal.cpp
  BoucWenInfill.cpp
)
