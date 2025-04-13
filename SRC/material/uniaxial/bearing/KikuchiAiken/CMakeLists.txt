#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Material
  PRIVATE
      KikuchiAikenHDR.cpp
      KikuchiAikenLRB.cpp
  PUBLIC
      KikuchiAikenHDR.h
      KikuchiAikenLRB.h
)

target_sources(OPS_Runtime PRIVATE dispatch.cpp)
