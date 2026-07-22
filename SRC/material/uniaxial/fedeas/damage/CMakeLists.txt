#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(XaraCmds_Modeling
  PRIVATE
    commands.cpp
    evolution.cpp
    distributions.c
    DegradingUniaxialWrapper.cpp
  PUBLIC
    DegradingUniaxialWrapper.h
)
