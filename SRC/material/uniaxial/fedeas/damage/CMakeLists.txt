#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Runtime
  PRIVATE
    commands.cpp
    evolution.cpp
    distributions.c
    DegradingUniaxialWrapper.cpp
  PUBLIC
    DegradingUniaxialWrapper.h
)
