#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Utilities
  PRIVATE
    Timer.cpp 
  PUBLIC
    Timer.h 
)

target_include_directories(OPS_Utilities PUBLIC ${CMAKE_CURRENT_LIST_DIR})
