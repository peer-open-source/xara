#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_SysOfEqn
  PRIVATE
    PARDISOGenLinSOE.cpp
    PARDISOGenLinSolver.cpp
#   PARDISOSymLinSOE.cpp
#   PARDISOSymLinSolver.cpp
)

target_include_directories(OPS_SysOfEqn PUBLIC ${CMAKE_CURRENT_LIST_DIR})

