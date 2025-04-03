#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Runtime PRIVATE dispatch.cpp)

target_sources(OPS_SysOfEqn
  PRIVATE
    MumpsSOE.cpp
    MumpsSolver.cpp
  PUBLIC
    MumpsSOE.h
    MumpsSolver.h

)

target_include_directories(OPS_SysOfEqn PUBLIC ${CMAKE_CURRENT_LIST_DIR})

#
# 
#
add_library(OPS_Mumps_Parallel OBJECT)

target_sources(OPS_Mumps_Parallel
  PRIVATE
    MumpsParallelSOE.cpp
    MumpsParallelSolver.cpp
)

target_include_directories(OPS_Mumps_Parallel PUBLIC ${CMAKE_CURRENT_LIST_DIR})
