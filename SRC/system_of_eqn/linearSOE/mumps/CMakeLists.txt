#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
find_package(MPI)
#find_package(Mumps)
target_sources(OPS_Runtime PRIVATE dispatch.cpp)
target_compile_definitions(OPS_Runtime PRIVATE XARA_USE_MUMPS)

target_sources(OPS_SysOfEqn
  PRIVATE
    MumpsSOE.cpp
    MumpsSolver.cpp
  PUBLIC
    MumpsSOE.h
    MumpsSolver.h
)

target_link_libraries(OPS_SysOfEqn PUBLIC MPI::MPI_C MUMPS)
target_include_directories(OPS_SysOfEqn PUBLIC 
  ${CMAKE_CURRENT_LIST_DIR}
  ${CMAKE_CURRENT_LIST_DIR}/include
)

#
# 
#
if (FALSE)
add_library(OPS_Mumps_Parallel OBJECT)

target_sources(OPS_Mumps_Parallel
  PRIVATE
    MumpsParallelSOE.cpp
    MumpsParallelSolver.cpp
)

target_include_directories(OPS_Mumps_Parallel PUBLIC ${CMAKE_CURRENT_LIST_DIR})
endif()