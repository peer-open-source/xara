#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#
# _PARALLEL_PROCESSING: OpenSeesSP
# _PARALLEL_INTERPRETERS: OpenSeesMP

find_package(MPI)

if (NOT MPI_FOUND)
  return()
endif()

add_subdirectory(balancer)

target_link_libraries(OPS_Parallel 
  PUBLIC 
  # OPS_Partition 
  OPS_Matrix
  MPI::MPI_CXX
)

# return()


target_include_directories(OPS_Parallel PRIVATE 
  ${OPS_Includes}
)

