#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# _PARALLEL_INTERPRETERS
# _PARALLEL_PROCESSING     # OpenSeesSP

return()


if (NOT MPI_FOUND)
  return()
endif()


add_library(OpenSeesRT_Parallel OBJECT EXCLUDE_FROM_ALL)

target_include_directories(OpenSeesRT_Parallel
  PRIVATE
  ${OPS_SRC_DIR}/domain/partitioner/
  ${OPS_SRC_DIR}/analysis/analysis/
)


target_link_libraries(OpenSeesRT_Parallel PRIVATE 
    OPS_Runtime 
    OPS_Parallel 
    OPS_Actor 
    MPI::MPI_CXX
#   ${MPI_CXX_LIBRARIES}
)

target_sources(OpenSeesRT_Parallel PRIVATE
  machine.cpp
)

if (MKL_FOUND)
  message(STATUS "MKL was found.")
  set(SCALAPACK_LIBRARIES ${MKL_LIBRARIES})
else()
  message(STATUS "MKL NOT found.")
  find_package(SCALAPACK)
  message(">>> ScaLAPACK: ${ScaLAPACK_LIBRARIES}")
endif()
  
target_link_libraries(OpenSeesRT_Parallel PRIVATE 
    OPS_Parallel 
    OPS_Actor
    MPI::MPI_CXX
#   /usr/lib/libscalapack.so
    ${MPI_Fortran_LIBRARIES}
)
