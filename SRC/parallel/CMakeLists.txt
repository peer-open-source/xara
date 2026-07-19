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



#
# OpenSeesMP
#

add_executable(OpenSeesMP EXCLUDE_FROM_ALL 
  ${OPS_SRC_DIR}/parallel/OpenSeesMP.cpp
  ${OPS_SRC_DIR}/domain/subdomain/ActorSubdomain.cpp
  ${OPS_SRC_DIR}/domain/subdomain/ShadowSubdomain.cpp
)

target_include_directories(OPS_Parallel PRIVATE 
  ${OPS_Includes}
)

target_include_directories(OpenSeesMP PRIVATE 
  ${OPS_Includes}
  ${MUMPS_DIR}/_deps/mumps-src/include 
  ${MPI_CXX_INCLUDE_DIRS}
)

target_compile_options(OpenSeesMP PRIVATE ${MPI_CXX_COMPILE_FLAGS})

# if (NOT DEFINED MUMPS_DIR)
#     add_dependencies(OpenSeesMP mumps)
# endif()

if (DEFINED OPENMPI)
    target_compile_definitions(OpenSeesMP 
    PUBLIC _PARALLEL_INTERPRETERS ${MUMPS_FLAG} _OPENMPI)
else()
    target_compile_definitions(OpenSeesMP 
    PUBLIC _PARALLEL_INTERPRETERS ${MUMPS_FLAG})
endif()

target_link_libraries(OpenSeesMP 
  OpenSeesRT_Static
  OPS_InterpTcl 
  OpenSeesLIB
  OPS_Reliability
  OPS_ReliabilityTcl
  OPS_Recorder      
  METIS
  SUPERLU_DIST
  OPS_Numerics
  OPS_Mumps_Parallel
  ${MUMPS_LIBRARIES}
  ${CMAKE_DL_LIBS} 
  # ${HDF5_LIBRARIES} 
  ${CONAN_LIBS}
  ${MPI_CXX_LIBRARIES}
  ${SCALAPACK_LIBRARIES}
  ${MPI_Fortran_LIBRARIES}
  ${MPI_CXX_LINK_FLAGS}
)



