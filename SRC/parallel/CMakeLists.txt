#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
find_package(MPI)

if (NOT MPI_FOUND)
  return()
endif()

add_subdirectory(balancer)

target_link_libraries(OPS_Parallel 
    PUBLIC 
    OPS_Partition 
    MPI::MPI_CXX
)

return()

#
# OpenSeesSP Tcl Parallel Interpreter
#

add_executable(OpenSeesSP EXCLUDE_FROM_ALL 
  ${OPS_SRC_DIR}/parallel/OpenSeesSP.cpp

  OPS_Partition
  ${OPS_SRC_DIR}/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp

  # OPS_Parallel
  ${OPS_SRC_DIR}/actor/machineBroker/MPI_MachineBroker.cpp

  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/diagonal/MPIDiagonalSOE.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/diagonal/MPIDiagonalSolver.cpp

  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsSOE.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsSolver.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsParallelSOE.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsParallelSolver.cpp

  ${OPS_SRC_DIR}/domain/subdomain/ActorSubdomain.cpp
  ${OPS_SRC_DIR}/domain/subdomain/ShadowSubdomain.cpp

  OPS_SuperLU_Dist
)


target_include_directories(OpenSeesSP PRIVATE 
	${MUMPS_DIR}/_deps/mumps-src/include ${MPI_CXX_INCLUDE_DIRS})

target_compile_options(OpenSeesSP PRIVATE ${MPI_CXX_COMPILE_FLAGS})

if (DEFINED OPENMPI)
    target_compile_definitions(OpenSeesSP 
    PUBLIC _PARALLEL_PROCESSING ${MUMPS_FLAG} _OPENMPI)
else()
    target_compile_definitions(OpenSeesSP 
    PUBLIC _PARALLEL_PROCESSING ${MUMPS_FLAG})
endif()

target_link_libraries(OpenSeesSP 
   OPS_InterpTcl 
   OpenSeesLIB
   OPS_Reliability
   OPS_ReliabilityTcl
   OPS_Recorder
   METIS
   SUPERLU_DIST
   OPS_Numerics
   ${MUMPS_LIBRARIES}
   ${CMAKE_DL_LIBS} 
   ${HDF5_LIBRARIES} 
   ${CONAN_LIBS}
   ${MPI_CXX_LIBRARIES}
   ${SCALAPACK_LIBRARIES}
   ${MPI_Fortran_LIBRARIES}       
   ${MPI_CXX_LINK_FLAGS}
)




#
# OpenSeesMP
#

add_executable(OpenSeesMP EXCLUDE_FROM_ALL 
  ${OPS_SRC_DIR}/tcl/OpenSeesMP.cpp

  ${OPS_SRC_DIR}/domain/domain/partitioned/PartitionedDomain.cpp
  ${OPS_SRC_DIR}/domain/domain/partitioned/PartitionedDomainEleIter.cpp
  ${OPS_SRC_DIR}/domain/domain/partitioned/PartitionedDomainSubIter.cpp
  ${OPS_SRC_DIR}/actor/machineBroker/MPI_MachineBroker.cpp
  ${OPS_SRC_DIR}/actor/objectBroker/FEM_ObjectBrokerAllClasses.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/diagonal/MPIDiagonalSOE.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/diagonal/MPIDiagonalSolver.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsSOE.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsSolver.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsParallelSOE.cpp
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps/MumpsParallelSolver.cpp
  ${OPS_SRC_DIR}/domain/subdomain/ActorSubdomain.cpp
  ${OPS_SRC_DIR}/domain/subdomain/ShadowSubdomain.cpp
)

target_include_directories(OpenSeesMP PRIVATE 
	 ${MUMPS_DIR}/_deps/mumps-src/include ${MPI_CXX_INCLUDE_DIRS})

target_compile_options(OpenSeesMP PRIVATE ${MPI_CXX_COMPILE_FLAGS})

if (NOT DEFINED MUMPS_DIR)
    add_dependencies(OpenSeesMP mumps)
endif()

if (DEFINED OPENMPI)
    target_compile_definitions(OpenSeesMP 
    PUBLIC _PARALLEL_INTERPRETERS ${MUMPS_FLAG} _OPENMPI)
else()
    target_compile_definitions(OpenSeesMP 
    PUBLIC _PARALLEL_INTERPRETERS ${MUMPS_FLAG})
endif()

target_link_libraries(OpenSeesMP 
   OPS_InterpTcl 
   OpenSeesLIB
   OPS_Reliability
   OPS_ReliabilityTcl
   OPS_Recorder      
   METIS
   SUPERLU_DIST
   OPS_Numerics
   ${MUMPS_LIBRARIES}
   ${CMAKE_DL_LIBS} 
   ${HDF5_LIBRARIES} 
   ${CONAN_LIBS}
   ${MPI_CXX_LIBRARIES}
   ${SCALAPACK_LIBRARIES}
   ${MPI_Fortran_LIBRARIES}       
   ${MPI_CXX_LINK_FLAGS}
)

