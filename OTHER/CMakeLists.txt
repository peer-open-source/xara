#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#     (c) Copyright 1999-2021 The Regents of the University of California
#                             All Rights Reserved
# (Copyright and Disclaimer @ http://www.berkeley.edu/OpenSees/copyright.html)
#
#------------------------------------------------------------------------------
# Pre-compile OpenSees numeric libraries
#==============================================================================

cmake_minimum_required(VERSION 3.16)
project(OpenSeesNumPack)

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

add_subdirectory(AMD/)
add_subdirectory(UMFPACK/)
add_subdirectory(ITPACK/)
if (FALSE)
  add_subdirectory(ARPACK/)
else()
  set(arpack_SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/arpack-ng)
  add_subdirectory(arpack-ng/)

  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(arpack PRIVATE -Wno-unused-variable
#   $<$<COMPILE_LANG_AND_ID:Fortran,GNU>:"-Wno-unused-variable">
    )
  endif()

  add_library(ARPACK INTERFACE)
  target_link_libraries(ARPACK INTERFACE arpack)
endif()

add_subdirectory(SuperLU_5.1.1/)
add_subdirectory(CSPARSE/)
add_subdirectory(METIS/)

if(MPI_FOUND)
  message("adding MUMPS and SUPERLU_DIST")
# add_subdirectory(MUMPS_5.4.1/)
  add_subdirectory(SuperLU_DIST_4.3/SRC)
endif()

