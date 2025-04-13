#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#                             External Libraries
#
#==============================================================================


# Synopsis
# - opensees_load(<PACKAGE> [BUILD|FIND|SEARCH|PATHS] [<PATHS>])
#
# Options:
# - BUILD:  Build OpenSees provided library
# - FIND:   Use CMake to find library, fail if not found
# - SEARCH: Try finding library with CMake, build OpenSees
#           Version if not found.
# - PATHS:  Provide specific paths for library.
#
#==============================================================================
find_package(LAPACK)

set(CONDA_PREFIX $ENV{CONDA_PREFIX})
# set(CMAKE_PREFIX_PATH "${CONDA_PREFIX}/lib/")

# message("CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")

# opensees_load(TCL                                          FIND
#   #LIBRARY ${CONDA_DIR}/Library/lib/tcl86t.lib
#   #INCLUDE ${CONDA_DIR}/Library/include 
# )

#set(TCL_LIBRARY ${TCL_LIBRARIES})

# opensees_load(BLAS                                         FIND
#   #LIBRARY ${CONDA_PREFIX}/Library/lib/blas.lib
#   #INCLUDE ${CONDA_PREFIX}/Library/include/
# )

# opensees_load(CBLAS                                         FIND
#   #LIBRARY ${CONDA_PREFIX}/Library/lib/cblas.lib
#   #INCLUDE ${CONDA_PREFIX}/Library/include/
# )

# opensees_load(LAPACK                                       FIND
#   #  LIBRARY ${CONDA_PREFIX}/lib/liblapack.so
#   #LIBRARY ${CONDA_PREFIX}/Library/lib/lapack.lib
# )

# set(ENV{SUPERLU_DIR})

# opensees_load(SUPERLU                                       #SEARCH
#     BUNDLED ${OPS_BUNDLED_DIR}/SuperLU_5.1.1/
# )


# opensees_load(METIS                                        SEARCH)

# opensees_load(HDF5                                           FIND)

# opensees_load(MySQL                                          FIND
#   #LIBRARY ${CONDA_PREFIX}/Library/lib/libmysql.lib
#   #INCLUDE ${CONDA_PREFIX}/Library/include/mysql
# )

# set(MYSQL_INCLUDE_DIR "${CONDA_PREFIX}/Library/include/mysql/")

