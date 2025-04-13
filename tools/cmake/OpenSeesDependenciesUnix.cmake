#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#                            External Libraries
# - BLAS_LIBRARIES
# - BLAS_INCLUDE_DIRS
#
# - LAPACK_LIBRARIES
# - LAPACK_INCLUDE_DIRS
#
# - ARPACK_LIBRARIES
#
# - SUPERLU_LIBRARIES
# - SUPERLU_INCLUDE_DIRS
#
#==============================================================================
# Synopsis
# - opensees_load(<PACKAGE> [FLAGS] [<PATHS>])
#
# Flags:
# - FIND:   Use CMake to find library on system, fail if not found
#
# Keyword arguments:
#   Provide specific paths for library.
#
# - BUNDLED <path/to/OTHER/LIB/>
#   Provide path to OpenSees bundled library directory containing a 
#   CMakeLists.txt file.
#
# - LIBRARY <path/to/lib.a> INCLUDE <path/to/include/>
#
# - CONAN <conan-package/version>
#   Point to a conan package.
#
#----------------------------------------------------------------

set(CMAKE_PREFIX_PATH "/usr/lib")
include(OpenSeesFunctions)

#find_package(MPI)

opensees_load(BLAS                                           FIND)

opensees_load(LAPACK                                         FIND)

opensees_load(METIS                                          FIND)

opensees_load(HDF5                                           FIND)

opensees_load(MySQL                                          FIND)

# Integrated exteral libraries
opensees_load(FEDEAS_Uniaxial LIBRARY FALSE)


