# intel_fortran_toolchain.cmake

set(CMAKE_SYSTEM_NAME Windows)
set(CMAKE_GENERATOR "Ninja")

set(CONDA_PREFIX $ENV{CONDA_PREFIX})
# Fortran compiler
set(CMAKE_Fortran_COMPILER "${CONDA_PREFIX}/Library/bin/ifx.exe")
set(CMAKE_Fortran_PREPROCESS ON)

# Compiler flags
set(CMAKE_Fortran_FLAGS "/fpp /libs:static")
set(CMAKE_Fortran_FLAGS_INIT "-IC:${CONDA_PREFIX}/opt/compiler/include/intel64/")
add_link_options("/NODEFAULTLIB:ifconsol.lib")

# Library paths
set(CMAKE_PREFIX_PATH "${CONDA_PREFIX}/Library/lib/")
set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)

#set(Dependencies "Conda")


# BLAS and LAPACK
add_link_options("/LIBPATH:${CONDA_PREFIX}/Library/lib")
set(BLAS_LIBRARIES "$ENV{CONDA_PREFIX}/Library/lib/mkl_rt.lib" CACHE STRING "BLAS library")
set(LAPACK_LIBRARIES "$ENV{CONDA_PREFIX}/Library/lib/mkl_rt.lib" CACHE STRING "LAPACK library")
set(BLAS_FOUND TRUE CACHE BOOL "Force BLAS to be found")
set(LAPACK_FOUND TRUE CACHE BOOL "Force LAPACK to be found")

