#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================


#set(USE_STATIC_MKL TRUE)
set(BLA_STATIC ON)
set(BLA_VENDOR Intel10_64lp)
# set(BLA_VENDOR Intel10_64ilp)
set (MKL_LPATH ${MKL_ROOT}/lib/intel64)
find_package(MKL  CONFIG REQUIRED)
find_package(BLAS)
find_package(LAPACK)


# set(BUNDLE_LIBS "${PROJECT_SOURCE_DIR}/Win64/lib/debug/")
# opensees_load(UMFPACK
#   #LIBRARY "${BUNDLE_LIBS}/umfpackC.lib"
#   LIBRARY ${OPS_BUNDLED_DIR}/bin/UMFPACK/Debug/UMFPACK.lib
# )
# 

# #set(ENV{SUPERLU_DIR})
# opensees_load(SUPERLU
#   #BUNDLED ${OPS_BUNDLED_DIR}/SuperLU_5.1.1/
#   LIBRARY ${OPS_BUNDLED_DIR}/bin/SuperLU_5.1.1/Debug/SUPERLU.lib
# )
# 

