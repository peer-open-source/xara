#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_link_libraries(OPS_SysOfEqn PRIVATE OPS_Logging)
target_include_directories(OPS_SysOfEqn 
  PUBLIC 
    ${CMAKE_CURRENT_LIST_DIR}
    ${OPS_SRC_DIR}/analysis/model
    ${OPS_SRC_DIR}/system_of_eqn/linearSOE
    ${OPS_SRC_DIR}/system_of_eqn/eigenSOE
    ${OPS_SRC_DIR}/system_of_eqn/linearSOE/bandGEN
    ${OPS_SRC_DIR}/system_of_eqn/linearSOE/bandSPD
    ${OPS_SRC_DIR}/system_of_eqn/linearSOE/cg
    ${OPS_SRC_DIR}/system_of_eqn/linearSOE/diagonal
    ${OPS_SRC_DIR}/system_of_eqn/linearSOE/fullGEN
    ${OPS_SRC_DIR}/system_of_eqn/linearSOE/itpack
    ${OPS_SRC_DIR}/system_of_eqn/linearSOE/mumps
    ${OPS_SRC_DIR}/system_of_eqn/linearSOE/pardiso
    ${OPS_SRC_DIR}/system_of_eqn/linearSOE/petsc
    ${OPS_SRC_DIR}/system_of_eqn/linearSOE/profileSPD
    ${OPS_SRC_DIR}/system_of_eqn/linearSOE/sparseGEN
    ${OPS_SRC_DIR}/system_of_eqn/linearSOE/sparseSYM
    ${OPS_SRC_DIR}/system_of_eqn/linearSOE/umfGEN
)

add_subdirectory(linearSOE)
add_subdirectory(eigenSOE)

