#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

add_subdirectory(algorithm)
add_subdirectory(numberer)
add_subdirectory(dof_grp)
add_subdirectory(fe_ele)
add_subdirectory(handler)
add_subdirectory(model)
add_subdirectory(integrator)
add_subdirectory(criteria)
add_subdirectory(analysis)


target_link_libraries(OPS_Analysis 
  PRIVATE 
    OPS_ConvergenceTest
  PUBLIC
    OPS_Algorithm
    OPS_Domain
    OPS_Logging
)

# Includes
target_include_directories(OPS_Analysis
  PUBLIC 
  criteria
  integragor
  model
)
target_include_directories(OPS_Algorithm
  PUBLIC 
  criteria
  integragor
  model
)
target_include_directories(OPS_ConvergenceTest
  PUBLIC 
  criteria
  integragor
  model
  algorithm/equiSolnAlgo
)

