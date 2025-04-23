#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
add_library(OPS_Reliability OBJECT)

target_include_directories(OPS_Reliability PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(./analysis/)
add_subdirectory(./domain/)
add_subdirectory(./FEsensitivity/)
add_subdirectory(./tcl/)

  # ${OPS_SRC_DIR}/reliability/analysis/analysis
  # ${OPS_SRC_DIR}/reliability/analysis/analysis/system
  # ${OPS_SRC_DIR}/reliability/analysis/convergenceCheck
  # ${OPS_SRC_DIR}/reliability/analysis/curvature
  # ${OPS_SRC_DIR}/reliability/analysis/designPoint
  # ${OPS_SRC_DIR}/reliability/analysis/direction
  # ${OPS_SRC_DIR}/reliability/analysis/gFunction
  # ${OPS_SRC_DIR}/reliability/analysis/gradient
  # ${OPS_SRC_DIR}/reliability/analysis/hessian
  # ${OPS_SRC_DIR}/reliability/analysis/hessianApproximation
  # ${OPS_SRC_DIR}/reliability/analysis/meritFunction
  # ${OPS_SRC_DIR}/reliability/analysis/misc
  # ${OPS_SRC_DIR}/reliability/analysis/randomNumber
  # ${OPS_SRC_DIR}/reliability/analysis/rootFinding
  # ${OPS_SRC_DIR}/reliability/analysis/sensitivity
  # ${OPS_SRC_DIR}/reliability/analysis/stepSize
  # ${OPS_SRC_DIR}/reliability/analysis/transformation
  # ${OPS_SRC_DIR}/reliability/domain/components
  # ${OPS_SRC_DIR}/reliability/domain/distributions
  # ${OPS_SRC_DIR}/reliability/domain/filter
  # ${OPS_SRC_DIR}/reliability/domain/functionEvaluator
  # ${OPS_SRC_DIR}/reliability/domain/modulatingFunction
  # ${OPS_SRC_DIR}/reliability/domain/performanceFunction
  # ${OPS_SRC_DIR}/reliability/domain/spectrum
  # ${OPS_SRC_DIR}/reliability/domain/storage
  # ${OPS_SRC_DIR}/reliability/FEsensitivity
  # ${OPS_SRC_DIR}/reliability/tcl
