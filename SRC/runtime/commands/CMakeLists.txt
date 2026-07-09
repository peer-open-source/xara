#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#

target_sources(OPS_Runtime PRIVATE
# Interpreter
    "interpreter.cpp"
    "pragma.cpp"
    "packages.cpp"
    "parallel/sequential.cpp"

# Modeling
    "modeling/model.cpp"
    "modeling/nodes.cpp"
    "modeling/constraint.cpp"
    "modeling/geomTransf.cpp"
    
    "modeling/section.cpp"
    "modeling/section/shell.cpp"
    "modeling/section/plane.cpp"
    "modeling/section/frame.cpp"
    "modeling/section/truss.cpp"

    "modeling/printing.cpp"

    "modeling/utilities/blockND.cpp"
    "modeling/utilities/Block2D.cpp"
    "modeling/utilities/Block3D.cpp"

# Analysis
    "analysis/algorithm.cpp"
    "analysis/integrator.cpp"
    "analysis/transient.cpp"
    "analysis/analysis.cpp"
    "analysis/numberer.cpp"
    "analysis/ctest.cpp"
    "analysis/solver.cpp"
    "analysis/solver.hpp"
    "analysis/sensitivity.cpp"

# Utilities
    "utilities/utilities.cpp"
    "utilities/progress.cpp"
    "utilities/formats.cpp"
)

add_subdirectory(domain)
add_subdirectory(parallel)
add_subdirectory(analysis/modal)
add_subdirectory(modeling)

target_include_directories(OPS_Runtime PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_include_directories(OPS_Runtime PUBLIC ${CMAKE_CURRENT_LIST_DIR}/modeling)
