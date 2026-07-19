#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#

add_library(Xara_CommandInterface INTERFACE)

target_include_directories(Xara_CommandInterface INTERFACE 
  ${OPS_SRC_DIR}/runtime/commands
  ${OPS_SRC_DIR}/runtime/parsing
  ${OPS_SRC_DIR}/runtime/runtime
  ${OPS_SRC_DIR}/runtime/interpreter
  ${OPS_SRC_DIR}/logging/
  ${OPS_SRC_DIR}
)

add_library(XaraCmds_Analysis OBJECT)


target_sources(OPS_Runtime PRIVATE
# Interpreter
    "interpreter.cpp"
    "pragma.cpp"
    "packages.cpp"
    "parallel/sequential.cpp"
    # "parallel/machine.cpp"

# Modeling
    "modeling/model.cpp"
    "modeling/nodes.cpp"
    "modeling/constraint.cpp"
    
    "modeling/section.cpp"
    "modeling/section/shell.cpp"
    "modeling/section/plane.cpp"
    "modeling/section/frame.cpp"
    "modeling/section/truss.cpp"

    "modeling/printing.cpp"

    "modeling/utilities/blockND.cpp"
    "modeling/utilities/Block2D.cpp"
    "modeling/utilities/Block3D.cpp"


# Utilities
    "utilities/utilities.cpp"
    "utilities/progress.cpp"
    "utilities/formats.cpp"
)

add_subdirectory(domain)
add_subdirectory(parallel)
add_subdirectory(analysis)
add_subdirectory(modeling)

target_include_directories(OPS_Runtime PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_include_directories(OPS_Runtime PUBLIC ${CMAKE_CURRENT_LIST_DIR}/modeling)


target_link_libraries(OPS_Runtime PUBLIC XaraCmds_Analysis)
