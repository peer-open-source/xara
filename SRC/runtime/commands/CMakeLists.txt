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



target_sources(XaraCmds_General PRIVATE
# Interpreter
    "general.cpp"
    "pragma.cpp"
    "packages.cpp"
    "parallel/sequential.cpp"
    # "parallel/machine.cpp"

# Utilities
    "utilities/utilities.cpp"
    "utilities/progress.cpp"
    "utilities/formats.cpp"
)

add_subdirectory(domain)
add_subdirectory(general)
add_subdirectory(analysis)
add_subdirectory(modeling)

