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
# add_subdirectory(parallel)
add_subdirectory(analysis)
add_subdirectory(modeling)

# target_include_directories(OPS_Runtime PUBLIC ${CMAKE_CURRENT_LIST_DIR})
# target_include_directories(OPS_Runtime PUBLIC ${CMAKE_CURRENT_LIST_DIR}/modeling)


target_link_libraries(XaraCmdsRT 
PRIVATE 
  Xara_InterpreterCore
PUBLIC 
  XaraCmds_General
  XaraCmds_Modeling
  XaraCmds_Analysis
  XaraCmds_AnalysisRT
)

target_link_libraries(XaraCmdsMP
PRIVATE 
  Xara_InterpreterCore
PUBLIC 
  XaraCmds_General
  XaraCmds_Modeling
  XaraCmds_Analysis
  XaraCmds_AnalysisMP
)
