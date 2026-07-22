#
# The interpreter core implements the infrastructure for interpreters.
#
# This includes:
# - The opaque Xara::Interpreter object that wraps the Tcl_Interp
# - Command interface (Parsing.h):
#   - ArgSize
#   - StrList
#   - ObjList


add_library(Xara_Parsing OBJECT)

target_sources(Xara_InterpreterCore PRIVATE 
  # CommandManager.cpp
  OpenSeesVersion.cpp
)

target_include_directories(Xara_InterpreterCore 
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
)

target_link_libraries(Xara_InterpreterCore 
  PUBLIC
    OPS_Recorder
    G3
    Xara_Parsing
    OPS_Logging
)

add_subdirectory(legacy)
add_subdirectory(parsing)