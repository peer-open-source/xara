#
# The interpreter core implements the infrastructure for interpreters.
#
# This includes:
# - The opaque Xara::Interpreter object that wraps the Tcl_Interp
# - Command interface (Parsing.h):
#   - ArgSize
#   - StrList
#   - ObjList

add_library(Xara_InterpreterCore OBJECT)

target_sources(Xara_InterpreterCore PRIVATE 
  CommandManager.cpp
  OpenSeesVersion.cpp
)

target_include_directories(Xara_InterpreterCore 
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
)

target_link_libraries(Xara_InterpreterCore 
  PUBLIC 
    Xara_Parsing
    OPS_Logging
)