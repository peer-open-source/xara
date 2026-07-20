#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================


add_library(OpenSeesRT SHARED)
add_library(XaraSharedMP SHARED)

add_library(XaraStaticRT STATIC)
add_library(XaraStaticMP STATIC)

add_library(Xara_Parallel INTERFACE)


target_sources(OpenSeesRT   PRIVATE "XaraInit.cpp")
target_sources(XaraSharedMP PRIVATE "XaraInit.cpp")
target_sources(XaraStaticRT PRIVATE "XaraInit.cpp")
target_sources(XaraStaticMP PRIVATE "XaraInit.cpp")

target_compile_definitions(OpenSeesRT PUBLIC USE_TCL_STUBS)
set_property(TARGET OpenSeesRT PROPERTY POSITION_INDEPENDENT_CODE 1)

# target_compile_definitions(OPS_Runtime PRIVATE "$<$<CXX_COMPILER_ID:MSVC>:_CRT_SECURE_NO_WARNINGS>")


target_include_directories(Xara_InterpreterCore PUBLIC
  ${OPS_SRC_DIR}/runtime/commands
  ${OPS_SRC_DIR}/runtime/parsing/
  ${OPS_SRC_DIR}/runtime/commands/modeling
  ${OPS_SRC_DIR}/runtime/runtime/
  ${OPS_Includes}
)

# target_link_libraries(OPS_Runtime 
# PRIVATE 
#   G3 
#   OPS_Algorithm 
#   # Xara::FastMath commented for git push
#   Xara_InterpreterCore
# )

add_subdirectory(interpreter)
add_subdirectory(commands)
add_subdirectory(runtime)


if (DEFINED XARA_VERSION)
  set_property(
    SOURCE XaraInit.cpp
    APPEND PROPERTY COMPILE_DEFINITIONS 
           "XARA_VERSION=\"${XARA_VERSION}\""
  )
endif()


#
# Tcl
#
if (${TCL_INCLUDE_PATH})
  target_include_directories(OpenSeesRT  PUBLIC ${TCL_INCLUDE_PATH})
  # target_include_directories(OPS_Runtime PUBLIC ${TCL_INCLUDE_PATH})
endif()

target_link_libraries(OpenSeesRT 
PRIVATE 
  OPS_Renderer 
  OPS_Algorithm
  XaraCmds_General
  XaraCmds_Analysis
  XaraCmds_AnalysisRT
  XaraCmds_Modeling
  Xara_InterpreterRT
  G3 
  # Xara::FastMath commented for git push
  Xara_InterpreterCore
  ${TCL_STUB_LIBRARY}
)

target_link_libraries(XaraStaticRT
PUBLIC 
  OPS_Renderer 
  OPS_Algorithm
  XaraCmds_General
  XaraCmds_Analysis
  XaraCmds_AnalysisRT
  XaraCmds_Modeling
  Xara_InterpreterRT
  G3 
  # Xara::FastMath commented for git push
  Xara_InterpreterCore
)


target_link_libraries(XaraStaticMP
PUBLIC 
  OPS_Renderer 
  OPS_Algorithm
  OPS_Parallel
  XaraCmds_General
  XaraCmds_Analysis
  XaraCmds_AnalysisMP
  XaraCmds_Modeling
  Xara_InterpreterMP
  G3 
  # Xara::FastMath commented for git push
  Xara_InterpreterCore
)

target_link_libraries(XaraSharedMP
PUBLIC 
  OPS_Renderer 
  OPS_Algorithm
  OPS_Parallel
  XaraCmds_General
  XaraCmds_Analysis
  XaraCmds_AnalysisMP
  XaraCmds_Modeling
  Xara_InterpreterMP
  G3
  # Xara::FastMath commented for git push
  Xara_InterpreterCore
)

# find_package(pybind11 CONFIG QUIET)

# if(pybind11_FOUND)
#     message(STATUS "Found pybind11 ${pybind11_VERSION}, building Python bindings")
#     add_subdirectory(python)
# else()
#     message(STATUS "pybind11 not found, skipping Python bindings")
# endif()


if (EXISTS ${CMAKE_CURRENT_LIST_DIR}/packages/UCSD_CompGeoMech/CMakeLists.txt)
  add_subdirectory(packages/UCSD_CompGeoMech)
endif ()


# Get git hash for versioning
set(GIT_HASH "unknown")
find_package(Git QUIET)
if(Git_FOUND)
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" rev-parse --short HEAD
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    OUTPUT_VARIABLE GIT_HASH
    RESULT_VARIABLE GIT_HASH_RESULT
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
  if(NOT GIT_HASH_RESULT EQUAL 0)
    set(GIT_HASH "unknown")
  endif()
  set_property(
    SOURCE OpenSeesRT.cpp
    APPEND PROPERTY COMPILE_DEFINITIONS 
           "XARA_COMMIT_HASH=\"${GIT_HASH}\""
  )
endif()

