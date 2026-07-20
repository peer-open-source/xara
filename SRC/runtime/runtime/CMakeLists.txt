#=====================================================================
# 
#   OpenSees -- Open System For Earthquake Engineering Simulation
#           Pacific Earthquake Engineering Research Center
#
#=====================================================================

target_sources(Xara_InterpreterCore
  PRIVATE
    BasicAnalysisBuilder.cpp
    ModelRegistry.cpp
    InterpreterResponse.cpp
    # analysis/EigenSolve.cpp
    XaraClassBroker.cpp

  PUBLIC
    BasicAnalysisBuilder.h
    ModelRegistry.h
    XaraClassBroker.h
)

target_sources(Xara_InterpreterRT PRIVATE ProcessContext.cpp)
target_sources(Xara_InterpreterMP PRIVATE ProcessContext.cpp)

foreach (target Xara_InterpreterRT Xara_InterpreterMP)
  target_link_libraries(${target} PUBLIC Xara_InterpreterCore)
endforeach()

if (MPI_FOUND)
  target_compile_definitions(Xara_InterpreterMP PUBLIC XARA_ENABLE_MPI)
  target_link_libraries(Xara_InterpreterMP PUBLIC 
    OPS_Parallel 
    MPI::MPI_CXX
  )
endif()

# OPS_Section_Repres
add_subdirectory(modeling/section)
# OPS_Renderer
add_subdirectory(legacy/Renderer)


