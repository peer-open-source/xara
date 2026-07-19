#=====================================================================
# 
#   OpenSees -- Open System For Earthquake Engineering Simulation
#           Pacific Earthquake Engineering Research Center
#
#=====================================================================

target_sources(OPS_Runtime
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

if (TRUE)
  target_sources(OPS_Runtime PRIVATE ProcessContext.cpp)
  target_link_libraries(OPS_Runtime PUBLIC OPS_Parallel)
endif()

add_subdirectory(modeling/section)
add_subdirectory(legacy/Renderer)


