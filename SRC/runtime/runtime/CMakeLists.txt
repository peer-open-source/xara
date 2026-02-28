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
      analysis/EigenSolve.cpp
      TclPackageClassBroker.cpp

    PUBLIC
      BasicAnalysisBuilder.h
      ModelRegistry.h
      TclPackageClassBroker.h
)

add_subdirectory(modeling/section)
add_subdirectory(legacy/Renderer)


