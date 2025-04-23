#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Analysis
    PRIVATE
      SolutionAlgorithm.cpp
    PUBLIC
      SolutionAlgorithm.h
)

add_subdirectory(domainDecompAlgo)
add_subdirectory(equiSolnAlgo)
add_subdirectory(eigenAlgo)
