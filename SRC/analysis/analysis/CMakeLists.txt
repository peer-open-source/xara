#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Analysis
    PRIVATE
      Analysis.cpp

      DomainDecompositionAnalysis.cpp
      StaticDomainDecompositionAnalysis.cpp 
      TransientDomainDecompositionAnalysis.cpp 
    PUBLIC
      Analysis.h

      DomainDecompositionAnalysis.h
      StaticDomainDecompositionAnalysis.h 
      TransientDomainDecompositionAnalysis.h 
)
