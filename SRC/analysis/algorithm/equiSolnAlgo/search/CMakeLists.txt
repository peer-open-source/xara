#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Algorithm
    PRIVATE
      LineSearch.cpp 
      SecantLineSearch.cpp 
      RegulaFalsiLineSearch.cpp 
      BisectionLineSearch.cpp
      InitialInterpolatedLineSearch.cpp 
    PUBLIC
      LineSearch.h 
      InitialInterpolatedLineSearch.h 
      SecantLineSearch.h 
      RegulaFalsiLineSearch.h 
      BisectionLineSearch.h
)
target_include_directories(OPS_Algorithm PUBLIC ${CMAKE_CURRENT_LIST_DIR})
