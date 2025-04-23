#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Analysis
    PRIVATE
        ArcLength1.cpp
        ArcLength.cpp
        DisplacementControl.cpp
        DistributedDisplacementControl.cpp
        EQPath.cpp
        HarmonicSteadyState.cpp
        HSConstraint.cpp
        LoadControl.cpp
        LoadPath.cpp
        DisplacementPath.cpp
        MinUnbalDispNorm.cpp
    PUBLIC
        ArcLength1.h
        ArcLength.h
        DisplacementControl.h
        DistributedDisplacementControl.h
        EQPath.h
        HarmonicSteadyState.h
        HSConstraint.h
        LoadControl.h
        LoadPath.h
        DisplacementPath.h
        MinUnbalDispNorm.h
)
target_include_directories(OPS_Analysis PUBLIC ${CMAKE_CURRENT_LIST_DIR})

