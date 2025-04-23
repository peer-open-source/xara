#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Domain
    PRIVATE
        EarthquakePattern.cpp
        LoadPattern.cpp
        LoadPatternIter.cpp
        MultiSupportPattern.cpp
        UniformExcitation.cpp
    PUBLIC
        EarthquakePattern.h
        LoadPattern.h
        LoadPatternIter.h
        MultiSupportPattern.h
        UniformExcitation.h
)

target_sources(OPS_Domain
    PRIVATE
        FireLoadPattern.cpp
    PUBLIC
        FireLoadPattern.h
)


target_include_directories(OPS_Domain PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(drm)
add_subdirectory(series)


