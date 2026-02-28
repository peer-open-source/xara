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
        MultiSupportPattern.cpp
        UniformExcitation.cpp
        LoadCase.cpp
        StaticPattern.cpp
    PUBLIC
        EarthquakePattern.h
        LoadPattern.h
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


