#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Domain
    PRIVATE
        LoadPattern.cpp
        MultiSupportPattern.cpp
        UniformExcitation.cpp
        LoadCase.cpp
        StaticPattern.cpp
    PUBLIC
        LoadPattern.h
        MultiSupportPattern.h
        UniformExcitation.h
)


target_include_directories(OPS_Domain PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(drm)
add_subdirectory(series)


