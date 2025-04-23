#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Runtime PRIVATE dispatch.cpp)

target_sources(OPS_Element
    PRIVATE
        DispBeamColumn2dInt.cpp
    PUBLIC
        DispBeamColumn2dInt.h
)
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

