#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Transform
    PRIVATE
        CorotFrameTransf3d.cpp
        CorotFrameTransf3d03.cpp
        LinearFrameTransf3d.cpp
        BasicFrameTransf.cpp
        PDeltaFrameTransf3d.cpp
    PUBLIC
        CorotFrameTransf3d.h
        CorotFrameTransf3d03.h
        LinearFrameTransf3d.h
        PDeltaFrameTransf3d.h
)
target_include_directories(OPS_Transform PUBLIC ${CMAKE_CURRENT_LIST_DIR})

