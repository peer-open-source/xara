#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
  # target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_sources(OPS_Runtime PRIVATE dispatch.cpp)
target_sources(OPS_Element
    PRIVATE
        GradientInelasticBeamColumn2d.cpp
        GradientInelasticBeamColumn3d.cpp
        #TclGradientInelasticBeamColumnCommand.cpp
    PUBLIC
        GradientInelasticBeamColumn2d.h
        GradientInelasticBeamColumn3d.h
)

