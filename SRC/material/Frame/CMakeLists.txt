#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Material
  PRIVATE
    FrameFiberSection3d.cpp
    FrameSolidSection3d.cpp
    ElasticLinearFrameSection3d.cpp
  PUBLIC
    FrameFiberSection3d.h
    ElasticLinearFrameSection3d.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
