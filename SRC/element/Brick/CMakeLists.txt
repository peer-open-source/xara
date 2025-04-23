#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Element
  PRIVATE
    Brick.cpp
    Twenty_Node_Brick.cpp
    BrickUP.cpp
    BbarBrick.cpp
    BBarBrickUP.cpp
    BbarBrickWithSensitivity.cpp
    Twenty_Eight_Node_BrickUP.cpp
  PUBLIC
    Brick.h
    Twenty_Node_Brick.h
    BrickUP.h
    BbarBrick.h
    BBarBrickUP.h
    BbarBrickWithSensitivity.h
    Twenty_Eight_Node_BrickUP.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
