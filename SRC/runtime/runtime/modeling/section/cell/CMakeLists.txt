#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Section_Repres
    PRIVATE
        CircSectionCell.cpp
        QuadFiberCell.cpp
    PUBLIC
        FiberCell.h
        CircSectionCell.h
        QuadFiberCell.h
)
target_include_directories(OPS_Section_Repres PUBLIC ${CMAKE_CURRENT_LIST_DIR})

