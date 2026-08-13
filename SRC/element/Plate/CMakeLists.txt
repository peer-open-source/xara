#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_link_libraries(XaraCmds_Modeling PRIVATE OPS_Element)

target_sources(OPS_Element
    PRIVATE
      HeterosisPlate.cpp
)
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

