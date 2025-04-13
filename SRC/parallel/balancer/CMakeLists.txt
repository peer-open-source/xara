#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Partition
  PRIVATE
#   ReleaseHeavierToLighterNeighbours.h
#   SwapHeavierToLighterNeighbours.h
    LoadBalancer.cpp
#   ReleaseHeavierToLighterNeighbours.cpp
    ShedHeaviest.cpp
#   SwapHeavierToLighterNeighbours.cpp
    PUBLIC
    LoadBalancer.h
    ShedHeaviest.h
)

target_include_directories(OPS_Partition PUBLIC ${CMAKE_CURRENT_LIST_DIR})
