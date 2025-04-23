#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Partition
  PRIVATE
    DomainPartitioner.cpp
  PUBLIC
    DomainPartitioner.h
)

target_include_directories(OPS_Partition PUBLIC ${CMAKE_CURRENT_LIST_DIR})
