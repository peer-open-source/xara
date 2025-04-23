#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Actor
    PRIVATE
      ChannelAddress.cpp
      SocketAddress.cpp
    PUBLIC
      ChannelAddress.h
      SocketAddress.h
)

target_sources(OPS_Parallel
  PRIVATE
    MPI_ChannelAddress.cpp
  PUBLIC
    MPI_ChannelAddress.h
)

target_include_directories(OPS_Actor    PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_include_directories(OPS_Parallel PUBLIC ${CMAKE_CURRENT_LIST_DIR})

