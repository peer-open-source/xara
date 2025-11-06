#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(G3_ObjectBroker
  PRIVATE
    FEM_ObjectBroker.cpp
  PUBLIC
    FEM_ObjectBroker.h
)

target_include_directories(G3_ObjectBroker
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}
)