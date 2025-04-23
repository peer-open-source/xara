#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
add_library(graph OBJECT)
target_link_libraries(graph PRIVATE OPS_Logging OPS_Analysis)

add_subdirectory(graph)
add_subdirectory(numberer)
add_subdirectory(partitioner)
