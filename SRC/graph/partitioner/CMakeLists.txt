#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_include_directories(OPS_Partition PUBLIC 
	${OPS_SRC_DIR}/logging
	${CMAKE_CURRENT_LIST_DIR}
)

target_sources(OPS_Partition
    PRIVATE
      Metis.cpp
    PUBLIC
      Metis.h
      GraphPartitioner.h
)
