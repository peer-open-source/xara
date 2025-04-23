#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

# TODO: SuperLU must be linked publically because
# SuperLU.h must include it. This could be made
# cleaner if SuperLU.h stored pointers to SuperLU 
# types, and not the types directly.
target_link_libraries(OPS_SysOfEqn PUBLIC SuperLU)
target_sources(OPS_SysOfEqn
  PRIVATE 
    SparseGenColLinSOE.cpp
    SparseGenColLinSolver.cpp
    SparseGenRowLinSOE.cpp
    SparseGenRowLinSolver.cpp
    SuperLU.cpp
  PUBLIC
    SparseGenColLinSOE.h
    SparseGenColLinSolver.h
    SparseGenRowLinSOE.h
    SparseGenRowLinSolver.h
    SuperLU.h
)

# target_sources(OPS_Parallel
#     PRIVATE 
#         DistributedSuperLU.cpp
#         DistributedSparseGenColLinSOE.cpp
#         DistributedSparseGenRowLinSOE.cpp
#     PUBLIC
#         DistributedSuperLU.h
#         DistributedSparseGenColLinSOE.h
#         DistributedSparseGenRowLinSOE.h
# )
target_include_directories(OPS_SysOfEqn PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_include_directories(OPS_Parallel PUBLIC ${CMAKE_CURRENT_LIST_DIR})

