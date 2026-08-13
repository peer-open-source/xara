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
target_link_libraries(OPS_SysOfEqn     PUBLIC SuperLU)

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
target_include_directories(OPS_SysOfEqn     PUBLIC ${CMAKE_CURRENT_LIST_DIR})

#
# SuperLU for distrubuted memory parallelism
#

add_library(OPS_SuperLU_Dist OBJECT EXCLUDE_FROM_ALL)

target_link_libraries(OPS_SuperLU_Dist PUBLIC SUPERLU_DIST OPS_SysOfEqn)
target_sources(OPS_SuperLU_Dist
    PRIVATE 
        DistributedSuperLU.cpp
        DistributedSparseGenColLinSOE.cpp
        DistributedSparseGenRowLinSOE.cpp
    PUBLIC
        DistributedSuperLU.h
        DistributedSparseGenColLinSOE.h
        DistributedSparseGenRowLinSOE.h
)

target_include_directories(OPS_SuperLU_Dist PUBLIC 
  ${CMAKE_CURRENT_LIST_DIR}
  ${OPS_Includes}
  ${OPS_BUNDLED_DIR}/SuperLU_DIST_4.3/SRC
)

