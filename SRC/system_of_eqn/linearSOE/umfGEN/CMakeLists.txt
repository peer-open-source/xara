#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

# opensees_library(OPS_SysOfEqn_UMF
#   LINK ${UMFPACK_LIBRARIES} ${AMD_LIBRARIES})
target_link_libraries(OPS_SysOfEqn PRIVATE UMFPACK AMD)

target_sources(OPS_Runtime PRIVATE dispatch.cpp)

target_sources(OPS_SysOfEqn
  PRIVATE
    UmfpackGenLinSOE.cpp
    UmfpackGenLinSolver.cpp

# PUBLIC
#   UmfpackGenLinSOE.h
#   UmfpackGenLinSolver.h
)

target_include_directories(OPS_SysOfEqn PUBLIC ${CMAKE_CURRENT_LIST_DIR})
# target_include_directories(OPS_SysOfEqn_UMF PUBLIC ${CMAKE_CURRENT_LIST_DIR})

