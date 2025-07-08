#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================


target_link_libraries(OPS_SysOfEqn PRIVATE UMFPACK AMD)

target_sources(OPS_Runtime PRIVATE dispatch.cpp)

target_sources(OPS_SysOfEqn
  PRIVATE
    UmfpackGenLinSOE.cpp
    UmfpackGenLinSolver.cpp
)

target_include_directories(OPS_SysOfEqn PUBLIC ${CMAKE_CURRENT_LIST_DIR})

