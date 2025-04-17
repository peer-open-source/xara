#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Runtime PRIVATE commands.cpp)

target_sources(OPS_Analysis
    PRIVATE
       EigenIntegrator.cpp                  # Jun Peng
       IncrementalIntegrator.cpp
       Integrator.cpp
       StaticIntegrator.cpp
       TransientIntegrator.cpp
    PUBLIC
       EigenIntegrator.h
       IncrementalIntegrator.h
       Integrator.h
       StaticIntegrator.h
       TransientIntegrator.h
)


add_subdirectory(Dynamic)
add_subdirectory(Static)
add_subdirectory(Hybrid)
target_include_directories(OPS_Analysis PUBLIC ${CMAKE_CURRENT_LIST_DIR})

