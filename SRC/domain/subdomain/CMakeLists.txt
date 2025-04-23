#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Domain
  PRIVATE
    Subdomain.cpp
    SubdomainNodIter.cpp 
    ActorSubdomain.cpp
  PUBLIC
    Subdomain.h
    SubdomainNodIter.h 
    ActorSubdomain.h
)

target_sources(OPS_Domain
  PRIVATE
    ShadowSubdomain.cpp 
  PUBLIC
    ShadowSubdomain.h 
)

target_link_libraries(OPS_Domain PRIVATE OPS_Algorithm)

target_include_directories(OPS_Domain PUBLIC "${OPS_SRC_DIR}/analysis/analysis" ${CMAKE_CURRENT_LIST_DIR})

