#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
# 
#==============================================================================

if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  target_compile_options(OPS_Domain PRIVATE -Wno-unused-parameter)
endif()

add_subdirectory(component)
add_subdirectory(domain)
add_subdirectory(pattern)
add_subdirectory(load)
add_subdirectory(node)
add_subdirectory(subdomain)
add_subdirectory(constraints)
add_subdirectory(groundMotion)
add_subdirectory(region)
add_subdirectory(partitioner)
add_subdirectory(peri)
