#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

add_executable(OpenSees tclMain.cpp)

if (WIN32)
  target_link_libraries(OpenSees PRIVATE Netapi32.lib)
  target_link_libraries(OpenSees PRIVATE ${TCL_LIBRARIES} OpenSeesRT)
else()
  target_link_libraries(OpenSees PRIVATE ${TCL_LIBRARY} OpenSeesRT)
endif()

