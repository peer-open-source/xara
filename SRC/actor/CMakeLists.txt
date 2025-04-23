#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# The objectBroker recipe uses its own target, so OPS_Actor
# won't see its directory
target_include_directories(OPS_Actor PUBLIC objectBroker)
target_link_libraries(OPS_Actor PRIVATE OPS_Logging)

add_subdirectory(address)
add_subdirectory(channel)
add_subdirectory(machineBroker)
add_subdirectory(message)
add_subdirectory(objectBroker)
add_subdirectory(shadow)
add_subdirectory(actor)




