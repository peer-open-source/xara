#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

# 1)
# Standalone executables that implements a Tcl shell, and
# loads the Tcl extension library compiled by target OpenSeesRT.
add_executable(xara_repl main.cpp)
set_target_properties(xara_repl PROPERTIES 
  OUTPUT_NAME "xara"
)

if (WIN32)
  target_link_libraries(xara_repl PRIVATE Netapi32.lib)
  target_link_libraries(xara_repl PRIVATE ${TCL_LIBRARIES} OpenSeesRT)
else()
  target_link_libraries(xara_repl PRIVATE ${TCL_LIBRARY} OpenSeesRT)
endif()

# 2)
# Standalone executables that implements a Tcl shell, and
# and directly sets up xara by calling Openseesrt_Init(),
# which is staticly linked from target XaraStaticRT
add_executable(XaraRT XaraMain.cpp)
set_target_properties(XaraRT PROPERTIES 
  OUTPUT_NAME "xaras"
)

target_link_libraries(XaraRT PUBLIC ${TCL_LIBRARY} XaraStaticRT)



add_executable(XaraMP XaraMain.cpp)
set_target_properties(XaraMP PROPERTIES 
  OUTPUT_NAME "xaramp"
)

target_link_libraries(XaraMP PUBLIC ${TCL_LIBRARY} XaraStaticMP)
