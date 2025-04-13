#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
add_executable(trussExample 
  # ${}/truss.cpp
  main.cpp
)

target_link_libraries(trussExample #OPS_Element_truss 
  G3_API G3)

add_test(TrussTest trussExample COMMAND trussExample)

