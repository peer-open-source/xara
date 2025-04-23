#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

add_library(OPS_Element_Masonry OBJECT)
target_sources(OPS_Element_Masonry
  PUBLIC
    BeamGT.cpp
    MasonPan12.cpp
    MasonPan3D.cpp
  PRIVATE
    BeamGT.h
    MasonPan12.h
    MasonPan3D.h
)    

if (UNIX)
  set_property(SOURCE
      BeamGT.cpp
      MasonPan12.cpp
      MasonPan3D.cpp
      APPEND PROPERTY COMPILE_FLAGS "-w"
  )
endif()

target_link_libraries(OPS_Element PRIVATE OPS_Element_Masonry)
