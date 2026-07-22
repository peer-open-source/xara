#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Material
  PRIVATE
    stressDensity.cpp
    StressDensityModel2D.cpp 
    StressDensityModel3D.cpp
  PUBLIC
    stressDensity.h
)


target_sources(OPS_MaterialFortran
  PRIVATE
    SDM-UC.f
    SDM3D.f90
    SDM2D.f90
    MODULE_ALL_INTERFACES_3D.f90
    MODULE_ALL_INTERFACES_2D.f90
)

if (APPLE) 
  target_link_options(OPS_MaterialFortran
    PUBLIC
    -Wl,-no_compact_unwind
  )
endif()


target_include_directories(OPS_Material PRIVATE ${CMAKE_CURRENT_LIST_DIR})

