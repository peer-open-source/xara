#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Matrix
  PRIVATE
    cbdi.cpp
    shp3d.cpp
    shp3dv.cpp
    nurbs.cpp
  PUBLIC
    cbdi.h
    shp3d.h
    shp3dv.h
    nurbs.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_include_directories(OPS_Matrix PUBLIC ${CMAKE_CURRENT_LIST_DIR})


