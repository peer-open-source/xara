#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Material
  PRIVATE
    FeapMaterial.cpp
    FeapMaterial01.cpp
    FeapMaterial02.cpp
    FeapMaterial03.cpp
  PUBLIC
    FeapMaterial.h
    FeapMaterial01.h
    FeapMaterial02.h
    FeapMaterial03.h
)

add_library(OPS_Material_nD_Feap_f)
set_property(TARGET OPS_Material_nD_Feap_f PROPERTY POSITION_INDEPENDENT_CODE 1)
target_sources(OPS_Material_nD_Feap_f
  PRIVATE
    feapCommon.f
    matl01.f
    matl02.f
    matl03.f
)

target_link_libraries(OPS_Material PUBLIC OPS_Material_nD_Feap_f)
target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

