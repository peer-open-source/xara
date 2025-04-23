#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Element
    PRIVATE
        ForceBeamColumn2d.cpp
        ForceBeamColumn3d.cpp
        ForceBeamColumn2dThermal.cpp
        ElasticForceBeamColumn2d.cpp
        ElasticForceBeamColumn3d.cpp
        ElasticForceBeamColumnWarping2d.cpp
        ForceBeamColumnCBDI2d.cpp
        ForceBeamColumnCBDI3d.cpp	
        ForceBeamColumnWarping2d.cpp
    PUBLIC
        ForceBeamColumn2d.h
        ForceBeamColumn3d.h
        ForceBeamColumn2dThermal.h
        ElasticForceBeamColumn2d.h
        ElasticForceBeamColumn3d.h
        ElasticForceBeamColumnWarping2d.h
        ForceBeamColumnCBDI2d.h
        ForceBeamColumnCBDI3d.h	
        ForceBeamColumnWarping2d.h
)

#target_link_libraries(OPS_Element_forceBeamColumn PUBLIC OPS_Element_forceBeamColumn_f)
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

