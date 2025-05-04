#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_sources(OPS_Material
    PRIVATE
        CycLiqCPPlaneStrain.cpp
        CycLiqCPSPPlaneStrain.cpp
        J2PlaneStrain.cpp
        J2PlaneStress.cpp
        PlaneStressUserMaterial.cpp
        PlaneStressSimplifiedJ2.cpp
        PlasticDamageConcretePlaneStress.cpp
        # PlasticDamageConcretePlaneStressThermal.cpp
    PUBLIC
        CycLiqCPPlaneStrain.h
        CycLiqCPSPPlaneStrain.h
        J2PlaneStrain.h
        J2PlaneStress.h
        PlasticDamageConcretePlaneStress.h
        # PlasticDamageConcretePlaneStressThermal.h
)

add_subdirectory(./reinforcedConcretePlaneStress/)

