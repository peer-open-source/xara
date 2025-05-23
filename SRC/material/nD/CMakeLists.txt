#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Material
  PRIVATE
    NDMaterial.cpp
    LinearCap.cpp
    FSAM.cpp
    AcousticMedium.cpp
#   ElasticPlaneStress.cpp
    VonPapaDamage.cpp
    ElasticOrthotropicMaterial.cpp
    ElasticOrthotropicPlaneStress.cpp

    J2AxiSymm.cpp

    CycLiqCP.cpp
    CycLiqCP3D.cpp
    CycLiqCPSP.cpp
    CycLiqCPSP3D.cpp
    # WrapperNDMaterial.cpp
    MinMaxNDMaterial.cpp
    ConcreteMcftNonLinear5.cpp
    ConcreteMcftNonLinear7.cpp
    CapPlasticity.cpp
    ConcreteS.cpp
    LinearCap.cpp
    AcousticMedium.cpp
    DruckerPrager3DThermal.cpp
    DruckerPragerThermal.cpp
    J2PlasticityThermal.cpp
    OrthotropicMaterial.cpp
    Series3DMaterial.cpp
    Parallel3DMaterial.cpp
    PressureDependentElastic3D.cpp

  PUBLIC
    NDMaterial.h
    LinearCap.h
    FSAM.h
    AcousticMedium.h
    ElasticOrthotropicMaterial.h
    ElasticOrthotropicPlaneStress.h
#   ElasticPlaneStress.h
    VonPapaDamage.h
    PressureDependentElastic3D.h

    J2AxiSymm.h

    CycLiqCP.h
    CycLiqCP3D.h
    CycLiqCPSP.h
    CycLiqCPSP3D.h
    # WrapperNDMaterial.h
    MinMaxNDMaterial.h
    ConcreteMcftNonLinear5.h
    ConcreteMcftNonLinear7.h
    CapPlasticity.h
    ConcreteS.h
    LinearCap.h
    AcousticMedium.h
    DruckerPrager3DThermal.h
    DruckerPragerThermal.h
    J2PlasticityThermal.h
    OrthotropicMaterial.h
    Series3DMaterial.h
    Parallel3DMaterial.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(soil)
add_subdirectory(feap)
add_subdirectory(cyclicSoil)
add_subdirectory(UWmaterials)
add_subdirectory(ASDPlasticMaterial3D)
# add_subdirectory(matCMM)
add_subdirectory(stressDensityModel)
add_subdirectory(UANDESmaterials)
add_subdirectory(OrthotropicRotatingAngleConcreteT2DMaterial01)
add_subdirectory(SmearedSteelDoubleLayerT2DMaterial01)

