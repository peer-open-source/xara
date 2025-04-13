#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_sources(OPS_Material
  PRIVATE
    NDMaterial.cpp
    PlasticDamageConcrete3d.cpp
    LinearCap.cpp
    FSAM.cpp
    AcousticMedium.cpp
    ASDConcrete3DMaterial.cpp
    ElasticPlaneStress.cpp
    VonPapaDamage.cpp
    ElasticIsotropicAxiSymm.cpp
    ElasticIsotropicMaterial.cpp
    ElasticOrthotropicMaterial.cpp
    ElasticOrthotropicPlaneStress.cpp
    PressureDependentElastic3D.cpp
    J2AxiSymm.cpp
    J2Plasticity.cpp

    J2PlateFiber.cpp
    J2PlateFibre.cpp
    PlateFiberMaterial.cpp
    ElasticIsotropicPlateFiber.cpp
    PlateFromPlaneStressMaterial.cpp
    PlateRebarMaterial.cpp
    PlateFiberMaterialThermal.cpp
    PlateFromPlaneStressMaterialThermal.cpp

    CycLiqCP.cpp
    CycLiqCP3D.cpp
    CycLiqCPSP.cpp
    CycLiqCPSP3D.cpp
    # WrapperNDMaterial.cpp
    MinMaxNDMaterial.cpp
    InitStressNDMaterial.cpp
    InitStrainNDMaterial.cpp
    ConcreteMcftNonLinear5.cpp
    ConcreteMcftNonLinear7.cpp
    CapPlasticity.cpp
    ConcreteS.cpp
    SimplifiedJ2.cpp
    LinearCap.cpp
    AcousticMedium.cpp
    DruckerPrager3DThermal.cpp
    DruckerPragerThermal.cpp
    ElasticIsotropic3DThermal.cpp
    ElasticIsotropicMaterialThermal.cpp
    J2PlasticityThermal.cpp
    UVCmultiaxial.cpp
    UVCplanestress.cpp
    PlateRebarMaterialThermal.cpp
    OrthotropicMaterial.cpp
    Series3DMaterial.cpp
    Parallel3DMaterial.cpp
  PUBLIC
    NDMaterial.h
    PlasticDamageConcrete3d.h
    LinearCap.h
    FSAM.h
    AcousticMedium.h
    ASDConcrete3DMaterial.h
    ElasticIsotropicAxiSymm.h

    ElasticIsotropicMaterial.h
    ElasticIsotropicPlateFiber.h
    ElasticOrthotropicMaterial.h
    ElasticOrthotropicPlaneStress.h
    ElasticPlaneStress.h
    VonPapaDamage.h
    PressureDependentElastic3D.h
    J2AxiSymm.h
    J2Plasticity.h
    J2PlateFiber.h
    J2PlateFibre.h
    PlateFiberMaterial.h
    PlateRebarMaterial.h
    CycLiqCP.h
    CycLiqCP3D.h
    CycLiqCPSP.h
    CycLiqCPSP3D.h
    # WrapperNDMaterial.h
    MinMaxNDMaterial.h
    InitStressNDMaterial.h
    InitStrainNDMaterial.h
    ConcreteMcftNonLinear5.h
    ConcreteMcftNonLinear7.h
    CapPlasticity.h
    ConcreteS.h
    SimplifiedJ2.h
    LinearCap.h
    AcousticMedium.h
    DruckerPrager3DThermal.h
    DruckerPragerThermal.h
    ElasticIsotropic3DThermal.h
    ElasticIsotropicMaterialThermal.h
    J2PlasticityThermal.h
    PlateFiberMaterialThermal.h
    UVCmultiaxial.h
    UVCplanestress.h
    PlateRebarMaterialThermal.h
    PlateFromPlaneStressMaterial.h
    PlateFromPlaneStressMaterialThermal.h
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

