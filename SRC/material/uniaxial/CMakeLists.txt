#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_sources(OPS_Material
  PRIVATE
    APDFMD.cpp
    APDMD.cpp
    APDVFD.cpp
    Bilin.cpp
    Bilin02.cpp
    DowelType.cpp
    Elastic2Material.cpp
    ElasticBilin.cpp
    ElasticMaterial.cpp
    ElasticMultiLinear.cpp
    ElasticPPMaterial.cpp
    ElasticPowerFunc.cpp
    IMKBilin.cpp
    IMKPeakOriented.cpp
    IMKPinching.cpp
    JankowskiImpact.cpp
    HystereticSMMaterial.cpp
    HystereticSmooth.cpp
    HystereticAsym.cpp
    MultiplierMaterial.cpp
    ResilienceLow.cpp
    Ratchet.cpp
    UVCuniaxial.cpp
    UniaxialMaterial.cpp
  PUBLIC
    APDMD.h
    APDVFD.h
    Bilin.h
    Bilin02.h
    DowelType.h
    Elastic2Material.h
    ElasticBilin.h
    ElasticMaterial.h
    ElasticMultiLinear.h
    ElasticPPMaterial.h
    ElasticPowerFunc.h
    IMKBilin.h
    IMKPeakOriented.h
    IMKPinching.h
    JankowskiImpact.h
    HystereticSMMaterial.h
    HystereticSmooth.h
    HystereticAsym.h
    MultiplierMaterial.h
    ResilienceLow.h
    Ratchet.h
    UVCuniaxial.h
    UniaxialMaterial.h
)


target_sources(OPS_Material
    PRIVATE
      ASD_SMA_3K.cpp
      ASDConcrete1DMaterial.cpp
      BackboneMaterial.cpp
      BarSlipMaterial.cpp
      BilinearOilDamper.cpp
      Bond_SP01.cpp
      CFSSSWP.cpp
      CFSWSWP.cpp
      CableMaterial.cpp
      Cast.cpp
      CubicSpline.cpp
      DamperMaterial.cpp
      DrainMaterial.cpp
      ECC01.cpp
      ENTMaterial.cpp
      EPPGapMaterial.cpp
      ElasticBDMaterial.cpp
      ElasticMaterialThermal.cpp
      FatigueMaterial.cpp
      #FedeasMaterial.cpp
      GNGMaterial.cpp
      HardeningMaterial.cpp
      #HardeningMaterial2.cpp
      HookGap.cpp
      HystereticMaterial.cpp
      HystereticPoly.cpp
      ImpactMaterial.cpp
      Maxwell.cpp                      # Uses ops_Dt
      MinMaxMaterial.cpp
      ModIMKPeakOriented.cpp
      ModIMKPeakOriented02.cpp
      ModIMKPinching.cpp
      ModIMKPinching02.cpp
      MultiLinear.cpp
      Neoprene.cpp
      OOHystereticMaterial.cpp
      PathIndependentMaterial.cpp
      Pinching4Material.cpp
      ResilienceMaterialHR.cpp
      SAWSMaterial.cpp
      SPSW02.cpp
      SelfCenteringMaterial.cpp
      ShearPanelMaterial.cpp
      SimpleFractureMaterial.cpp
      SLModel.cpp
      TensionOnlyMaterial.cpp
      UniaxialJ2Plasticity.cpp
      ViscousDamper.cpp
      ViscousMaterial.cpp
      ViscoelasticGap.cpp
      # WrapperUniaxialMaterial.cpp
      pyUCLA.cpp
      SMAMaterial.cpp
      DuctileFracture.cpp
)


target_sources(OPS_Material
  PRIVATE
   Masonry.cpp
   Masonryt.cpp
   Trilinwp.cpp
   Trilinwp2.cpp
   Trilinwpd.cpp
  PUBLIC
   Masonry.h
   Masonryt.h
   Trilinwp.h
   Trilinwp2.h
   Trilinwpd.h
)


add_library(OPS_Material_f)
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  target_compile_options(OPS_Material_f PRIVATE "-w" "-std=legacy")
endif()

# target_sources(OPS_Material_f PRIVATE
#      #DoddRestrepo.f
#      #STEELDR.f
#      #c14-SK-M.f
# )

target_link_libraries(OPS_Material PUBLIC OPS_Material_f)
add_subdirectory(abutment)
add_subdirectory(BoucWen)
add_subdirectory(BoucWenMG)
add_subdirectory(AxialSpring)
add_subdirectory(AlgebraicHysteresis)
add_subdirectory(backbone)
add_subdirectory(bearing)
add_subdirectory(limitState)
add_subdirectory(concrete)
add_subdirectory(steel)

add_subdirectory(PY)
add_subdirectory(snap)
add_subdirectory(fedeas)
add_subdirectory(drain)
add_subdirectory(legacy)

