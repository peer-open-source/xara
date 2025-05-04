#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Material 
  PRIVATE
    ElasticIsotropicPlaneStrain2D.cpp
    ElasticIsotropicThreeDimensional.cpp
    ElasticIsotropicPlaneStress2D.cpp
    ElasticIsotropicBeamFiber2d.cpp
    ElasticIsotropicBeamFiber.cpp
    ElasticIsotropicAxiSymm.cpp
    ElasticIsotropicMaterial.cpp
    ElasticIsotropicPlateFiber.cpp
    ElasticIsotropic3DThermal.cpp
    ElasticIsotropicMaterialThermal.cpp
  PUBLIC
#   ElasticIsotropicThreeDimensional.h
    ElasticIsotropicPlaneStrain2D.h
    ElasticIsotropicPlaneStress2D.h
    ElasticIsotropicBeamFiber2d.h
    ElasticIsotropicBeamFiber.h
)

target_include_directories(OPS_Material
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
)
