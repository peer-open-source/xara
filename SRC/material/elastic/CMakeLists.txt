
target_sources(OPS_Runtime PRIVATE commands.cpp isotropy.cpp)

target_sources(OPS_Material 
  PRIVATE
    ElasticIsotropicPlaneStrain2D.cpp
    ElasticIsotropicPlaneStress2D.cpp
    ElasticIsotropicBeamFiber2d.cpp
    ElasticIsotropicBeamFiber.cpp
  PUBLIC
    ElasticIsotropicPlaneStrain2D.h
    ElasticIsotropicPlaneStress2D.h
    ElasticIsotropicBeamFiber2d.h
    ElasticIsotropicBeamFiber.h
)

target_include_directories(OPS_Material
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
)
