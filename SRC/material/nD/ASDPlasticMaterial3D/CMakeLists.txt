#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
find_package(Eigen3 NO_MODULE)
if (TARGET Eigen3::Eigen)

target_compile_definitions(OPS_Material PRIVATE OPS_USE_ASDPlasticMaterials)
target_link_libraries(OPS_Material PRIVATE Eigen3::Eigen)
target_include_directories(OPS_Material PRIVATE "${OPS_BUNDLED_DIR}/eigenAPI")

target_sources(OPS_Material
    PRIVATE
#     ASDPlasticMaterialGlobals.cpp
      OPS_AllASDPlasticMaterials.cpp
    PUBLIC
      AllASDPlasticMaterials.h     
      ASDPlasticMaterial.h        
      ElasticityBase.h    
      YieldFunctionBase.h
      ASDPlasticMaterialGlobals.h
      ASDPlasticMaterialTraits.h  
      PlasticFlowBase.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(ElasticityModels)
add_subdirectory(PlasticFlowDirections)
add_subdirectory(YieldFunctions)

endif (TARGET Eigen3::Eigen)

