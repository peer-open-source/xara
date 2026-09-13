
target_sources(OPS_Element
PRIVATE 
  ${CMAKE_CURRENT_LIST_DIR}/element/cpp/Macroelement3d/Macroelement3d/Macroelement3d.cpp
  ${CMAKE_CURRENT_LIST_DIR}/element/cpp/Macroelement3d/Macroelement3d/NoTensionSection3d.cpp
  ${CMAKE_CURRENT_LIST_DIR}/element/cpp/Macroelement3d/Macroelement3d/CohesiveSurface.cpp
  ${CMAKE_CURRENT_LIST_DIR}/element/cpp/Macroelement3d/Macroelement3d/GambarottaLagomarsinoModel.cpp
  ${CMAKE_CURRENT_LIST_DIR}/element/cpp/Macroelement3d/Macroelement3d/GenericDamagePlasticityShear.cpp
)

target_sources(OPS_Material
PRIVATE 
  ${CMAKE_CURRENT_LIST_DIR}/element/cpp/Macroelement3d/Macroelement3d/GenericDamagePlasticityShear.h
  ${CMAKE_CURRENT_LIST_DIR}/material/cpp/FrictionCohesion3d/FrictionCohesion3d/FrictionCohesion3d.cpp
  ${CMAKE_CURRENT_LIST_DIR}/material/cpp/CompressionDamage1d/CompressionDamage1d/CompressionDamage1d.cpp
)