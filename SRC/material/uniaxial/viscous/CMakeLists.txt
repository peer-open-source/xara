

target_sources(OPS_Material
  PRIVATE
      Maxwell.cpp                      # Uses ops_Dt
      ViscousDamper.cpp
      ViscousMaterial.cpp
)


target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})
