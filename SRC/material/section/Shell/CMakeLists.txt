
target_sources(OPS_Material
  PRIVATE
    LayeredShellFiberSection.cpp
    DoubleMembranePlateFiberSection.cpp
    ElasticMembranePlateSection.cpp
#   ElasticPlateSection.cpp
    MembranePlateFiberSection.cpp
    # Thermal
    thermal/LayeredShellFiberSectionThermal.cpp
    thermal/MembranePlateFiberSectionThermal.cpp
  PUBLIC
    LayeredShellFiberSection.h
    DoubleMembranePlateFiberSection.h
    ElasticMembranePlateSection.h
#   ElasticPlateSection.h
    MembranePlateFiberSection.h
    thermal/LayeredShellFiberSectionThermal.h
    thermal/MembranePlateFiberSectionThermal.h
)
target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})
