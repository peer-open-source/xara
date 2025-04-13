
target_sources(OPS_Material
  PRIVATE
    LayeredShellFiberSection.cpp
    LayeredShellFiberSectionThermal.cpp
    DoubleMembranePlateFiberSection.cpp
    ElasticMembranePlateSection.cpp
#   ElasticPlateSection.cpp
    MembranePlateFiberSection.cpp
    MembranePlateFiberSectionThermal.cpp
  PUBLIC
    LayeredShellFiberSection.h
    LayeredShellFiberSectionThermal.h
    DoubleMembranePlateFiberSection.h
    ElasticMembranePlateSection.h
#   ElasticPlateSection.h
    MembranePlateFiberSection.h
    MembranePlateFiberSectionThermal.h
)
target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})
