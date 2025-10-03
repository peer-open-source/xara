

target_include_directories(UCSD_FiniteDeformation PRIVATE
  .
  ./fdEvolution
  ./fdFlow
  ./fdYield
)

target_sources(UCSD_FiniteDeformation
    PRIVATE
        BJMaterial.cpp
        FDdecoupledElastic3D.cpp
        FDEPState.cpp
        FiniteDeformationElastic3D.cpp
        FiniteDeformationEP3D.cpp
        LogWEnergy.cpp
        MooneyRivlinSimoWEnergy.cpp
        MooneyRivlinWEnergy.cpp
        NeoHookeanCompressible3D.cpp
        NeoHookeanWEnergy.cpp
        OgdenSimoWEnergy.cpp
        OgdenWEnergy.cpp
        SimoPisterWEnergy.cpp
        W.cpp
    PUBLIC
        BJMaterial.h
        FDdecoupledElastic3D.h
        FDEPState.h
        FiniteDeformationElastic3D.h
        FiniteDeformationEP3D.h
        LogWEnergy.h
        MooneyRivlinSimoWEnergy.h
        MooneyRivlinWEnergy.h
        NeoHookeanCompressible3D.h
        NeoHookeanWEnergy.h
        OgdenSimoWEnergy.h
        OgdenWEnergy.h
        SimoPisterWEnergy.h
        W.h
)

add_subdirectory(./fdEvolution/)
add_subdirectory(./fdFlow/)
add_subdirectory(./fdYield/)

target_link_libraries(UCSD_FiniteDeformation PRIVATE UCSD_Array)

add_executable(fd_test fd_test.cpp)
target_include_directories(fd_test PRIVATE
  .
  ./fdEvolution
  ./fdFlow
  ./fdYield
)

target_link_libraries(fd_test PRIVATE UCSD_FiniteDeformation UCSD_Array OpenSeesRT)
