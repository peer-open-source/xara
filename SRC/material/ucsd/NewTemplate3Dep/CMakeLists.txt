add_library(UCSD_Template3D)

target_link_libraries(UCSD_Template3D PRIVATE UCSD_Array)
target_include_directories(UCSD_Template3D PRIVATE ${OPS_Includes})

target_sources(UCSD_Template3D
  PRIVATE
    AF_Eij.cpp
    CC_Ev.cpp
    CC_PF.cpp
    CC_YF.cpp
    DM04_alpha_Eij.cpp
    DM04_Elastic.cpp
    DM04_PF.cpp
    DM04_YF.cpp
    DM04_z_Eij.cpp
    DP_PF.cpp
    DP_YF.cpp
    ElasticState.cpp
    elnp_Elastic.cpp
    Isotropic_Elastic.cpp
    Linear_Eeq.cpp
    Linear_Eij.cpp
    MaterialParameter.cpp
    NewTemplate3Dep.cpp
    PlasticFlow.cpp
    PressureDependent_Elastic.cpp
    RMC_PF.cpp
    RMC_YF.cpp
    ScalarEvolution.cpp
    TensorEvolution.cpp
    VM_PF.cpp
    VM_YF.cpp
    YieldFunction.cpp
  PUBLIC
)

target_link_libraries(UCSD_Template3D 
  PUBLIC UCSD_Array OPS_Material)

add_executable(FT3Dep_test EXCLUDE_FROM_ALL FT3Dep_test.cpp)
target_link_libraries(FT3Dep_test 
  PRIVATE UCSD_Template3D OpenSeesRT OPS_Numerics)
