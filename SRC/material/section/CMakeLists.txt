#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Material
  PRIVATE
    ASDCoupledHinge3D.cpp
    BiaxialHysteretic.cpp
    Bidirectional.cpp
#   ElasticBDShearSection2d.cpp
    ElasticSection2d.cpp
    ElasticSection3d.cpp
    ElasticShearSection2d.cpp 
#   ElasticShearSection3d.cpp
#   ElasticTubeSection3d.cpp
#   ElasticWarpingShearSection2d.cpp
    Elliptical2.cpp
    FiberSection2d.cpp
    FiberSection2dInt.cpp
    FiberSection2dThermal.cpp

    FiberSection3d.cpp
    FiberSectionWarping3d.cpp    
    FiberSectionAsym3d.cpp
    FiberSection3dThermal.cpp
#   FiberSectionGJ.cpp
#   FiberSectionGJThermal.cpp    
    GenericSection1d.cpp
    Isolator2spring.cpp
    NDFiberSection2d.cpp
    NDFiberSection3d.cpp
    NDFiberSectionWarping2d.cpp
    ParallelSection.cpp
    SectionAggregator.cpp
    SectionForceDeformation.cpp
#   TimoshenkoSection3d.cpp
  PUBLIC
    ASDCoupledHinge3D.h
    BiaxialHysteretic.h
    Bidirectional.h
#   ElasticBDShearSection2d.h
    ElasticSection2d.h 
    ElasticSection3d.h
    ElasticShearSection2d.h 
#   ElasticShearSection3d.h
#   ElasticTubeSection3d.h
#   ElasticWarpingShearSection2d.h
    Elliptical2.h
    FiberSection2d.h
    FiberSection2dInt.h
    FiberSection2dThermal.h
    FiberSection3d.h
    FiberSectionWarping3d.h    
    FiberSectionAsym3d.h
    FiberSection3dThermal.h
#   FiberSectionGJ.h
#   FiberSectionGJThermal.h
    GenericSection1d.h
    Isolator2spring.h
    NDFiberSection2d.h
    NDFiberSection3d.h
    NDFiberSectionWarping2d.h
    ParallelSection.h
    SectionAggregator.h
    SectionForceDeformation.h
#   TimoshenkoSection3d.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(yieldSurface)
add_subdirectory(Shell)
add_subdirectory(LayeredMembraneSection)

