#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_sources(OPS_Material
  PRIVATE
    CreepMaterial.cpp
    ContinuumUniaxial.cpp
    InitStrainMaterial.cpp
    InitStressMaterial.cpp
    InitStressNDMaterial.cpp
    InitStrainNDMaterial.cpp
    ParallelMaterial.cpp
    SeriesMaterial.cpp
    OrthotropicMaterial.cpp
    # Plate
    PlateFiberMaterial.cpp
    PlateFromPlaneStressMaterial.cpp
    PlateRebarMaterial.cpp
    PlateFiberMaterialThermal.cpp
    PlateFromPlaneStressMaterialThermal.cpp
    PlateRebarMaterialThermal.cpp
    # Plane
    PlaneStrainMaterial.cpp
    PlaneStressLayeredMaterial.cpp
    PlaneStressMaterial.cpp
    PlaneStressRebarMaterial.cpp
    # Frame
    BeamFiberMaterial.cpp
    BeamFiberMaterial2d.cpp
    BeamFiberMaterial2dPS.cpp
    # FrameWrapper.cpp
  PUBLIC
    ParallelMaterial.h
    ContinuumUniaxial.h
    CreepMaterial.h
)
