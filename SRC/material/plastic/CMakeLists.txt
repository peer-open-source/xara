#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Material 
  PRIVATE
    J2Plasticity.cpp
    J2ThreeDimensional.cpp
    J2AxiSymm.cpp
    J2BeamFiber2d.cpp
    J2BeamFiber3d.cpp
    J2BeamThread3d.cpp
    J2PlateFiber.cpp
    J2PlateFibre.cpp
    UVCmultiaxial.cpp
    UVCplanestress.cpp
    SimplifiedJ2.cpp
    NonlinearJ2.cpp
  PUBLIC
    J2Plasticity.h
    SimplifiedJ2.h
    J2BeamFiber2d.h
    J2BeamFiber3d.h
)

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

if (EXISTS "${CMAKE_CURRENT_LIST_DIR}/GeneralizedJ2.cpp")
  target_sources(OPS_Material
    PRIVATE
      GeneralizedJ2.cpp
    PUBLIC
      GeneralizedJ2.h
  )
  target_compile_definitions(OPS_Runtime PRIVATE XARA_HAVE_GENERALIZEDJ2)
endif()
