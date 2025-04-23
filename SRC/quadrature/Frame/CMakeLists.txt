#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
add_library(OPS_FrameQuadrature)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_include_directories(OPS_FrameQuadrature PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_sources(OPS_Runtime PRIVATE commands.cpp)

target_sources(OPS_FrameQuadrature
  PRIVATE
      BeamIntegration.cpp
      ConcentratedCurvatureBeamIntegration.cpp
      ConcentratedPlasticityBeamIntegration.cpp
      LobattoBeamIntegration.cpp
      LegendreBeamIntegration.cpp
      RadauBeamIntegration.cpp
      NewtonCotesBeamIntegration.cpp
      HingeMidpointBeamIntegration.cpp
      HingeEndpointBeamIntegration.cpp
      HingeRadauBeamIntegration.cpp
      HingeRadauTwoBeamIntegration.cpp
      UserDefinedBeamIntegration.cpp
      SimpsonBeamIntegration.cpp
      UserDefinedHingeIntegration.cpp
      DistHingeIntegration.cpp
      RegularizedHingeIntegration.cpp
      TrapezoidalBeamIntegration.cpp
      CompositeSimpsonBeamIntegration.cpp
      FixedLocationBeamIntegration.cpp
      LowOrderBeamIntegration.cpp
      MidDistanceBeamIntegration.cpp
      #GaussQBeamIntegration.cpp
      ChebyshevBeamIntegration.cpp
  PUBLIC
      BeamIntegration.h
      ConcentratedCurvatureBeamIntegration.h
      ConcentratedPlasticityBeamIntegration.h
      LobattoBeamIntegration.h
      LegendreBeamIntegration.h
      RadauBeamIntegration.h
      NewtonCotesBeamIntegration.h
      HingeMidpointBeamIntegration.h
      HingeEndpointBeamIntegration.h
      HingeRadauBeamIntegration.h
      HingeRadauTwoBeamIntegration.h
      UserDefinedBeamIntegration.h
      SimpsonBeamIntegration.h
      UserDefinedHingeIntegration.h
      DistHingeIntegration.h
      RegularizedHingeIntegration.h
      TrapezoidalBeamIntegration.h
      CompositeSimpsonBeamIntegration.h
      FixedLocationBeamIntegration.h
      LowOrderBeamIntegration.h
      MidDistanceBeamIntegration.h
      #GaussQBeamIntegration.h
      ChebyshevBeamIntegration.h
)
