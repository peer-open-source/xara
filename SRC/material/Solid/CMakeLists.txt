#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_sources(OPS_Material
    PRIVATE
        ElasticOrthotropicThreeDimensional.cpp
        IncrementalElasticIsotropicThreeDimensional.cpp
        J2ThreeDimensionalThermal.cpp
    PUBLIC
        ElasticOrthotropicThreeDimensional.h
        IncrementalElasticIsotropicThreeDimensional.h
        J2ThreeDimensionalThermal.h
)

