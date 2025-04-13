#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_sources(OPS_Material
    PRIVATE
        ElasticIsotropicThreeDimensional.cpp
        ElasticOrthotropicThreeDimensional.cpp
        IncrementalElasticIsotropicThreeDimensional.cpp
        J2ThreeDimensional.cpp
        J2ThreeDimensionalThermal.cpp
    PUBLIC
        ElasticIsotropicThreeDimensional.h
        ElasticOrthotropicThreeDimensional.h
        IncrementalElasticIsotropicThreeDimensional.h
        J2ThreeDimensional.h
        J2ThreeDimensionalThermal.h
)

