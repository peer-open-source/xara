#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Runtime PRIVATE
    invoke.cpp
    invoke_uniaxial.cpp
    invoke_section.cpp
    invoke_stress.cpp
    invoke_material.cpp
)