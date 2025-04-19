#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================


target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_include_directories(OPS_Runtime  PUBLIC
  ${OPS_SRC_DIR}/material/uniaxial/stiffness
  ${OPS_SRC_DIR}/material/uniaxial/strength
  ${OPS_SRC_DIR}/material/uniaxial/unloading
  ${OPS_SRC_DIR}/material/uniaxial/backbone
)

target_sources(OPS_Material 
  PRIVATE
    MaterialState.cpp
  PUBLIC
    MaterialState.h
)

target_sources(OPS_Runtime PRIVATE
    "backbone.cpp"
    "StiffnessDegradationCommand.cpp"
    "StrengthDegradationCommand.cpp"
    "UnloadingRuleCommand.cpp"
)

add_subdirectory(stiffness)
add_subdirectory(strength)
add_subdirectory(unloading)