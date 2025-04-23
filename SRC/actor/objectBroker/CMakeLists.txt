#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# target_link_libraries(OPS_ObjectBroker PUBLIC G3_ObjectBroker OPS_Material OPS_Element)
# target_include_directories(OPS_ObjectBroker PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_sources(G3_ObjectBroker
  PRIVATE
    FEM_ObjectBroker.cpp
  PUBLIC
    FEM_ObjectBroker.h
)

add_library(OPS_ParallelBroker OBJECT)
target_sources(OPS_ParallelBroker
  PRIVATE
    FEM_ObjectBrokerAllClasses.cpp
  PUBLIC
    FEM_ObjectBrokerAllClasses.h
)
target_include_directories(OPS_ParallelBroker PRIVATE
 ${OPS_SRC_DIR}/coordTransformation/
 ${OPS_SRC_DIR}/element/Other/
 ${OPS_SRC_DIR}/damping/

 ${OPS_SRC_DIR}/material/uniaxial/backbone/ 
 ${OPS_SRC_DIR}/material/uniaxial/bearing/ 
 ${OPS_SRC_DIR}/material/uniaxial/concrete/ 
 ${OPS_SRC_DIR}/material/uniaxial/drain/ 
 ${OPS_SRC_DIR}/material/uniaxial/fedeas/ 
 ${OPS_SRC_DIR}/material/uniaxial/limitState/ 
 ${OPS_SRC_DIR}/material/uniaxial/PY/ 
 ${OPS_SRC_DIR}/material/uniaxial/snap/ 
 ${OPS_SRC_DIR}/material/uniaxial/steel/ 
 ${OPS_SRC_DIR}/material/uniaxial/stiffness/ 
 ${OPS_SRC_DIR}/material/uniaxial/strength/ 
 ${OPS_SRC_DIR}/material/uniaxial/unloading/ 
 ${OPS_SRC_DIR}/material/uniaxial/viscous/ 
 ${OPS_SRC_DIR}/material/uniaxial/wrapper/
)

target_link_libraries(OPS_ParallelBroker PUBLIC
  OPS_Material
  OPS_Element
  OPS_Analysis
)
