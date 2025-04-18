#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_link_libraries(OPS_Element 
  PUBLIC  OPS_Domain 

  PRIVATE OPS_Transform
          OPS_FrameQuadrature
          OPS_Logging
)

target_include_directories(OPS_Element PUBLIC
  ${OPS_SRC_DIR}/element/Boundary
  ${OPS_SRC_DIR}/element/Plane
  ${OPS_SRC_DIR}/element/tetrahedron
  ${OPS_SRC_DIR}/element/adapter
  ${OPS_SRC_DIR}/element/Brick
  ${OPS_SRC_DIR}/element/Joint
  ${OPS_SRC_DIR}/element/Cable
  ${OPS_SRC_DIR}/element/Shell
  ${OPS_SRC_DIR}/element/Truss

  ${OPS_SRC_DIR}/element/Frame/UpdatedLagrange
  ${OPS_SRC_DIR}/element/Frame/GradientInelastic
  ${OPS_SRC_DIR}/element/Frame/ForceInterp
  ${OPS_SRC_DIR}/element/Frame/MixedInterp
  ${OPS_SRC_DIR}/element/Frame/DisplInterp
  ${OPS_SRC_DIR}/element/Frame/LHMYS
  ${OPS_SRC_DIR}/element/componentElement

  ${OPS_SRC_DIR}/element/ASDEA/CEqElement

  ${OPS_SRC_DIR}/element/community
  ${OPS_SRC_DIR}/element/community/HUelements

  ${OPS_SRC_DIR}/element/Bearing/elastomeric
  ${OPS_SRC_DIR}/element/Bearing/elastomeric
  ${OPS_SRC_DIR}/element/Bearing/friction
  ${OPS_SRC_DIR}/element/Bearing/friction/frictionModel

  ${OPS_SRC_DIR}/element/feap
  ${OPS_SRC_DIR}/element/masonry
  ${OPS_SRC_DIR}/element/mvlem
  ${OPS_SRC_DIR}/element/Other/PML
  ${OPS_SRC_DIR}/element/Other/generic
  ${OPS_SRC_DIR}/element/Other/pyMacro
  ${OPS_SRC_DIR}/element/surfaceLoad
  ${OPS_SRC_DIR}/element/Link
)

target_sources(OPS_Element
    PRIVATE
      Element.cpp
      ElementalLoad.cpp
      Other/WrapperElement.cpp
    PUBLIC
      Element.h
      ElementalLoad.h
      Other/WrapperElement.h
)

target_sources(OPS_Utilities
    PRIVATE
      Information.cpp
    PUBLIC
      Information.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

add_subdirectory(Link)
add_subdirectory(tetrahedron)
add_subdirectory(Plane)
add_subdirectory(Shell)
add_subdirectory(Brick)
add_subdirectory(Truss)
add_subdirectory(Cable)
add_subdirectory(Joint)
add_subdirectory(Point)
add_subdirectory(Frame)
add_subdirectory(Frame/LHMYS)
add_subdirectory(componentElement)

add_subdirectory(adapter)
add_subdirectory(Boundary)

add_subdirectory(Bearing/friction)
add_subdirectory(Bearing/elastomeric)

add_subdirectory(ASDEA/CEqElement)
add_subdirectory(ASDEA/absorbentBoundaries)

add_subdirectory(Other/PML)
add_subdirectory(Other/generic)
add_subdirectory(Other/pyMacro)

add_subdirectory(surfaceLoad)

add_subdirectory(feap)

add_subdirectory(mvlem)
add_subdirectory(masonry)
# add_subdirectory(IGA)

add_subdirectory(community/HUelements)
add_subdirectory(community/XMUelements)
add_subdirectory(community/UWelements)

add_subdirectory(mefi)
