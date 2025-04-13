#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_link_libraries(OPS_Material PRIVATE OPS_Domain OPS_Logging)
target_sources(OPS_Runtime PRIVATE commands.cpp)

target_sources(OPS_Runtime PRIVATE commands.cpp)

target_sources(OPS_Material
    PRIVATE
      Material.cpp
    PUBLIC
      Material.h
)

target_include_directories(OPS_Material
  PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
  PRIVATE
    ${CMAKE_CURRENT_LIST_DIR}/Solid
    ${CMAKE_CURRENT_LIST_DIR}/nD/cyclicSoil
    ${CMAKE_CURRENT_LIST_DIR}/nD/feap
    ${CMAKE_CURRENT_LIST_DIR}/nD/matCMM
    ${CMAKE_CURRENT_LIST_DIR}/nD/soil
    ${CMAKE_CURRENT_LIST_DIR}/section/yieldSurface
    ${CMAKE_CURRENT_LIST_DIR}/Yield/evolution
    ${CMAKE_CURRENT_LIST_DIR}/Yield/plasticHardeningMaterial
    ${CMAKE_CURRENT_LIST_DIR}/Yield/yieldSurfaceBC
    ${CMAKE_CURRENT_LIST_DIR}/uniaxial/backbone
    ${CMAKE_CURRENT_LIST_DIR}/uniaxial/drain
    ${CMAKE_CURRENT_LIST_DIR}/uniaxial/fedeas
    ${CMAKE_CURRENT_LIST_DIR}/uniaxial/limitState
    ${CMAKE_CURRENT_LIST_DIR}/uniaxial/limitState/limitCurve
    ${CMAKE_CURRENT_LIST_DIR}/uniaxial/PY
    ${CMAKE_CURRENT_LIST_DIR}/uniaxial/snap
    ${CMAKE_CURRENT_LIST_DIR}/uniaxial/stiffness
    ${CMAKE_CURRENT_LIST_DIR}/uniaxial/strength
    ${CMAKE_CURRENT_LIST_DIR}/uniaxial/unloading
    ${CMAKE_CURRENT_LIST_DIR}/nD/stressDensityModel
    ${CMAKE_CURRENT_LIST_DIR}/nD/UANDESmaterials
    ${CMAKE_CURRENT_LIST_DIR}/nD/UWmaterials
)
add_subdirectory(Frame)
add_subdirectory(Shell)
add_subdirectory(uniaxial)
add_subdirectory(section)
add_subdirectory(elastic)
add_subdirectory(plastic)
add_subdirectory(Solid)
add_subdirectory(Plane)
add_subdirectory(Yield)
add_subdirectory(nD)

