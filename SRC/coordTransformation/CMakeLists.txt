#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_link_libraries(OPS_Transform PRIVATE OPS_Logging OPS_Recorder)

target_sources(OPS_Transform
  PRIVATE
    CrdTransf.cpp
    LinearCrdTransf2d.cpp
    LinearCrdTransf3d.cpp
    PDeltaCrdTransf2d.cpp
    PDeltaCrdTransf3d.cpp
    CorotCrdTransf2d.cpp
    CorotCrdTransf3d.cpp
    CorotCrdTransfWarping2d.cpp
    CorotCrdTransfWarping3d.cpp
    LinearCrdTransf2dInt.cpp
#   LinearCrdTransf2d02.cpp
#   LinearCrdTransf3d02.cpp

  PUBLIC
    CrdTransf.h
    LinearCrdTransf2d.h
    LinearCrdTransf3d.h
    PDeltaCrdTransf2d.h
    PDeltaCrdTransf3d.h
    CorotCrdTransf2d.h
    CorotCrdTransf3d.h
    CorotCrdTransfWarping2d.h
    LinearCrdTransf2dInt.h
#   LinearCrdTransf2d02.h
#   LinearCrdTransf3d02.h

    ASDEICR.h
    Shell/ASDShellQ4CorotationalTransformation.h
    Shell/ASDShellQ4LocalCoordinateSystem.h
    Shell/ASDShellQ4Transformation.h
)

target_include_directories(OPS_Transform 
  PUBLIC 
    ${CMAKE_CURRENT_LIST_DIR}
    ${CMAKE_CURRENT_LIST_DIR}/Shell
)

add_subdirectory(Frame)

