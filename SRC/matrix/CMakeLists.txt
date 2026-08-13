#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_link_libraries(OPS_Matrix PRIVATE Xara::FastMath)

target_sources(OPS_Matrix
  PRIVATE
    ID.cpp
    Matrix.cpp
    Vector.cpp
    R3vectors.cpp
    TriMatrix.cpp
  PUBLIC
    ID.h
    Matrix.h
    Vector.h
    R3vectors.h
    TriMatrix.h
)

add_subdirectory(routines)

target_include_directories(OPS_Matrix 
  PUBLIC ${CMAKE_CURRENT_LIST_DIR}
  PRIVATE ${OPS_Includes}
)
