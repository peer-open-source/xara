#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_sources(OPS_Runtime PRIVATE commands.cpp)

target_sources(OPS_Element
  PRIVATE
    Truss.cpp
    TrussSection.cpp
    CorotTruss.cpp
    CorotTrussSection.cpp
  PUBLIC
    Truss.h
    TrussSection.h
    CorotTruss.h
    CorotTrussSection.h
)

target_sources(OPS_Element
  PRIVATE
    InertiaTruss.cpp
  PUBLIC
    InertiaTruss.h
)

# These elements include specific material headers
#add_library(OPS_ConcreteTruss OBJECT)
target_include_directories(OPS_Element PUBLIC "${OPS_SRC_DIR}/material/uniaxial/concrete")
target_sources(OPS_Element
  PRIVATE
    Truss2.cpp
    N4BiaxialTruss.cpp
    CorotTruss2.cpp
  PUBLIC
    N4BiaxialTruss.h
    CorotTruss2.h
    Truss2.h
)

target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

