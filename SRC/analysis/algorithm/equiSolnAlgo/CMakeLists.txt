#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
add_library(OPS_Algorithm OBJECT)
target_link_libraries(OPS_Algorithm PRIVATE OPS_Logging)
target_include_directories(OPS_Analysis  PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_include_directories(OPS_Algorithm PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_sources(OPS_Algorithm
    PRIVATE
      EquiSolnAlgo.cpp 
      ExpressNewton.cpp
      Linear.cpp 
      NewtonRaphson.cpp
      ModifiedNewton.cpp 
      NewtonLineSearch.cpp 
      Broyden.cpp 
      BFGS.cpp
      KrylovNewton.cpp 
      PeriodicNewton.cpp
      AcceleratedNewton.cpp
      NewtonHallM.cpp
    PUBLIC
      EquiSolnAlgo.h 
      ExpressNewton.h
      Linear.h 
      NewtonRaphson.h
      ModifiedNewton.h 
      NewtonLineSearch.h 
      Broyden.h 
      BFGS.h
      KrylovNewton.h 
      PeriodicNewton.h 
      AcceleratedNewton.h
      NewtonHallM.h
)

add_subdirectory(accelerator)
add_subdirectory(search)

