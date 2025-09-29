#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_link_libraries(OPS_SysOfEqn PRIVATE ARPACK)

target_sources(OPS_SysOfEqn
  PRIVATE
    ArpackSOE.cpp
    ArpackSolver.cpp
    EigenSOE.cpp
    EigenSolver.cpp
    FullGenEigenSOE.cpp
    FullGenEigenSolver.cpp
    SymBandEigenSOE.cpp
    SymBandEigenSolver.cpp
    SymmGeneralizedEigenSolver.cpp
    SymmGeneralizedEigenSOE.cpp
  PUBLIC
    ArpackSOE.h
    ArpackSolver.h
    EigenSOE.h
    EigenSolver.h
    FullGenEigenSOE.h
    FullGenEigenSolver.h
    SymBandEigenSOE.h
    SymBandEigenSolver.h
)
