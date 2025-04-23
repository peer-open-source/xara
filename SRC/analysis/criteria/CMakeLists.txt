#==============================================================================
#
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_link_libraries(OPS_ConvergenceTest PRIVATE OPS_Algorithm OPS_Logging)

target_sources(OPS_ConvergenceTest
  PRIVATE
    ConvergenceTest.cpp
    CTestNormUnbalance.cpp
    CTestNormDispIncr.cpp
    CTestEnergyIncr.cpp
    CTestRelativeNormUnbalance.cpp
    CTestRelativeNormDispIncr.cpp
    CTestRelativeEnergyIncr.cpp
    CTestFixedNumIter.cpp
    CTestRelativeTotalNormDispIncr.cpp
    NormDispAndUnbalance.cpp
    NormDispOrUnbalance.cpp
  PUBLIC
    ConvergenceTest.h
    CTestNormUnbalance.h
    CTestNormDispIncr.h
    CTestEnergyIncr.h
    CTestRelativeNormUnbalance.h
    CTestRelativeNormDispIncr.h
    CTestRelativeEnergyIncr.h
    CTestFixedNumIter.h
    CTestRelativeTotalNormDispIncr.h
    NormDispAndUnbalance.h
    NormDispOrUnbalance.h
)


target_include_directories(OPS_ConvergenceTest PUBLIC ${CMAKE_CURRENT_LIST_DIR})

