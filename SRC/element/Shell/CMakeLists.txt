#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_link_libraries(OPS_Runtime PRIVATE OPS_Element)

target_sources(OPS_Element
    PRIVATE
      ASDShellQ4.cpp
      ASDShellT3.cpp

      ShellANDeS.cpp
      ShellDKGQ.cpp
      ShellDKGT.cpp
      ShellMITC4.cpp
      ShellMITC4Thermal.cpp
      ShellMITC9.cpp
#     IGAShellMITC9.cpp
      ShellNLDKGQ.cpp
      ShellNLDKGQThermal.cpp
      ShellNLDKGT.cpp
    PUBLIC
      ASDShellQ4.h
      ASDShellT3.h

      ShellANDeS.h
      ShellDKGQ.h
      ShellDKGT.h
      ShellMITC4.h
      ShellMITC4Thermal.h
      ShellMITC9.h
      ShellNLDKGQ.h
      ShellNLDKGQThermal.h
      ShellNLDKGT.h
)
target_include_directories(OPS_Element PUBLIC ${CMAKE_CURRENT_LIST_DIR})

