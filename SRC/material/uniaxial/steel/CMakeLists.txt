#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_sources(OPS_Material_f PRIVATE STEELDR.f DoddRestrepo.f)

add_subdirectory(ReinforcingSteel)

target_sources(OPS_Material
    PRIVATE
      DoddRestr.cpp
      Dodd_Restrepo.cpp
      PrestressedSteelMaterial.cpp
      RambergOsgoodSteel.cpp
      StainlessECThermal.cpp
      Steel01.cpp
      Steel01Thermal.cpp
      Steel02.cpp
      Steel02Fatigue.cpp
      Steel02Thermal.cpp
      Steel03.cpp
      Steel2.cpp
      Steel4.cpp
      SteelBRB.cpp
      SteelDRC.cpp
      SteelECThermal.cpp
      SteelFractureDI.cpp
      SteelMP.cpp
      SteelMPF.cpp
    PUBLIC
      DoddRestr.h
      Dodd_Restrepo.h
      PrestressedSteelMaterial.h
      RambergOsgoodSteel.h
      StainlessECThermal.h
      Steel01.h
      Steel01Thermal.h
      Steel02Fatigue.h
      Steel02.h
      Steel02Thermal.h
      Steel03.h
      Steel2.h
      Steel4.h
      SteelBRB.h
      SteelDRC.h
      SteelECThermal.h
      SteelFractureDI.h
      SteelMPF.h
      SteelMP.h
)

