#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================

target_include_directories(OPS_Material PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_sources(OPS_Runtime PRIVATE
  dispatch.cpp
)

target_sources(OPS_Material_f PRIVATE c14-SK-M.f)

target_sources(OPS_Material
  PRIVATE
    Concrete01.cpp
    Concrete01WithSITC.cpp
    Concrete02.cpp
    Concrete02IS.cpp
    Concrete02Thermal.cpp
    Concrete04.cpp
    Concrete05.cpp
    Concrete06.cpp
    Concrete07.cpp
    ConcreteCM.cpp
    ConcreteD.cpp
    ConcreteECThermal.cpp
    ConcreteSakaiKawashima.cpp
    ConcretewBeta.cpp
    ConfinedConcrete01.cpp
    FRCC.cpp
    FRPConfinedConcrete02.cpp
    FRPConfinedConcrete.cpp
    GMG_CyclicReinforcedConcrete.cpp
    SecantConcrete.cpp
    SmoothPSConcrete.cpp
    TDConcrete.cpp
    TDConcreteEXP.cpp
    TDConcreteMC10.cpp
    TDConcreteMC10NL.cpp
    TDConcreteNL.cpp
  PUBLIC
    Concrete01.h
    Concrete01WithSITC.h
    Concrete02.h
    Concrete02IS.h
    Concrete02Thermal.h
    Concrete04.h
    Concrete05.h
    Concrete06.h
    Concrete07.h
    ConcreteCM.h
    ConcreteD.h
    ConcreteECThermal.h
    ConcreteSakaiKawashima.h
    ConcretewBeta.h
    ConfinedConcrete01.h
    FRCC.h
    FRPConfinedConcrete02.h
    FRPConfinedConcrete.h
    GMG_CyclicReinforcedConcrete.h
    SecantConcrete.h
    SmoothPSConcrete.h
    TDConcrete.h
    TDConcreteEXP.h
    TDConcreteMC10.h
    TDConcreteMC10NL.h
    TDConcreteNL.h
)

