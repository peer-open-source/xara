
target_sources(OPS_Analysis
    PRIVATE
        CollocationHSFixedNumIter.cpp
        CollocationHSIncrLimit.cpp
        CollocationHSIncrReduct.cpp
        HHTGeneralizedExplicit_TP.cpp
        HHTHSFixedNumIter.cpp
        HHTHSFixedNumIter_TP.cpp
        HHTHSIncrLimit.cpp
        HHTHSIncrLimit_TP.cpp
        HHTHSIncrReduct.cpp
        HHTHSIncrReduct_TP.cpp
        NewmarkHSFixedNumIter.cpp
        NewmarkHSIncrLimit.cpp
        NewmarkHSIncrReduct.cpp
    PUBLIC
        CollocationHSFixedNumIter.h
        CollocationHSIncrLimit.h
        CollocationHSIncrReduct.h
        HHTGeneralizedExplicit_TP.h
        HHTHSFixedNumIter.h
        HHTHSFixedNumIter_TP.h
        HHTHSIncrLimit.h
        HHTHSIncrLimit_TP.h
        HHTHSIncrReduct.h
        HHTHSIncrReduct_TP.h
        NewmarkHSFixedNumIter.h
        NewmarkHSIncrLimit.h
        NewmarkHSIncrReduct.h
)
target_include_directories(OPS_Analysis PUBLIC ${CMAKE_CURRENT_LIST_DIR})

