#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_include_directories(OPS_Domain PUBLIC ${CMAKE_CURRENT_LIST_DIR})

target_sources(OPS_Domain
    PRIVATE
        ConstantSeries.cpp
        LinearSeries.cpp
        PathSeries.cpp
        PathTimeSeries.cpp
        PulseSeries.cpp
        RectangularSeries.cpp
        SimpsonTimeSeriesIntegrator.cpp
        TimeSeries.cpp
        TimeSeriesIntegrator.cpp
        TrapezoidalTimeSeriesIntegrator.cpp
        TriangleSeries.cpp
        TrigSeries.cpp
        PathTimeSeriesThermal.cpp
    PUBLIC
        ConstantSeries.h
        LinearSeries.h
        PathSeries.h
        PathTimeSeries.h
        PulseSeries.h
        RectangularSeries.h
        SimpsonTimeSeriesIntegrator.h
        TimeSeries.h
        TimeSeriesIntegrator.h
        TrapezoidalTimeSeriesIntegrator.h
        TriangleSeries.h
        TrigSeries.h
        PathTimeSeriesThermal.h
)

target_sources(OPS_Domain
    PRIVATE
        PeerMotion.cpp
        PeerNGAMotion.cpp
    PUBLIC
        PeerMotion.h
        PeerNGAMotion.h
)

