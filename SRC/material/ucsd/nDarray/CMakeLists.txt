add_library(UCSD_Array)
target_include_directories(UCSD_Array PUBLIC ${CMAKE_CURRENT_LIST_DIR} ${OPS_Includes})
target_sources(UCSD_Array
    PRIVATE
        basics.cpp
        BJmatrix.cpp
        BJtensor.cpp
        BJvector.cpp
#       Cosseratstraint.cpp
#       Cosseratstresst.cpp
        nDarray.cpp
        straint.cpp
        stresst.cpp
    PUBLIC
        basics.h
        BJmatrix.h
        BJtensor.h
        BJvector.h
#       Cosseratstraint.h
#       Cosseratstresst.h
        nDarray.h
        straint.h
        stresst.h
)

