#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
target_link_libraries(OPS_Database PRIVATE OPS_Logging)
target_sources(OPS_Database
    PRIVATE
        FE_Datastore.cpp
        FileDatastore.cpp
#       MySqlDatastore.cpp
#       OracleDatastore.cpp
#       BerkeleyDbDatastore.cpp
    PUBLIC
        FE_Datastore.h
        FileDatastore.h
#       MySqlDatastore.h
#       OracleDatastore.h
#       BerkeleyDbDatastore.h
)
target_include_directories(OPS_Database PUBLIC ${CMAKE_CURRENT_LIST_DIR})

