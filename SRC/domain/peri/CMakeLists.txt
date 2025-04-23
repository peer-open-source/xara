target_sources(OPS_Runtime PRIVATE commands.cpp)

target_sources(OPS_Domain
  PRIVATE
    PeriDomainBase.cpp
  PUBLIC
    PeriDomain.h
)
target_include_directories(OPS_Domain PUBLIC ${CMAKE_CURRENT_LIST_DIR})
target_include_directories(OPS_Domain PUBLIC 
  ${CMAKE_CURRENT_LIST_DIR}/nosb
  ${CMAKE_CURRENT_LIST_DIR}/mate
)
