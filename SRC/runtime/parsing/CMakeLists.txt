
target_sources(OPS_Runtime PRIVATE
  InterpreterAPI.cpp
)

add_library(Xara_Parsing OBJECT)
target_sources(Xara_Parsing PRIVATE ParseFloats.cpp)
target_include_directories(Xara_Parsing PUBLIC ${CMAKE_CURRENT_LIST_DIR})
