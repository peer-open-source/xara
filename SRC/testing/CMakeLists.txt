add_executable(test_matrix EXCLUDE_FROM_ALL test_matrix.cpp)
target_link_libraries(test_matrix PRIVATE OpenSeesRT) # G3 OPS_Runtime)
