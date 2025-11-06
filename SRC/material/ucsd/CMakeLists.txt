# cmake_minimum_required(VERSION 3.24)
# project(UCSD)

add_subdirectory(nDarray)

add_library(UCSD_FiniteDeformation)
target_include_directories(UCSD_FiniteDeformation 
  PUBLIC "${CMAKE_CURRENT_LIST_DIR}/nDarray/" ${OPS_Includes}
)

target_link_libraries(UCSD_FiniteDeformation
  PUBLIC UCSD_Array OPS_Material
)
add_subdirectory(FiniteDeformation)
add_subdirectory(NewTemplate3Dep)


add_executable(BJmath_tst   tests/BJmath_tst.cpp)
add_executable(elast_tst    tests/elast_tst.cpp)
add_executable(el-pl-test   tests/el-pl-test.cpp)
add_executable(math_tst     tests/math_tst.cpp)

target_link_libraries(BJmath_tst   PRIVATE UCSD_Array OpenSeesRT)
target_link_libraries(elast_tst    PRIVATE UCSD_Array OpenSeesRT)
target_link_libraries(el-pl-test   PRIVATE UCSD_Array OpenSeesRT)
target_link_libraries(math_tst     PRIVATE UCSD_Array OpenSeesRT)

