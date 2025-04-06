cmake_minimum_required(VERSION 3.5)

# project(arpack C Fortran)

set(arpack_ng_MAJOR_VERSION 3)
set(arpack_ng_MINOR_VERSION 9)
set(arpack_ng_PATCH_VERSION 0)
set(arpack_ng_VERSION ${arpack_ng_MAJOR_VERSION}.${arpack_ng_MINOR_VERSION}.${arpack_ng_PATCH_VERSION})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# set C99 standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)

# Adopted from https://github.com/feymark/arpack.git

if (POLICY CMP0042)
  # enable MACOSX_RPATH by default
  cmake_policy (SET CMP0042 NEW)
endif ()

option(BUILD_SHARED_LIBS "Build shared libraries instead of static libraries" OFF)
option(MPI "Enable parallel support" OFF)
option(ICB "Enable support for *[ae]upd_c with ISO_C_BINDING" ON)
option(PYTHON3 "Enable python3 support" OFF)
set(BOOST_PYTHON_LIBSUFFIX "" CACHE STRING "suffix to add to custom boost python libs")
option(EXAMPLES "Compile ARPACK examples" OFF)
option(TESTS "Compile ARPACK tests" OFF)

# Suffixes: LIBSUFFIX modify ONLY libraries names, ITF64SUFFIX modify BOTH libraries AND include directory names.
set(LIBSUFFIX ""
    CACHE STRING "suffix to add to ARPACK libraries names")
set(ITF64SUFFIX ""
    CACHE STRING "suffix to add to ARPACK include directory and libraries names (use with INTERFACE64)")
set(SYMBOLSUFFIX ""
    CACHE STRING "suffix to add to ARPACK, BLAS and LAPACK function names")
option(INTERFACE64 "use the 64-bit integer interface (ILP64) for ARPACK, BLAS and LAPACK")

# # 'make install' to the correct location, and also define
# # paths for target_include_directories and target_link_libraries
# include(GNUInstallDirs)
# set(ARPACK_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}/arpack${ITF64SUFFIX}")
# set(ARPACK_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/arpackng${LIBSUFFIX}${ITF64SUFFIX}")

#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

function(prefixlist list_name prefix)
    set(${list_name}_TMP)
    foreach(l ${${list_name}})
      list(APPEND ${list_name}_TMP ${prefix}${l} )
    endforeach()
    set(${list_name} "${${list_name}_TMP}" PARENT_SCOPE)
endfunction(prefixlist)

function(examples list_name)
    foreach(l ${${list_name}})
        get_filename_component(lwe ${l} NAME_WE)
        add_executable(${lwe} ${arpackexample_DIR}/${l} ${examples_EXTRA_SRCS})
        target_link_libraries(${lwe} arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
        add_test(NAME "${lwe}_ex" COMMAND ${lwe})
    endforeach()
endfunction(examples)

function(pexamples list_name)
    foreach(l ${${list_name}})
        get_filename_component(lwe ${l} NAME_WE)
        add_executable(${lwe} ${parpackexample_DIR}/${l} )
        target_link_libraries(${lwe} parpack arpack MPI::MPI_Fortran)
        add_test(NAME "${lwe}_ex" COMMAND mpiexec -n 2 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${lwe})
    endforeach()
endfunction(pexamples)


# Enable language(s) before any find_package (in particular before MPI find_package).
if (ICB)
    enable_language(C CXX) # For testing binding with c/c++.

    set(CMAKE_CXX_STANDARD 11) # OK, since cmake-3.1 only.
#   set_property(TARGET tgt PROPERTY CXX_STANDARD 11)

    file(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/PROG_ICB.f90
         "
          PROGRAM PROG_ICB
          USE iso_c_binding
          IMPLICIT NONE
          INTEGER(C_INT) :: a
          a = 1
          END PROGRAM PROG_ICB
         ")
    try_compile(COMPILE_ICB ${CMAKE_BINARY_DIR} ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/PROG_ICB.f90)
    if(NOT ${COMPILE_ICB})
      message(FATAL_ERROR "-- Fortran compiler does not support iso_c_binding.")
    else()
      message("-- Fortran compiler does support iso_c_binding.")
    endif()

else()
    # ICB saves you from old-fashion-boring-cumbersome-fortran/C crap... For arpack symbols (only).

    if (SYMBOLSUFFIX)
        if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
          set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -ffixed-line-length-none")
        elseif ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "LLVMFlang")
          set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -ffixed-line-length=none")
        elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
          set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpp -extend-source")
        else ()
          message(WARNING "build script does not know how to enable your Fortran compiler's preprocessor and support for lines longer than 72 characters: set them manually via FFLAGS.")
        endif ()

        set(SCALARFUNS naitr napps naup2 naupd nconv neigh neupd ngets statn saitr sapps saup2 saupd sconv seigt seupd sgets stats getv0 sortc sortr sesrt stqrb)
        set(COMPLEXFUNS, naitr napps naup2 naupd neigh neupd ngets statn getv0 sortc)

        foreach (f IN LISTS SCALARFUNS)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Ds${f}=s${f}${SYMBOLSUFFIX} -Dd${f}=d${f}${SYMBOLSUFFIX}")
        endforeach ()

        foreach (f IN LISTS COMPLEXFUNS)
            set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Dc${f}=c${f}${SYMBOLSUFFIX} -Dz${f}=z${f}${SYMBOLSUFFIX}")
        endforeach ()

        set(CFUNS snaupd sneupd dnaupd dneupd)
        foreach (f IN LISTS CFUNS)
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D${f}=${f}${SYMBOLSUFFIX}")
            list(APPEND CFUNS_SUFFIXED ${f}${SYMBOLSUFFIX})
        endforeach ()

        include(FortranCInterface)
        FortranCInterface_HEADER(FCMangle.h SYMBOLS ${CFUNS_SUFFIXED})

        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DINCLUDE_FCMANGLE")

        FortranCInterface_VERIFY()
    endif ()
endif ()

# Find MPI

if (MPI)
    if (NOT TARGET MPI::MPI_Fortran) # Search only if not already found by upper CMakeLists.txt
        include(FindMPI)
        find_package(MPI REQUIRED COMPONENTS Fortran)

        # MPI::MPI_* target was already created at this point by FindMPI.cmake if cmake version >= 3.9
        if (NOT TARGET MPI::MPI_Fortran) # Create target "at hand" to ensure compatibility if cmake version < 3.9
            add_library(MPI::MPI_Fortran INTERFACE IMPORTED)
            set_target_properties(MPI::MPI_Fortran PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${MPI_Fortran_INCLUDE_DIRS}")
            set_target_properties(MPI::MPI_Fortran PROPERTIES INTERFACE_LINK_LIBRARIES      "${MPI_Fortran_LIBRARIES}")
        endif()
    endif()

    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAG}")
    if(CMAKE_SYSTEM_NAME MATCHES "Windows" AND CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-invalid-boz")
    endif()

    # Check if we can use ISO_C_BINDING provided by MPI.
    file(WRITE ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/PROG_ICB.f90
         "
          PROGRAM PROG_ICB
          USE :: mpi_f08
          IMPLICIT NONE
          type(MPI_Comm)    comm
          type(MPI_Status)  status
          END PROGRAM PROG_ICB
         ")
         try_compile(COMPILE_ICB ${CMAKE_BINARY_DIR} ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/PROG_ICB.f90 LINK_LIBRARIES MPI::MPI_Fortran)
    if(NOT ${COMPILE_ICB})
      message("-- MPI library does not support iso_c_binding.")
      set(HAVE_MPI_ICB 0)
    else()
      message("-- MPI library does support iso_c_binding.")
      set(HAVE_MPI_ICB 1)
      add_compile_definitions(HAVE_MPI_ICB=1)
    endif()

    # As MPI can be used with or without ISO_C_BINDING (#ifdef), we need to preprocess code before compiling.
    if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp")
    elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpp")
    else ()
      message(WARNING "build script does not know how to preprocess Fortran code: set it manually via FFLAGS.")
    endif ()

    if (ICB)
        if (NOT TARGET MPI::MPI_C) # Search only if not already found by upper CMakeLists.txt
            include(FindMPI)
            find_package(MPI REQUIRED COMPONENTS C)

            if (NOT TARGET MPI::MPI_C) # Create target "at hand" to ensure compatibility if cmake version < 3.9
                add_library(MPI::MPI_C INTERFACE IMPORTED)
                set_target_properties(MPI::MPI_C PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${MPI_C_INCLUDE_DIRS}")
                set_target_properties(MPI::MPI_C PROPERTIES INTERFACE_LINK_LIBRARIES      "${MPI_C_LIBRARIES}")
            endif()
        endif()

        if (NOT TARGET MPI::MPI_CXX) # Search only if not already found by upper CMakeLists.txt
            include(FindMPI)
            find_package(MPI REQUIRED COMPONENTS CXX)

            if (NOT TARGET MPI::MPI_CXX) # Create target "at hand" to ensure compatibility if cmake version < 3.9
                add_library(MPI::MPI_CXX INTERFACE IMPORTED)
                set_target_properties(MPI::MPI_CXX PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${MPI_CXX_INCLUDE_DIRS}")
                set_target_properties(MPI::MPI_CXX PROPERTIES INTERFACE_LINK_LIBRARIES      "${MPI_CXX_LIBRARIES}")
            endif()
        endif()

        include(CheckSymbolExists)
        check_symbol_exists(MPI_Comm_c2f "${MPI_C_INCLUDE_DIRS}/mpi.h" MPI_Comm_c2f_FOUND)
        if(NOT ${MPI_Comm_c2f_FOUND})
          message(FATAL_ERROR "symbol MPI_Comm_c2f does not exist")
        endif()
    endif()
endif()

# Find BLAS

if (NOT TARGET BLAS::BLAS) # Search only if not already found by upper CMakeLists.txt
    find_package(BLAS REQUIRED)

    # BLAS::BLAS target was already created at this point by FindBLAS.cmake if cmake version >= 3.18
    if (NOT TARGET BLAS::BLAS) # Create target "at hand" to ensure compatibility if cmake version < 3.18
        add_library(BLAS::BLAS INTERFACE IMPORTED)
        set_target_properties(BLAS::BLAS PROPERTIES INTERFACE_LINK_LIBRARIES "${BLAS_LIBRARIES}")
    endif()
endif()

# Find LAPACK

if (NOT TARGET LAPACK::LAPACK) # Search only if not already found by upper CMakeLists.txt
    find_package(LAPACK REQUIRED)

    # LAPACK::LAPACK target was already created at this point by FindLAPACK.cmake if cmake version >= 3.18
    if (NOT TARGET LAPACK::LAPACK) # Create target "at hand" to ensure compatibility if cmake version < 3.18
        add_library(LAPACK::LAPACK INTERFACE IMPORTED)
        set_target_properties(LAPACK::LAPACK PROPERTIES INTERFACE_LINK_LIBRARIES "${LAPACK_LIBRARIES}")
    endif()
endif()

# As BLAS/LAPACK does not provide ICB, we may have to deal with symbols the old-fashion-boring-cumbersome-fortran/C way...

if (SYMBOLSUFFIX)
    if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -ffixed-line-length-none")
    elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fpp -extend-source")
    else ()
      message(WARNING "build script does not know how to enable your Fortran compiler's preprocessor and support for lines longer than 72 characters: set them manually via FFLAGS.")
    endif ()

    set(BLASFUNS1 axpy copy gemv geqr2 lacpy lae2 lahqr lanhs larnv lartg lascl laset lasrt scal trevc trmm trsen gbmv gbtrf gbtrs gttrf gttrs pttrf pttrs)
    set(BLASFUNS2 dot ger labad laev2 lamch lanst lanv2 lapy2 larf larfg lasr nrm2 orm2r rot steqr swap)
    set(BLASFUNS3 dotc geru unm2r)
    set(BLASFUNS4 COPY LABAD LAMCH LANHS LANV2 LARFG ROT GEMV)
    set(BLASFUNS5 scnrm2 dznrm2 csscal zdscal)

    foreach (f IN LISTS BLASFUNS1 BLASFUNS2)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Ds${f}=s${f}${SYMBOLSUFFIX} -Dd${f}=d${f}${SYMBOLSUFFIX}")
    endforeach ()

    foreach (f IN LISTS BLASFUNS1 BLASFUNS3)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Dc${f}=c${f}${SYMBOLSUFFIX} -Dz${f}=z${f}${SYMBOLSUFFIX}")
    endforeach ()

    foreach (f IN LISTS BLASFUNS4)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DS${f}=S${f}${SYMBOLSUFFIX} -DD${f}=D${f}${SYMBOLSUFFIX}")
    endforeach ()

    foreach (f IN LISTS BLASFUNS5)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -D${f}=${f}${SYMBOLSUFFIX}")
    endforeach ()

    set(CFUNS sgemm cheev)
    foreach (f IN LISTS CFUNS)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -D${f}=${f}${SYMBOLSUFFIX}")
        list(APPEND CFUNS_SUFFIXED ${f}${SYMBOLSUFFIX})
    endforeach ()

    include(FortranCInterface)
    FortranCInterface_HEADER(FCMangle.h SYMBOLS ${CFUNS_SUFFIXED})

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DINCLUDE_FCMANGLE")

    FortranCInterface_VERIFY()
endif ()

if (MPI)
    set(parpackutil_DIR ${arpack_SOURCE_DIR}/PARPACK/UTIL/)
    set(parpacksrc_DIR ${arpack_SOURCE_DIR}/PARPACK/SRC/)
endif()

file(GLOB arpacksrc_STAT_SRCS ${arpack_SOURCE_DIR}/dbgini.f ${arpack_SOURCE_DIR}/staini.f ${arpack_SOURCE_DIR}/SRC/*.f)

set(arpacksrc_ICB "")
set(parpacksrc_ICB "")
if (ICB)
  file(GLOB arpacksrc_ICB 
       ${arpack_SOURCE_DIR}/SRC/icba*.F90 
       ${arpack_SOURCE_DIR}/ICB/debug_icb.F90 
       ${arpack_SOURCE_DIR}/ICB/stat_icb.F90
  )
  file(GLOB parpacksrc_ICB 
       ${arpack_SOURCE_DIR}/PARPACK/SRC/MPI/icbp*.F90 
       ${arpack_SOURCE_DIR}/ICB/debug_icb.F90 
       ${arpack_SOURCE_DIR}/ICB/stat_icb.F90
  )
endif()

set(arpackutil_STAT_SRCS
    ${arpack_SOURCE_DIR}/UTIL/icnteq.f
    ${arpack_SOURCE_DIR}/UTIL/icopy.f
    ${arpack_SOURCE_DIR}/UTIL/iset.f
    ${arpack_SOURCE_DIR}/UTIL/iswap.f
    ${arpack_SOURCE_DIR}/UTIL/ivout.f
    ${arpack_SOURCE_DIR}/UTIL/second_NONE.f
    ${arpack_SOURCE_DIR}/UTIL/svout.f
    ${arpack_SOURCE_DIR}/UTIL/smout.f
    ${arpack_SOURCE_DIR}/UTIL/dvout.f
    ${arpack_SOURCE_DIR}/UTIL/dmout.f
    ${arpack_SOURCE_DIR}/UTIL/cvout.f
    ${arpack_SOURCE_DIR}/UTIL/cmout.f
    ${arpack_SOURCE_DIR}/UTIL/zvout.f
    ${arpack_SOURCE_DIR}/UTIL/zmout.f )


if (MPI)
    file(GLOB parpackutil_STAT_SRCS ${arpack_SOURCE_DIR}/PARPACK/UTIL/MPI/*.f)
    file(GLOB parpacksrc_STAT_SRCS ${arpack_SOURCE_DIR}/dbgini.f ${arpack_SOURCE_DIR}/staini.f ${arpack_SOURCE_DIR}/PARPACK/SRC/MPI/*.f)
endif()

# use -DBUILD_SHARED_LIBS=ON|OFF to control static/shared
add_library(arpack 
            ${arpackutil_STAT_SRCS}
            ${arpacksrc_STAT_SRCS}
            ${arpacksrc_ICB}
)

target_link_libraries(arpack
  PUBLIC
  $<INSTALL_INTERFACE:$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,,LAPACK::LAPACK>>
  $<INSTALL_INTERFACE:$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,,BLAS::BLAS>>
  $<BUILD_INTERFACE:LAPACK::LAPACK>
  $<BUILD_INTERFACE:BLAS::BLAS>
)

target_link_options(arpack PUBLIC "${EXTRA_LDFLAGS}")
set_target_properties(arpack PROPERTIES OUTPUT_NAME arpack${LIBSUFFIX}${ITF64SUFFIX})
set_target_properties(arpack PROPERTIES VERSION 2.1.0)
set_target_properties(arpack PROPERTIES SOVERSION 2)
target_include_directories(arpack
  PUBLIC
  # Exported location of headers
  $<INSTALL_INTERFACE:${ARPACK_INSTALL_INCLUDEDIR}>
  # Find arpackdef.h, arpackicb.h, stat*.h, debug*.h at build time
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  # For ICB interface
  $<BUILD_INTERFACE:${arpack_SOURCE_DIR}/ICB>
)

if (MPI)
    # use -DBUILD_SHARED_LIBS=ON|OFF to control static/shared
    add_library(parpack ${parpacksrc_STAT_SRCS} ${parpackutil_STAT_SRCS} ${parpacksrc_ICB})
    target_link_libraries(parpack
      PUBLIC
      arpack
      $<INSTALL_INTERFACE:$<IF:$<BOOL:${BUILD_SHARED_LIBS}>,,MPI::MPI_Fortran>>
      $<BUILD_INTERFACE:MPI::MPI_Fortran>
    )
    set_target_properties(parpack PROPERTIES OUTPUT_NAME parpack${LIBSUFFIX}${ITF64SUFFIX})
    set_target_properties(parpack PROPERTIES VERSION 2.1.0)
    set_target_properties(parpack PROPERTIES SOVERSION 2)
endif ()

if (INTERFACE64)
    set(INTERFACE64 1)
    if ("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_Fortran_COMPILER_ID}" STREQUAL "LLVMFlang")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-integer-8")
    elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -i8")
    else ()
      message(WARNING "build script does not know how to make your Fortran compiler use 64-bit integers: set it manually via FFLAGS.")
    endif ()

    # TODO: this needs full re-write of parpack to support ILP64...
    if (MPI)
      message(FATAL_ERROR "Parallel arpack does not support ILP64.")
    endif()
else ()
    set(INTERFACE64 0)
endif ()

include_directories("${PROJECT_BINARY_DIR}") # Find arpackdef.h and arpackicb.h
set(CMAKE_INCLUDE_CURRENT_DIR ON)

#---------------------------
# EXAMPLES
#---------------------------
add_subdirectory(examples)

#---------------------------
# TEST
#---------------------------
function(build_tests)  
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/TESTS)

  add_executable(dnsimp_test TESTS/dnsimp.f TESTS/mmio.f TESTS/debug.h)
  set_target_properties( dnsimp_test PROPERTIES OUTPUT_NAME  dnsimp )
  target_link_libraries(dnsimp_test arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
  add_custom_command(TARGET dnsimp_test POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/TESTS/testA.mtx testA.mtx
  )
  add_test(dnsimp_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/dnsimp)

  if (ICB)
      add_executable(bug_1315_single TESTS/bug_1315_single.c)
      target_include_directories(bug_1315_single PUBLIC ${PROJECT_SOURCE_DIR}/ICB) # Get arpack.h
      target_link_libraries(bug_1315_single arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
      add_test(bug_1315_single_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bug_1315_single)

      add_executable(bug_1315_double TESTS/bug_1315_double.c)
      target_include_directories(bug_1315_double PUBLIC ${PROJECT_SOURCE_DIR}/ICB) # Get arpack.h
      target_link_libraries(bug_1315_double arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
      add_test(bug_1315_double_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bug_1315_double)
  endif()

  add_executable(bug_1323 TESTS/bug_1323.f)
  target_link_libraries(bug_1323 arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
  add_test(bug_1323_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bug_1323)

  add_executable(bug_58_double TESTS/bug_58_double.f)
  target_link_libraries(bug_58_double arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
  add_test(bug_58_double_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bug_58_double)

  add_executable(bug_79_double_complex TESTS/bug_79_double_complex.f)
  target_link_libraries(bug_79_double_complex arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
  add_test(bug_79_double_complex_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bug_79_double_complex)

  add_executable(bug_142 TESTS/bug_142.f)
  target_link_libraries(bug_142 arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
  add_test(bug_142_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bug_142)

  add_executable(bug_142_gen TESTS/bug_142_gen.f)
  target_link_libraries(bug_142_gen arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
  add_test(bug_142_gen_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/bug_142_gen)

  if (MPI)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/PARPACK/TESTS/MPI)

    add_executable(issue46 PARPACK/TESTS/MPI/issue46.f)
    target_link_libraries(issue46 parpack arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
    add_test(issue46_tst mpiexec -n 2 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/issue46)
  endif()

  if (ICB)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/TESTS)

    add_executable(icb_arpack_c TESTS/icb_arpack_c.c)
    target_include_directories(icb_arpack_c PUBLIC ${PROJECT_SOURCE_DIR}/ICB) # Get arpack.h
    target_link_libraries(icb_arpack_c arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
    add_test(icb_arpack_c_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/icb_arpack_c)

    add_executable(icb_arpack_cpp TESTS/icb_arpack_cpp.cpp)
    target_include_directories(icb_arpack_cpp PUBLIC ${PROJECT_SOURCE_DIR}/ICB) # Get arpack.hpp
    target_link_libraries(icb_arpack_cpp arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS})
    add_test(icb_arpack_cpp_tst ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/icb_arpack_cpp)

    if (MPI)
      set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/PARPACK/TESTS/MPI)

      add_executable(icb_parpack_c PARPACK/TESTS/MPI/icb_parpack_c.c)
      target_include_directories(icb_parpack_c PUBLIC ${PROJECT_SOURCE_DIR}/ICB MPI::MPI_C) # Get parpack.h mpi.h
      target_link_libraries(icb_parpack_c parpack arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS} MPI::MPI_C)
      add_test(icb_parpack_c_tst mpiexec -n 2 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/icb_parpack_c)

      add_executable(icb_parpack_cpp PARPACK/TESTS/MPI/icb_parpack_cpp.cpp)
      target_include_directories(icb_parpack_cpp PUBLIC ${PROJECT_SOURCE_DIR}/ICB MPI::MPI_CXX) # Get parpack.hpp mpi.h
      target_link_libraries(icb_parpack_cpp parpack arpack BLAS::BLAS LAPACK::LAPACK ${EXTRA_LDFLAGS} MPI::MPI_CXX)
      add_test(icb_parpack_cpp_tst mpiexec -n 2 ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/icb_parpack_cpp)
    endif()
  endif()
endfunction(build_tests)

if (TESTS)
    enable_testing()
    set(CMAKE_CTEST_COMMAND ctest -V)   
    build_tests()
endif()

############################
# install
############################

# Convert variable names to those expected by the .pc file.
set(prefix ${CMAKE_INSTALL_PREFIX})
set(exec_prefix \${prefix})
set(libdir ${CMAKE_INSTALL_FULL_LIBDIR})
set(includedir ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR})

# Convert (LAPACK|BLAS)_LIBRARIES to -l flags, store in ARPACK_PC_LIBS_PRIVATE.
# LAPACK_LIBRARIES contains transitive deps, no need to parse BLAS_LIBRARIES.
set(ARPACK_PC_LIBS_PRIVATE)
foreach(lib ${LAPACK_LIBRARIES})
    get_filename_component(libname ${lib} NAME)
    string(REGEX REPLACE "^(lib)?([^.]+).*$" "-l\\2" libname ${libname})
    list(APPEND ARPACK_PC_LIBS_PRIVATE "${libname}")
endforeach()
string(REPLACE ";" " " ARPACK_PC_LIBS_PRIVATE "${ARPACK_PC_LIBS_PRIVATE}")

set(PARPACK_PC_LIBS_PRIVATE)
foreach(lib ${LAPACK_LIBRARIES} ${MPI_Fortran_LIBRARIES})
    get_filename_component(libname ${lib} NAME)
    string(REGEX REPLACE "^(lib)?([^.]+).*$" "-l\\2" libname ${libname})
    list(APPEND PARPACK_PC_LIBS_PRIVATE "${libname}")
endforeach()
string(REPLACE ";" " " PARPACK_PC_LIBS_PRIVATE "${PARPACK_PC_LIBS_PRIVATE}")

configure_file(pkg-config/arpack.pc.in "${PROJECT_BINARY_DIR}/arpack${LIBSUFFIX}${ITF64SUFFIX}.pc" @ONLY)
configure_file(pkg-config/parpack.pc.in "${PROJECT_BINARY_DIR}/parpack${LIBSUFFIX}${ITF64SUFFIX}.pc" @ONLY)

if (FALSE)
  install(TARGETS arpack
      EXPORT arpackngTargets
      ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
      LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
      RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
  install(
    FILES 
      "${PROJECT_BINARY_DIR}/arpack${LIBSUFFIX}${ITF64SUFFIX}.pc"
    DESTINATION 
      ${CMAKE_INSTALL_LIBDIR}/pkgconfig
  )
  if (MPI)
    install(TARGETS parpack
        EXPORT arpackngTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
    install(FILES "${PROJECT_BINARY_DIR}/parpack${LIBSUFFIX}${ITF64SUFFIX}.pc"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig)
  endif ()

  if(ICB)
    install(FILES ICB/arpack.h DESTINATION "${ARPACK_INSTALL_INCLUDEDIR}")
    install(FILES ICB/arpack.hpp DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})
    if (MPI)
      install(FILES ICB/parpack.h DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})
      install(FILES ICB/parpack.hpp DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})
    endif()
  endif()
  install(FILES debug.h DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})
  install(FILES debugF90.h DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})
  if(ICB)
    install(FILES ICB/debug_c.h DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})
    install(FILES ICB/debug_c.hpp DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})
  endif()

  install(FILES stat.h DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})
  install(FILES statF90.h DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})

  if (ICB)
    install(FILES ICB/stat_c.h DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})
    install(FILES ICB/stat_c.hpp DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})
  endif()
endif(FALSE)

  configure_file(arpackdef.h.in "${PROJECT_BINARY_DIR}/arpackdef.h" @ONLY)
  configure_file(arpackicb.h.in "${PROJECT_BINARY_DIR}/arpackicb.h" @ONLY)

if (FALSE)
  install(FILES "${PROJECT_BINARY_DIR}/arpackdef.h" DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})

  install(FILES "${PROJECT_BINARY_DIR}/arpackicb.h" DESTINATION ${ARPACK_INSTALL_INCLUDEDIR})

  install(EXPORT arpackngTargets
    DESTINATION "${ARPACK_INSTALL_CMAKEDIR}"
  )

# Provide find_package for arpack-ng to users.
  configure_file(cmake/arpackng-config.cmake.in "${PROJECT_BINARY_DIR}/arpackng-config.cmake" @ONLY)
  configure_file(cmake/arpackng-config-version.cmake.in "${PROJECT_BINARY_DIR}/arpackng-config-version.cmake" @ONLY)
  install(
    FILES
    "${PROJECT_BINARY_DIR}/arpackng-config.cmake"
    "${PROJECT_BINARY_DIR}/arpackng-config-version.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/arpackng${LIBSUFFIX}${ITF64SUFFIX}) # find_package(arpackng)

  configure_file(cmake/tstCMakeInstall.sh.in ${PROJECT_BINARY_DIR}/tstCMakeInstall.sh @ONLY)

  # Packaging: ease arpack-ng distribution (precompiled binaries and sources tarballs).
  set(CPACK_VERSION_MAJOR "${arpack_ng_MAJOR_VERSION}")
  set(CPACK_VERSION_MINOR "${arpack_ng_MINOR_VERSION}")
  set(CPACK_VERSION_PATCH "${arpack_ng_PATCH_VERSION}")
  set(CPACK_SOURCE_IGNORE_FILES "/BUILD/" "/Build/" "/build/" "/local/") # Do not ship (eventual) build or install directories in tarballs.
  set(CPACK_SOURCE_PACKAGE_FILE_NAME "arpack-ng-${CPACK_VERSION_MAJOR}.${CPACK_VERSION_MINOR}.${CPACK_VERSION_PATCH}")
  include(CPack)
endif(FALSE)


function(libsummary title include libraries)
    message("   -- ${title}:")
    foreach(inc ${include})
        message("      -- compile: ${inc}")
    endforeach()
    foreach(lib ${libraries})
        message("      -- link:    ${lib}")
    endforeach()
endfunction(libsummary)

function(cprsummary title compiler debug_flags minsizerel_flags release_flags relwithdebinfo_flags more_flags)
    message("   -- ${title}:      ${compiler}")
    if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug" OR "${CMAKE_BUILD_TYPE}" STREQUAL "DEBUG")
        message("   -- ${title}FLAGS: ${debug_flags} ${more_flags}")
    endif()
    if("${CMAKE_BUILD_TYPE}" STREQUAL "MinSizeRel" OR "${CMAKE_BUILD_TYPE}" STREQUAL "MINSIZEREL")
        message("   -- ${title}FLAGS: ${minsizerel_flags} ${more_flags}")
    endif()
    if("${CMAKE_BUILD_TYPE}" STREQUAL "Release" OR "${CMAKE_BUILD_TYPE}" STREQUAL "RELEASE")
        message("   -- ${title}FLAGS: ${release_flags} ${more_flags}")
    endif()
    if("${CMAKE_BUILD_TYPE}" STREQUAL "RelWithDebInfo" OR "${CMAKE_BUILD_TYPE}" STREQUAL "RELWITHDEBINFO")
        message("   -- ${title}FLAGS: ${relwithdebinfo_flags} ${more_flags}")
    endif()
endfunction(cprsummary)

message("-- Configuration summary for arpack-ng-${arpack_ng_VERSION}:")
message("   -- prefix: ${CMAKE_INSTALL_PREFIX}")
message("   -- MPI: ${MPI} (ICB provided ${HAVE_MPI_ICB})")
message("   -- ICB: ${ICB}")
message("   -- INTERFACE64: ${INTERFACE64}")
cprsummary("FC" "${CMAKE_Fortran_COMPILER}"
                "${CMAKE_Fortran_FLAGS_DEBUG}"
                "${CMAKE_Fortran_FLAGS_MINSIZEREL}"
                "${CMAKE_Fortran_FLAGS_RELEASE}"
                "${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}"
                "${CMAKE_Fortran_FLAGS}")
if (ICB)
    cprsummary("CC" "${CMAKE_C_COMPILER}"
                    "${CMAKE_C_FLAGS_DEBUG}"
                    "${CMAKE_C_FLAGS_MINSIZEREL}"
                    "${CMAKE_C_FLAGS_RELEASE}"
                    "${CMAKE_C_FLAGS_RELWITHDEBINFO}"
                    "${CMAKE_C_FLAGS}")
    cprsummary("CXX" "${CMAKE_CXX_COMPILER}"
                     "${CMAKE_CXX_FLAGS_DEBUG}"
                     "${CMAKE_CXX_FLAGS_MINSIZEREL}"
                     "${CMAKE_CXX_FLAGS_RELEASE}"
                     "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}"
                     "${CMAKE_CXX_FLAGS}")
endif()
if (MPI)
    libsummary("MPIFC" "${MPI_Fortran_INCLUDE_DIRS}" "${MPI_Fortran_LIBRARIES}")
    if (ICB)
        libsummary("MPICC" "${MPI_C_INCLUDE_DIRS}" "${MPI_C_LIBRARIES}")
        libsummary("MPICXX" "${MPI_CXX_INCLUDE_DIRS}" "${MPI_CXX_LIBRARIES}")
    endif()
endif()
libsummary("BLAS" "" "${BLAS_LIBRARIES}")
libsummary("LAPACK" "" "${LAPACK_LIBRARIES}")

