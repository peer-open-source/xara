#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
# Claudio Perez
#==============================================================================

cmake_minimum_required(VERSION 3.16)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(OpenSees 
        DESCRIPTION "Open System For Earthquake Engineering Simulation"
        LANGUAGES   C CXX Fortran)


option(Xara_MumpsSV_Sequential
  "Include sequential Mumps solver" OFF)

set(OPS_BUNDLED_DIR   "${CMAKE_CURRENT_LIST_DIR}/OTHER/")
set(OPS_SRC_DIR       "${CMAKE_CURRENT_LIST_DIR}/SRC/")
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}"
                      "${CMAKE_CURRENT_LIST_DIR}/tools/cmake")

enable_testing()


#==============================================================================
#                            OS Configuration
#==============================================================================
add_library(OPS_OS_Specific_libs INTERFACE)

if (UNIX)
  if (DEFINED Dependencies)
    include("OpenSeesDependencies${Dependencies}")
  else ()
    include(OpenSeesDependenciesUnix)
  endif ()

  # General debug / optimization options
  add_compile_options(
    $<$<CONFIG:DEBUG>:-g3>
    $<$<CONFIG:DEBUG>:-Og>
    $<$<CONFIG:DEBUG>:-O1>
  )
  if (FALSE) #DEFINED ENV{OPENSEESRT_BUILD} AND  STREQUAL "light")
    # For CI environments; minimize optimization and debug info
    # for efficiency
    add_compile_options(-O2 -g0)
  elseif (CMAKE_BUILD_TYPE STREQUAL "RELEASE")
    add_compile_options(-O3 -fno-stack-protector)
  endif()

  if (DEFINED ProfileBuild)
    # Profiling with gprof
    add_compile_options(-fno-omit-frame-pointer -g -p)
  endif()

  if (APPLE)
    message(":: MacOS")
    add_compile_definitions(_LINUX _UNIX)
    add_link_options("-Wl,-no_compact_unwind")
    set(CMAKE_OSX_DEPLOYMENT_TARGET 15.0)

    if (CMAKE_BUILD_TYPE STREQUAL "DEBUG")
      add_compile_options(
        -g
        -fno-omit-frame-pointer
#       -gdwarf-4
        # -fsanitize=address,undefined # , leak
      )
      # add_link_options(
      #   -fsanitize=address,undefined # , leak
      #   -shared-libsan
      # )
    else()
      add_compile_options(
        $<$<COMPILE_LANGUAGE:CXX>:-O3>
        -fno-stack-protector
        $<$<COMPILE_LANGUAGE:CXX>:-fno-math-errno>
        # -freciprocal-math
        # $<$<COMPILE_LANGUAGE:CXX>:-fno-trapping-math>
      )
    endif()

  else() # LINUX
    message(":: LINUX")
    add_compile_definitions(_LINUX _UNIX USE_TCL_STUBS _TCL85)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmax-errors=1 -ffloat-store")
    if (CMAKE_BUILD_TYPE STREQUAL "DEBUG")
      add_compile_options(
        $<$<CONFIG:DEBUG>:-ggdb>
        $<$<COMPILE_LANG_AND_ID:CXX,GNU>:-fsanitize=address,undefined,leak> #, leak
        $<$<COMPILE_LANG_AND_ID:CXX,GNU>:-fno-omit-frame-pointer>
      )
      add_link_options(
        $<$<LINK_LANG_AND_ID:CXX,GNU>:-fsanitize=address,undefined,leak> # , leak
      )
    elseif (CMAKE_BUILD_TYPE STREQUAL "RELEASE")
      if (CMAKE_SYSTEM_PROCESSOR MATCHES "(x86)|(X86)|(amd64)|(AMD64)")
        # Some x86_64 optimizations
        add_compile_options(
          $<$<CONFIG:RELEASE>:-march=haswell>
          $<$<CONFIG:RELEASE>:-msse4>
          $<$<CONFIG:RELEASE>:-mtune=native>
          $<$<COMPILE_LANG_AND_ID:CXX,GNU>:-mavx2>
        )
      endif()
    endif()
  endif()

elseif (WIN32) # NOTE: this will execute for both 32-bit and 64-bit builds.
  message(STATUS ":: Windows")
  include(OpenSeesDependenciesWin)
  add_compile_definitions(_WIN32 _TCL85)
  # Optimize for Architecture Specifics
  add_compile_options($<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/favor:INTEL64>)
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>DLL")
endif()

if (CMAKE_BUILD_TYPE STREQUAL "RELEASE")
  add_compile_definitions(XARA_RELEASE)
endif()

# Check for IPO and maybe use it; this can significantly slow down linking 
# and can be problematic on slow or memory constrained machines (e.g., my laptop)
include(CheckIPOSupported)
check_ipo_supported(RESULT has_ipo OUTPUT error)
if (NOT has_ipo)
  set(CMAKE_INTERPROCEDURAL_OPTIMIZATION FALSE)
endif ()

#==============================================================================
#                            General Setup
#
#==============================================================================
#----------------------------------------------------------------
# Compilers
#---------------------------------------------------------------- 

# Fortran
#--------------------------------------
enable_language(Fortran)
set(CMAKE_FORTRAN_STANDARD 08)

# C/C++
#--------------------------------------

set(CMAKE_C_STANDARD 11)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_COLOR_DIAGNOSTICS ON)
set(CMAKE_C_STANDARD_REQUIRED True)
set(CMAKE_CXX_STANDARD_REQUIRED True)

#
# Global compiler definitions
#
add_compile_definitions(
  "$<$<COMPILE_LANGUAGE:C,CXX>: OPS_USE_RUNTIME>"
  "$<$<COMPILE_LANGUAGE:C,CXX>: USE_TCL_STUBS>"
  "$<$<COMPILE_LANGUAGE:C,CXX>: _TCL85>"
) # _THREADS

#
# Global compiler warnings
#
if (NOT CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  add_compile_options(
    "$<$<COMPILE_LANGUAGE:CXX>:-Wall>"
    "$<$<COMPILE_LANGUAGE:CXX>:-Wno-missing-braces>"
    "$<$<COMPILE_LANGUAGE:CXX>:-Wno-unused-parameter>"
    "$<$<COMPILE_LANGUAGE:CXX>:-Wno-deprecated-declarations>"
    "$<$<COMPILE_LANGUAGE:CXX>:-Wno-unused-parameter>"
  )
  # if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    add_compile_options(
      "$<$<COMPILE_LANGUAGE:CXX>:-Wno-reorder-ctor>"
      "$<$<COMPILE_LANGUAGE:CXX>:-Wno-inconsistent-missing-override>"
      "$<$<COMPILE_LANG_AND_ID:CXX,AppleClang>:-Wno-reorder>"
      "$<$<COMPILE_LANG_AND_ID:CXX,AppleClang>:-Wno-deprecated-declarations>"
      "$<$<COMPILE_LANG_AND_ID:CXX,AppleClang>:-Wno-ignored-optimization-argument>"
      "$<$<COMPILE_LANG_AND_ID:CXX,AppleClang>:-Wno-invalid-utf8>"
    )
  # else()
    add_compile_options(
      "$<$<COMPILE_LANG_AND_ID:CXX,GNU>:-Wno-reorder>"
      "$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-Wno-unused-dummy-argument>"
    )
  # endif()
endif()

function(disable_rtti)
    cmake_parse_arguments(ARG "PUBLIC" "TARGET;FILE" "" ${ARGN})

    if(ARG_TARGET AND ARG_FILE)
        message(FATAL_ERROR "disable_rtti: TARGET and FILE are mutually exclusive")
    elseif(NOT ARG_TARGET AND NOT ARG_FILE)
        message(FATAL_ERROR "disable_rtti: must specify TARGET or FILE")
    endif()

    set(rtti_flags
        $<$<CXX_COMPILER_ID:MSVC>:/GR->
        $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:-fno-rtti>
    )

    if(ARG_TARGET)
        if(ARG_PUBLIC)
            set(scope PUBLIC)
        else()
            set(scope PRIVATE)
        endif()
        target_compile_options(${ARG_TARGET} ${scope} ${rtti_flags})
    else()
        set_source_files_properties(${ARG_FILE}
            PROPERTIES COMPILE_OPTIONS "${rtti_flags}"
        )
    endif()
endfunction()

#----------------------------------------------------------------
# Global Includes
#---------------------------------------------------------------- 
set(OPS_Includes
  ${OPS_SRC_DIR}
  ${OPS_SRC_DIR}/api
  ${OPS_SRC_DIR}/handler
  ${OPS_SRC_DIR}/actor/actor
  ${OPS_SRC_DIR}/actor/objectBroker
  ${OPS_SRC_DIR}/actor/machineBroker
  ${OPS_SRC_DIR}/actor/channel
  ${OPS_SRC_DIR}/actor/message
  ${OPS_SRC_DIR}/tagged
  ${OPS_SRC_DIR}/tagged/storage

  ${OPS_SRC_DIR}/database

  ${OPS_SRC_DIR}/domain/constraints
  ${OPS_SRC_DIR}/domain/component
  ${OPS_SRC_DIR}/domain/domain
  ${OPS_SRC_DIR}/domain/domain/single
  ${OPS_SRC_DIR}/domain/groundMotion
  ${OPS_SRC_DIR}/domain/load
  ${OPS_SRC_DIR}/domain/node
  ${OPS_SRC_DIR}/domain/pattern
  ${OPS_SRC_DIR}/domain/pattern/drm
  ${OPS_SRC_DIR}/domain/region
  ${OPS_SRC_DIR}/domain/subdomain

  ${OPS_SRC_DIR}/recorder

  ${OPS_SRC_DIR}/element
  ${OPS_SRC_DIR}/graph/graph
  ${OPS_SRC_DIR}/graph/numberer
  ${OPS_SRC_DIR}/graph/partitioner
  ${OPS_SRC_DIR}/material
  ${OPS_SRC_DIR}/material/nD
  ${OPS_SRC_DIR}/material/section
  ${OPS_SRC_DIR}/material/section/yieldSurface
  ${OPS_SRC_DIR}/material/uniaxial
  ${OPS_SRC_DIR}/matrix
  ${OPS_SRC_DIR}/system_of_eqn
  ${OPS_SRC_DIR}/system_of_eqn/linearSOE
  ${OPS_SRC_DIR}/system_of_eqn/eigenSOE

  ${OPS_SRC_DIR}/analysis/algorithm
  ${OPS_SRC_DIR}/analysis/dof_grp
  ${OPS_SRC_DIR}/analysis/fe_ele
  ${OPS_SRC_DIR}/analysis/handler
  ${OPS_SRC_DIR}/analysis/integrator
)


#==============================================================================
#                            Define Targets
#
#==============================================================================
# find_package(MKL)
set(BLA_STATIC ON)
add_subdirectory(OTHER)

find_package(TCL)
find_package(TclStub) # TODO: may not need to find TCL first
set(TCL_INCLUDE_PATH ${TCL_INCLUDE_DIRS})

include_directories(${TCL_INCLUDE_PATH})

# Print out a summary of build configuration
message("\n:: OpenSees Summary")
message("  :: CXX:       ${CMAKE_CXX_FLAGS}")
message("  :: BLAS:      ${BLAS_LIBRARIES}")
message("  :: LAPACK:    ${LAPACK_LIBRARIES}")
message("  :: Eigen:     ${Eigen_INCLUDE_DIRS}")
message("  :: SUPERLU:   ${SUPERLU_LIBRARIES}")
message("  :: ARPACK:    ${ARPACK_LIBRARIES}")
message("  :: UMFPACK:   ${UMFPACK_LIBRARIES}")
message("  :: AMD:       ${AMD_LIBRARIES}")
message("  :: TCL  a:    ${TCL_LIBRARY}  ")
message("          b:    ${TCL_STUB_LIBRARY}")
message("          c:    ${TCL_LIBRARIES}")
message("")
message("  :: IPO        ${CMAKE_INTERPROCEDURAL_OPTIMIZATION}")
message("  :: Arch       ${CMAKE_SYSTEM_PROCESSOR}")

# TODO: Probably not needed, link libs as needed
# in sys_of_eqn/linearSOE/*/CMakeLists.txt
add_library(OPS_Numerics           INTERFACE)
target_link_libraries(OPS_Numerics INTERFACE
  ${ARPACK_LIBRARIES} 
  ${SUPERLU_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES}
)

add_library(FastMath INTERFACE)
add_library(Xara::FastMath ALIAS FastMath)
target_compile_options(FastMath INTERFACE
  $<$<COMPILE_LANG_AND_ID:CXX,GNU>:-fno-math-errno>
  $<$<COMPILE_LANG_AND_ID:CXX,AppleClang>:-fno-math-errno>
  $<$<COMPILE_LANG_AND_ID:CXX,Clang>:-fno-math-errno>

  $<$<COMPILE_LANG_AND_ID:CXX,GNU>:-fno-signaling-nans>

  $<$<COMPILE_LANG_AND_ID:CXX,GNU>:-fno-trapping-math>
  $<$<COMPILE_LANG_AND_ID:CXX,AppleClang>:-fno-trapping-math>
  $<$<COMPILE_LANG_AND_ID:CXX,Clang>:-fno-trapping-math>
  $<$<CXX_COMPILER_ID:MSVC>:/fp:except->

  $<$<COMPILE_LANG_AND_ID:CXX,GNU>:-freciprocal-math>
  $<$<COMPILE_LANG_AND_ID:CXX,AppleClang>:-freciprocal-math>
  $<$<COMPILE_LANG_AND_ID:CXX,Clang>:-freciprocal-math>

  # $<$<COMPILE_LANG_AND_ID:CXX,GNU>:-ffp-contract=fast>
  # $<$<COMPILE_LANG_AND_ID:CXX,AppleClang>:-ffp-contract=fast>
  # $<$<COMPILE_LANG_AND_ID:CXX,Clang>:-ffp-contract=fast>

  $<$<COMPILE_LANG_AND_ID:CXX,GNU>:-fno-signed-zeros>
  $<$<COMPILE_LANG_AND_ID:CXX,AppleClang>:-fno-signed-zeros>
  $<$<COMPILE_LANG_AND_ID:CXX,Clang>:-fno-signed-zeros>
)

# Core OpenSees
add_library(OPS_Matrix             OBJECT)
# disable_rtti(TARGET OPS_Matrix)
add_library(OPS_Handler            OBJECT)
add_library(OPS_Logging            OBJECT)
add_library(OPS_Actor              OBJECT)
add_library(G3_ObjectBroker        OBJECT)
add_library(OPS_Tagged             OBJECT)
add_library(OPS_Utilities          OBJECT)
add_library(OPS_Domain             OBJECT)
add_library(OPS_Recorder           OBJECT)

add_library(OPS_ConvergenceTest    OBJECT)
add_library(OPS_Thermal            OBJECT)
add_library(OPS_ElementFortran     OBJECT)
add_library(OPS_Material           OBJECT)
add_library(OPS_MaterialFortran    OBJECT)
add_library(OPS_Transform "")
add_library(OPS_Element            OBJECT)

add_library(OPS_SysOfEqn           OBJECT)
# add_library(OPS_SysOfEqnMP         OBJECT)
add_library(OPS_Analysis           OBJECT)
add_library(OPS_Damage             OBJECT)
# add_library(OPS_Runtime            OBJECT)

add_library(Xara_InterpreterCore   OBJECT)
add_library(Xara_InterpreterRT     OBJECT)
add_library(Xara_InterpreterMP     OBJECT)
add_library(XaraCmdsRT             OBJECT)
add_library(XaraCmdsMP             OBJECT)
add_library(XaraCmds_Analysis      OBJECT)
add_library(XaraCmds_AnalysisRT    OBJECT)
add_library(XaraCmds_AnalysisMP    OBJECT)
add_library(XaraCmds_Modeling      OBJECT)
add_library(XaraCmds_General       OBJECT)

add_library(OPS_Database           OBJECT EXCLUDE_FROM_ALL)

# Optional Extensions
add_library(OPS_Parallel           OBJECT EXCLUDE_FROM_ALL)
add_library(OPS_Partition          OBJECT EXCLUDE_FROM_ALL)

# Packaged libraries
add_library(G3)

# set_property(TARGET OPS_Runtime           PROPERTY POSITION_INDEPENDENT_CODE 1)
set_property(TARGET OPS_MaterialFortran   PROPERTY POSITION_INDEPENDENT_CODE 1)
set_property(TARGET OPS_ElementFortran    PROPERTY POSITION_INDEPENDENT_CODE 1)
set_property(TARGET G3                    PROPERTY POSITION_INDEPENDENT_CODE 1)
target_include_directories(G3 PUBLIC 
  ${OPS_SRC_DIR}/handler
  ${OPS_SRC_DIR}/actor/actor
  ${OPS_SRC_DIR}/actor/objectBroker
  ${OPS_SRC_DIR}/actor/channel
  ${OPS_SRC_DIR}/actor/message
  ${OPS_Includes}
)

target_include_directories(OPS_Domain PRIVATE
  ${OPS_SRC_DIR}/recorder/response
  ${OPS_SRC_DIR}/material/uniaxial
  ${OPS_SRC_DIR}/material/elastic
  ${OPS_SRC_DIR}/material/nD
  ${OPS_SRC_DIR}/tagged
  ${OPS_SRC_DIR}/tagged/storage
  ${OPS_Includes}
)
target_include_directories(OPS_Element PRIVATE
  ${OPS_SRC_DIR}/tagged
  ${OPS_SRC_DIR}/runtime/runtime/legacy/Renderer
  ${OPS_Includes}
)
target_include_directories(OPS_Recorder PRIVATE
  ${OPS_SRC_DIR}/tagged
  ${OPS_SRC_DIR}/recorder/response
  ${OPS_Includes}
)
target_include_directories(OPS_Transform PRIVATE
  ${OPS_SRC_DIR}/recorder/response
  ${OPS_Includes}
)
target_include_directories(OPS_Damage PRIVATE
  ${OPS_SRC_DIR}/recorder/response
  ${OPS_Includes}
)

# Add sources to targets
add_subdirectory(${OPS_SRC_DIR})

#==============================================================================
#                            Configure targets
#
#==============================================================================
target_link_libraries(OPS_Element   PRIVATE ${OPS_Element_List})
target_link_libraries(OPS_Material  PRIVATE ${TCL_STUB_LIBRARY})


# G3
#----------------------------
target_link_libraries(G3
  PUBLIC
    OPS_Database
    OPS_Matrix
    OPS_Actor
    OPS_Analysis 
    OPS_Domain
    OPS_Thermal
    OPS_ConvergenceTest
    OPS_Element
    OPS_ElementFortran
    OPS_Material
    OPS_MaterialFortran
    OPS_Material_YieldSurface
    OPS_Section_Repres
    OPS_Section_YieldSurface
    OPS_Recorder
    OPS_Handler
    OPS_SysOfEqn
    OPS_Tagged
    OPS_Utilities
    graph
    OPS_Transform
    G3_ObjectBroker
    OPS_Numerics 
    ${OPS_Element_List}
    OPS_Material
    OPS_Damage
    G3_ObjectBroker
    OPS_Algorithm 
)


#==============================================================================
#                            Apply Options
#
#==============================================================================

# Extensions
#----------------------------
message(":: Configuring OpenSees extensions")
foreach(extension IN LISTS OPS_Element_List)
  string(TOUPPER "${extension}" ext_flag) 
  string(REGEX REPLACE "^OPS_" "OPSDEF_" ext_flag "${ext_flag}")
  add_compile_definitions(${ext_flag})
endforeach()

foreach(extension IN LISTS OPS_Exclude_List)
  string(TOUPPER "${extension}" ext_flag) 
  string(REGEX REPLACE "^OPS_" "OPS_EXCLUDE_" ext_flag "${ext_flag}")
  message("    Adding macro definition '${ext_flag}'")
  add_compile_definitions(${ext_flag})
endforeach()


# HDF5
#----------------------------
if (FALSE)
  find_package(HDF5)
  if(HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIR})
    set(_hdf5_libs hdf5 hdf5_cpp)
    add_compile_definitions(-D_H5DRM)
  else()
      message(STATUS "OPS >>> Could not find HDF5")
  endif()
endif()


# Include test suite
#add_subdirectory(EXAMPLES/)

get_target_property(OPS_Damage_COMPILE_OPTIONS OPS_Damage COMPILE_OPTIONS)
  string(REPLACE "-Wall" "" OPS_Damage_COMPILE_OPTIONS "${OPS_Damage_COMPILE_OPTIONS}")
  string(REPLACE "-Wextra" "" OPS_Damage_COMPILE_OPTIONS "${OPS_Damage_COMPILE_OPTIONS}")
  set_target_properties(OPS_Damage PROPERTIES COMPILE_OPTIONS "${OPS_Damage_COMPILE_OPTIONS}")

get_target_property(OPS_Material_COMPILE_OPTIONS OPS_Material COMPILE_OPTIONS)
  string(REPLACE "-Wall" "" OPS_Material_COMPILE_OPTIONS "${OPS_Material_COMPILE_OPTIONS}")
  string(REPLACE "-Wextra" "" OPS_Material_COMPILE_OPTIONS "${OPS_Material_COMPILE_OPTIONS}")
  set_target_properties(OPS_Material PROPERTIES COMPILE_OPTIONS "${OPS_Material_COMPILE_OPTIONS}")

get_target_property(OPS_MaterialFortran_COMPILE_OPTIONS OPS_MaterialFortran COMPILE_OPTIONS)
  string(REPLACE "-Wall" "" OPS_MaterialFortran_COMPILE_OPTIONS   "${OPS_MaterialFortran_COMPILE_OPTIONS}")
  string(REPLACE "-Wextra" "" OPS_MaterialFortran_COMPILE_OPTIONS "${OPS_MaterialFortran_COMPILE_OPTIONS}")
  set_target_properties(OPS_MaterialFortran PROPERTIES COMPILE_OPTIONS   "${OPS_MaterialFortran_COMPILE_OPTIONS}")
  target_compile_options(OPS_MaterialFortran PRIVATE
    "$<$<COMPILE_LANG_AND_ID:Fortran,GNU>:-Wno-unused-variable>"
  )

get_target_property(OPS_Element_COMPILE_OPTIONS OPS_Element COMPILE_OPTIONS)
  string(REPLACE "-Wall" "" OPS_Element_COMPILE_OPTIONS "${OPS_Element_COMPILE_OPTIONS}")
  string(REPLACE "-Wextra" "" OPS_Element_COMPILE_OPTIONS "${OPS_Element_COMPILE_OPTIONS}")
  set_target_properties(OPS_Element PROPERTIES COMPILE_OPTIONS "${OPS_Element_COMPILE_OPTIONS}")

if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
  target_compile_options(OPS_ElementFortran PRIVATE "-w")
endif()

#
# Install
#

set(install_targets "OpenSeesRT")
if (TARGET OpenSeesPyRT)
    list(APPEND install_targets OpenSeesPyRT)
endif()
if (TARGET LibOpenSeesMP)
    list(APPEND install_targets LibOpenSeesMP)
endif()
if (TARGET LibOpenSeesSP)
    list(APPEND install_targets LibOpenSeesSP)
endif()

install(
    TARGETS ${install_targets} 
    LIBRARY DESTINATION .
)
