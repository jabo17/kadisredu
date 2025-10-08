cmake_minimum_required(VERSION 3.25)
# project setup for  'KaDisRedu'
project(
        kadisredu
        DESCRIPTION "KaDisRedu - A distributed kernelization algorithm for the Maximum Weight Independent Set Problem."
        LANGUAGES C CXX

        VERSION 0.0.1
)

option(KADISREDU_BUILD_OFFLINE OFF)

set(FETCHCONTENT_BASE_DIR "${CMAKE_CURRENT_LIST_DIR}/_deps")
if (KADISREDU_BUILD_OFFLINE)
  set(FETCHCONTENT_FULLY_DISCONNECTED TRUE)
endif()
include(FetchContent)

include(CheckCXXCompilerFlag)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

set(CMAKE_CXX_STANDARD 20)
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    add_compile_options(-O1)
    add_compile_options(-g)
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-O3)
    add_definitions(-DNDEBUG)
    CHECK_CXX_COMPILER_FLAG(-fno-omit-frame-pointer COMPILER_SUPPORTS_OMIT_FRAME_POINTER)

    CHECK_CXX_COMPILER_FLAG(-march=native COMPILER_SUPPORTS_MARCH_NATIVE)
    if(COMPILER_SUPPORTS_MARCH_NATIVE)
        add_compile_options(-march=native)
    endif()

    CHECK_CXX_COMPILER_FLAG(-funroll-loops COMPILER_SUPPORTS_FUNROLL_LOOPS)
    if(COMPILER_SUPPORTS_FUNROLL_LOOPS)
        add_compile_options(-funroll-loops)
    endif()
endif ()

# project setup for  'KaDisRedu'
project(
        kadisredu
        DESCRIPTION "KaDisRedu - A distributed kernelization algorithm for the Maximum Weight Independent Set Problem."
        LANGUAGES CXX

        VERSION 0.0.1
)

# options
################
option(KADISREDU_BUILD_EXECUTABLES "build executables for KaDisRedu" ON)

if (KADISREDU_BUILD_EXECUTABLES)
    FetchContent_Declare(
            CLI11
            GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
            GIT_TAG v2.3.2
            SYSTEM
    )
    FetchContent_MakeAvailable(CLI11)
    add_subdirectory(extern/json)


    FetchContent_Declare(
      tomlplusplus
      GIT_REPOSITORY https://github.com/marzer/tomlplusplus.git
      GIT_TAG        v3.4.0 # or latest
    )
    FetchContent_MakeAvailable(tomlplusplus)
endif()


if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(KASSERT_ASSERTION_LEVEL 30)
    set(KAMINPAR_ASSERTION_LEVEL "normal")
    set(KAMPING_ASSERTION_LEVEL "normal")
elseif(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(KASSERT_ASSERTION_LEVEL 10)
    set(KAMINPAR_ASSERTION_LEVEL "light")
    set(KAMPING_ASSERTION_LEVEL "light")
endif ()

set(INSTALL_KAGEN OFF)
add_subdirectory(extern/KaGen)


# KaMinPar
################
# KaMinPar uses version of kassert that does not recognizes redefinitions of its targets,
#  but kassert in kamping does. Therefore, add KaMinPar before kamping
set(KAMINPAR_BUILD_DISTRIBUTED ON)
# Disable the building of apps/benchmarks to prevent KaMinPar to include is older dependency of KaGen (used for CLI)
# We include our own version of KaGen later for our CLI applications
set(KAMINPAR_BUILD_APPS OFF)
set(KAMINPAR_BUILD_BENCHMARKS OFF)
set(INSTALL_KAMINPAR OFF)
add_subdirectory(extern/KaMinPar)
#dd_library(growt::growt ALIAS growt)

# KaMPIng
################
FetchContent_Declare(
        kamping
        GIT_REPOSITORY https://github.com/kamping-site/kamping.git
        GIT_TAG e2c5e6e
        SYSTEM
)
FetchContent_MakeAvailable(kamping)
#add_subdirectory(extern/kamping)

# malloc_count
################
option(USE_MALLOC_COUNT "Use malloc count (override mallocs)" OFF)
FetchContent_Declare(
        malloc_count
        GIT_REPOSITORY https://github.com/bingmann/malloc_count.git
        GIT_TAG ddd0565
        SYSTEM
)
if(USE_MALLOC_COUNT)
    FetchContent_MakeAvailable(malloc_count)
    FetchContent_GetProperties(malloc_count)
    add_library(malloc_count ${malloc_count_SOURCE_DIR}/malloc_count.c)
    target_include_directories(malloc_count SYSTEM PUBLIC ${malloc_count_SOURCE_DIR})
    target_link_libraries(malloc_count PUBLIC dl)
endif()

# Abseil
################
FetchContent_Declare(
        absl
        GIT_REPOSITORY https://github.com/abseil/abseil-cpp.git
        GIT_TAG 4447c75
        SYSTEM
)
set(ABSL_PROPAGATE_CXX_STD ON)
FetchContent_MakeAvailable(absl)

# String formatting
################
FetchContent_Declare(
        fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG 10.1.1
        SYSTEM
)
FetchContent_MakeAvailable(fmt)

# message-queue
add_subdirectory(extern/message-queue)

# set KaMIS/wmis path for red2pack-lib and includes for app/
set(KAMIS_WMIS ${CMAKE_CURRENT_SOURCE_DIR}/extern/KaMIS/wmis)

# includes KaMIS
set(KAMIS_WMIS_INCLUDES
        ${KAMIS_WMIS}/extern/argtable3-3.0.3
        ${KAMIS_WMIS}/lib
        ${KAMIS_WMIS}/lib/data_structure
        ${KAMIS_WMIS}/lib/data_structure/priority_queues
        ${KAMIS_WMIS}/lib/mis
        ${KAMIS_WMIS}/lib/mis/ils
        ${KAMIS_WMIS}/lib/mis/initial_mis
        ${KAMIS_WMIS}/lib/mis/kernel
        ${KAMIS_WMIS}/lib/tools

)
set(KAMIS_WMIS_KAHIP_INCLUDES
        ${KAMIS_WMIS}/extern/KaHIP
        ${KAMIS_WMIS}/extern/KaHIP/app
        ${KAMIS_WMIS}/extern/KaHIP/interface
        ${KAMIS_WMIS}/extern/KaHIP/lib
        ${KAMIS_WMIS}/extern/KaHIP/lib/algorithms
        ${KAMIS_WMIS}/extern/KaHIP/lib/data_structure
        ${KAMIS_WMIS}/extern/KaHIP/lib/data_structure/matrix
        ${KAMIS_WMIS}/extern/KaHIP/lib/data_structure/priority_queues
        ${KAMIS_WMIS}/extern/KaHIP/lib/io
        ${KAMIS_WMIS}/extern/KaHIP/lib/parallel_mh
        ${KAMIS_WMIS}/extern/KaHIP/lib/parallel_mh/exchange
        ${KAMIS_WMIS}/extern/KaHIP/lib/parallel_mh/galinier_combine
        ${KAMIS_WMIS}/extern/KaHIP/lib/partition
        ${KAMIS_WMIS}/extern/KaHIP/lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement
        ${KAMIS_WMIS}/extern/KaHIP/lib/tools
)
# add KaMIS/wmis
add_subdirectory(${KAMIS_WMIS})
# this is the original KaHIP version used by KaMIS
add_library(lib_kamis_kahip STATIC $<TARGET_OBJECTS:libkaffpa2>)
target_include_directories(lib_kamis_kahip PUBLIC ${KAMIS_WMIS_KAHIP_INCLUDES})

# KaMIS only exposes its own header files
add_library(lib_kamis_wmis STATIC $<TARGET_OBJECTS:libsources>)
target_include_directories(lib_kamis_wmis PUBLIC ${KAMIS_WMIS_INCLUDES})
target_link_libraries(lib_kamis_wmis PRIVATE lib_kamis_kahip)
add_library(kamis::wmis ALIAS lib_kamis_wmis)

# KaDisRedu lib.
################

add_subdirectory(lib)

add_library(KaDisRedu::KaDisRedu ALIAS kadisredu)

# KaDisRedu exe.
################
if (KADISREDU_BUILD_EXECUTABLES)
    add_subdirectory(apps)
endif ()
