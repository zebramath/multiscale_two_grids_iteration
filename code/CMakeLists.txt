cmake_minimum_required(VERSION 3.16)
project(two_grids_iteration VERSION 2.8.0 LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Threads REQUIRED)
find_package(OpenMP QUIET)
include(CheckIPOSupported)

option(TGI_NATIVE_OPTIMIZATION
       "Compile for the host CPU (faster, less portable)" OFF)
option(TGI_ENABLE_IPO
       "Enable interprocedural optimization in Release builds" ON)

add_library(tgi_headers INTERFACE)
target_include_directories(tgi_headers INTERFACE src)
target_compile_features(tgi_headers INTERFACE cxx_std_17)
target_link_libraries(tgi_headers INTERFACE Threads::Threads)
if(OpenMP_CXX_FOUND)
    target_link_libraries(tgi_headers INTERFACE OpenMP::OpenMP_CXX)
endif()
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(
        tgi_headers INTERFACE -Wall -Wextra -Wpedantic)
    if(TGI_NATIVE_OPTIMIZATION)
        target_compile_options(tgi_headers INTERFACE -march=native)
    endif()
endif()

if(TGI_ENABLE_IPO)
    check_ipo_supported(RESULT TGI_IPO_SUPPORTED OUTPUT TGI_IPO_ERROR)
endif()

function(tgi_add_program target source)
    add_executable(${target} ${source})
    target_compile_definitions(
        ${target} PRIVATE TGI_RESULTS_DIR="${CMAKE_SOURCE_DIR}/results")
    target_link_libraries(${target} PRIVATE tgi_headers)
    if(TGI_ENABLE_IPO AND TGI_IPO_SUPPORTED)
        set_property(TARGET ${target} PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
    endif()
endfunction()

tgi_add_program(experiment1 experiments/experiment1.cpp)
tgi_add_program(experiment2 experiments/experiment2.cpp)
tgi_add_program(experiment3 experiments/experiment3.cpp)
tgi_add_program(experiment4 experiments/experiment4.cpp)
tgi_add_program(experiment5 experiments/experiment5.cpp)

add_custom_target(experiments
    DEPENDS experiment1 experiment2 experiment3 experiment4 experiment5)

include(CTest)
if(BUILD_TESTING)
    add_executable(unit_core tests/unit_core.cpp)
    target_link_libraries(unit_core PRIVATE tgi_headers)
    if(TGI_ENABLE_IPO AND TGI_IPO_SUPPORTED)
        set_property(TARGET unit_core PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
    endif()
    add_test(NAME unit_core COMMAND unit_core)
endif()
