cmake_minimum_required(VERSION 3.16)
project(two_grids_iteration VERSION 9.2.0 LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Threads REQUIRED)
find_package(OpenMP QUIET)

option(TGI_WARNINGS_AS_ERRORS "Treat compiler warnings as errors" ON)

add_library(tgi_headers INTERFACE)
target_include_directories(tgi_headers INTERFACE src)
target_compile_features(tgi_headers INTERFACE cxx_std_17)
target_link_libraries(tgi_headers INTERFACE Threads::Threads)
if(OpenMP_CXX_FOUND)
    target_link_libraries(tgi_headers INTERFACE OpenMP::OpenMP_CXX)
endif()
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(
        tgi_headers INTERFACE
        -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion
        -Wduplicated-cond -Wduplicated-branches -Wlogical-op
        -Wnull-dereference -Wformat=2)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    target_compile_options(
        tgi_headers INTERFACE
        -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wsign-conversion
        -Wnull-dereference -Wformat=2)
endif()
if(TGI_WARNINGS_AS_ERRORS AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(tgi_headers INTERFACE -Werror)
endif()

function(tgi_add_program target)
    add_executable(${target} ${ARGN})
    target_compile_definitions(
        ${target} PRIVATE TGI_RESULTS_DIR="${CMAKE_SOURCE_DIR}/results")
    target_link_libraries(${target} PRIVATE tgi_headers)
endfunction()

tgi_add_program(experiment1_finite_path_comparison
    experiments/experiment1_finite_path_comparison.cpp)
tgi_add_program(experiment2_spectral_path
    experiments/experiment2_spectral_path.cpp)
tgi_add_program(experiment3_scaling_propagation
    experiments/experiment3_scaling_propagation.cpp)
tgi_add_program(experiment4_local_diagnostic_export
    experiments/experiment4_local_diagnostic_export.cpp)
tgi_add_program(experiment5_robustness
    experiments/experiment5_robustness.cpp)
tgi_add_program(experiment6_endpoint_comparison
    experiments/experiment6_endpoint_comparison.cpp)
tgi_add_program(experiment7_multilevel_pilot
    experiments/experiment7_multilevel_pilot.cpp)
add_custom_target(experiments
    DEPENDS experiment1_finite_path_comparison experiment2_spectral_path
            experiment3_scaling_propagation
            experiment4_local_diagnostic_export experiment5_robustness
            experiment6_endpoint_comparison experiment7_multilevel_pilot)
