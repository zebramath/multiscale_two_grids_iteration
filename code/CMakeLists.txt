cmake_minimum_required(VERSION 3.16)
project(two_grids_iteration VERSION 7.1.0 LANGUAGES CXX)

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

tgi_add_program(experiment1_two_grid_comparison
    experiments/experiment1_two_grid_comparison.cpp)
tgi_add_program(experiment2_finite_path
    experiments/experiment2_finite_path.cpp)
tgi_add_program(experiment3_oracle_validation
    experiments/experiment3_oracle_validation.cpp)
tgi_add_program(experiment4_submission_robustness
    experiments/experiment4_submission_robustness.cpp)
tgi_add_program(experiment5_stopping_ablation
    experiments/experiment5_stopping_ablation.cpp)
tgi_add_program(experiment6_fixed_physical_refinement
    experiments/experiment6_fixed_physical_refinement.cpp)
tgi_add_program(experiment7_multilevel_pilot
    experiments/experiment7_multilevel_pilot.cpp)
add_custom_target(experiments
    DEPENDS experiment1_two_grid_comparison experiment2_finite_path
            experiment3_oracle_validation experiment4_submission_robustness
            experiment5_stopping_ablation
            experiment6_fixed_physical_refinement
            experiment7_multilevel_pilot)

include(CTest)
if(BUILD_TESTING)
    add_executable(unit_core tests/unit_core.cpp)
    target_link_libraries(unit_core PRIVATE tgi_headers)
    add_test(NAME unit_core COMMAND unit_core)
    add_executable(regression_core tests/regression_core.cpp)
    target_link_libraries(regression_core PRIVATE tgi_headers)
    add_test(NAME regression_core COMMAND regression_core)
endif()
