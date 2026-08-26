cmake_minimum_required(VERSION 3.16)
project(two_grids_iteration VERSION 4.4.0 LANGUAGES CXX)

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

function(tgi_add_program target)
    add_executable(${target} ${ARGN})
    target_compile_definitions(
        ${target} PRIVATE TGI_RESULTS_DIR="${CMAKE_SOURCE_DIR}/results")
    target_link_libraries(${target} PRIVATE tgi_headers)
    if(TGI_ENABLE_IPO AND TGI_IPO_SUPPORTED)
        set_property(TARGET ${target} PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
    endif()
endfunction()

tgi_add_program(experiment1_two_grid_comparison
    experiments/experiment1_two_grid_comparison.cpp
    src/experiment/studies/two_grid_comparison.cpp)
tgi_add_program(experiment2_finite_path
    experiments/experiment2_finite_path.cpp
    src/experiment/studies/finite_path_evidence.cpp)
tgi_add_program(experiment3_oracle_validation
    experiments/experiment3_oracle_validation.cpp
    src/experiment/studies/oracle_quality.cpp)
tgi_add_program(experiment4_submission_robustness
    experiments/experiment4_submission_robustness.cpp
    src/experiment/studies/submission_robustness.cpp)
add_custom_target(experiments
    DEPENDS experiment1_two_grid_comparison experiment2_finite_path
            experiment3_oracle_validation experiment4_submission_robustness)

include(CTest)
if(BUILD_TESTING)
    add_executable(unit_core tests/unit_core.cpp)
    target_link_libraries(unit_core PRIVATE tgi_headers)
    if(TGI_ENABLE_IPO AND TGI_IPO_SUPPORTED)
        set_property(TARGET unit_core PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
    endif()
    add_test(NAME unit_core COMMAND unit_core)
    add_executable(regression_v44 tests/regression_v44.cpp)
    target_link_libraries(regression_v44 PRIVATE tgi_headers)
    add_test(NAME regression_v44 COMMAND regression_v44)
endif()
