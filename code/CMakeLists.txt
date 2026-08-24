cmake_minimum_required(VERSION 3.16)
project(two_grids_iteration VERSION 3.7.0 LANGUAGES CXX)

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

tgi_add_program(experiment1_localization
    experiments/experiment1_localization.cpp
    src/experiment/studies/localization_radius.cpp
    src/experiment/studies/local_tolerance.cpp
    src/experiment/studies/global_pruning.cpp
    src/experiment/studies/support_strategy.cpp)
tgi_add_program(experiment2_pcg_path
    experiments/experiment2_pcg_path.cpp
    src/experiment/studies/pcg_budget.cpp
    src/experiment/studies/pcg_dense_scan.cpp
    src/experiment/studies/pilot_diagnostic.cpp)
tgi_add_program(experiment3_adaptive_oracle
    experiments/experiment3_adaptive_oracle.cpp
    src/experiment/studies/adaptive_trace.cpp
    src/experiment/studies/oracle_quality.cpp
    src/experiment/studies/selector_ablation.cpp)
tgi_add_program(experiment4_robustness
    experiments/experiment4_robustness.cpp
    src/experiment/studies/robustness_matrix.cpp)
tgi_add_program(experiment5_diagnostics
    experiments/experiment5_diagnostics.cpp
    src/experiment/studies/diagnostic_cases.cpp
    src/experiment/studies/stress_cases.cpp)
tgi_add_program(experiment6_workload
    experiments/experiment6_workload.cpp
    src/experiment/studies/workload.cpp)

add_custom_target(experiments
    DEPENDS experiment1_localization experiment2_pcg_path
            experiment3_adaptive_oracle experiment4_robustness
            experiment5_diagnostics experiment6_workload)

include(CTest)
if(BUILD_TESTING)
    add_executable(unit_core tests/unit_core.cpp)
    target_link_libraries(unit_core PRIVATE tgi_headers)
    if(TGI_ENABLE_IPO AND TGI_IPO_SUPPORTED)
        set_property(TARGET unit_core PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
    endif()
    add_test(NAME unit_core COMMAND unit_core)
    add_executable(regression_v37 tests/regression_v3.cpp)
    target_link_libraries(regression_v37 PRIVATE tgi_headers)
    add_test(NAME regression_v37 COMMAND regression_v37)
endif()
