cmake_minimum_required(VERSION 3.16)
project(two_grids_iteration VERSION 2.1.0 LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Threads REQUIRED)
find_package(OpenMP QUIET)

option(TGI_NATIVE_OPTIMIZATION
       "Compile for the host CPU (faster, less portable)" OFF)

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

function(tgi_add_program target source)
    add_executable(${target} ${source})
    target_compile_definitions(
        ${target} PRIVATE TGI_RESULTS_DIR="${CMAKE_SOURCE_DIR}/results")
    target_link_libraries(${target} PRIVATE tgi_headers)
endfunction()

set(TGI_PROGRAMS
    generate_test_fields
    experiment1
    experiment2
    experiment3
    experiment4)
foreach(program IN LISTS TGI_PROGRAMS)
    tgi_add_program(${program} test/${program}.cpp)
endforeach()
