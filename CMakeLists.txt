cmake_minimum_required(VERSION 3.10)

set(CMAKE_VERBOSE_MAKEFILE ON)

# set the project name
project(ValidateMe VERSION 0.1)

# set output directory of builds
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

option(BUILD_BENCHMARKS "Build Bechmarks" OFF)

# specify global compiler flags
add_compile_options(-std=c++2a -march=native -fopenmp -pedantic -Wall -Wextra -lstdc++fs -static -g -O3)
link_directories($ENV{CONDA_PREFIX}/lib)

# add the tests executable
add_executable(validateme src/validateme.cpp)

# configure benchmarks
#if (BUILD_BENCHMARKS)
    #add_executable(benchmarks src/benchmarks/benchmark_main.cc)
    #target_compile_options(benchmarks PUBLIC -O3 -DNDEBUG -mtune=native)
    #find_package(benchmark)
    #target_link_libraries(benchmarks benchmark::benchmark)
    #target_include_directories(benchmarks PRIVATE ./ext/libsimdpp ./include ./ext/libsimdpp ./ext/bit)
    #install(TARGETS benchmarks DESTINATION .)
#endif()

# specify the C++ standard
set_property(GLOBAL PROPERTY CXX_STANDARD 17)
set_property(GLOBAL PROPERTY CXX_STANDARD_REQUIRED true)

# specify global include directories
include_directories($ENV{CONDA_PREFIX}/include include ext/CLI11/include ext/cpptqdm)

# specify test-specific flags
target_link_libraries(validateme -fopenmp -lstdc++fs -lspdlog)

# specify test-specific libraries
#target_include_directories(validateme PRIVATE ./ext/catch)

install(TARGETS validateme DESTINATION .)
