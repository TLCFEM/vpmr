cmake_minimum_required(VERSION 3.24)
project(vpmr)

set(CMAKE_CXX_STANDARD 17)

include_directories(eigen)
include_directories(eigen/unsupported)
include_directories(mpreal)
include_directories(exprtk)
include_directories(exprtk-custom-types)

link_libraries(gmp mpfr tbb)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    # add coverage
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -fprofile-arcs -ftest-coverage")
endif ()

add_executable(vpmr src/VPMR.cpp)
