cmake_minimum_required(VERSION 3.24)
project(vpmr)

set(CMAKE_CXX_STANDARD 17)

include_directories(eigen)
include_directories(eigen/unsupported)
include_directories(mpreal)
include_directories(exprtk)
include_directories(exprtk-custom-types)

if (MINGW)
    find_package(TBB REQUIRED)
    link_libraries(TBB::tbb)
else ()
    link_libraries(tbb)
endif ()

link_libraries(gmp mpfr)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -fprofile-arcs -ftest-coverage")
endif ()

include(CheckIncludeFile)
CHECK_INCLUDE_FILE("Windows.h" HAVE_WINDOWS_H)
IF (HAVE_WINDOWS_H)
    add_definitions(-DHAVE_WINDOWS_H)
ENDIF ()

add_executable(vpmr src/VPMR.cpp)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND MINGW)
    target_compile_options(${PROJECT_NAME} PRIVATE "-Wa,-mbig-obj")
endif ()
