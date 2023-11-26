cmake_minimum_required(VERSION 3.24)
project(vpmr)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -DNDEBUG")

include_directories(eigen)
include_directories(eigen/unsupported)
include_directories(mpreal)
include_directories(exprtk)
include_directories(exprtk-custom-types)
include_directories(pybind11/include)

if (MINGW)
    find_package(TBB REQUIRED)
    link_libraries(TBB::tbb)
else ()
    link_libraries(tbb)
endif ()

link_libraries(gmp mpfr)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -fprofile-arcs -ftest-coverage")
elseif (UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto=auto")
endif ()

include(CheckIncludeFile)
check_include_file("Windows.h" HAVE_WINDOWS_H)
if (HAVE_WINDOWS_H)
    add_definitions(-DHAVE_WINDOWS_H)
endif ()

add_executable(vpmr src/VPMR.cpp)

find_package(Python3 REQUIRED COMPONENTS Development Interpreter)
message(STATUS "Python3: ${Python3_VERSION}")
include_directories(${Python3_INCLUDE_DIRS})
add_subdirectory(pybind11)

pybind11_add_module(_pyvpmr src/VPMR.cpp)
target_compile_definitions(_pyvpmr PRIVATE PYVPMR)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_sources(${PROJECT_NAME} PRIVATE resource/vpmr.rc)
    if (MINGW)
        target_compile_options(${PROJECT_NAME} PRIVATE "-Wa,-mbig-obj")
        target_compile_options(_pyvpmr PRIVATE "-Wa,-mbig-obj")
    endif ()
endif ()
