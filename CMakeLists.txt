cmake_minimum_required(VERSION 3.24)
project(BA)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_FLAGS "-pthread")
set(GCC_COVERAGE_COMPILE_FLAGS "")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

add_executable(BA src/cpp/main.cpp)
