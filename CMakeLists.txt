cmake_minimum_required(VERSION 3.10)

project(USP)

# set(CMAKE_BUILD_TYPE "Debug")
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

find_package(Eigen3 REQUIRED)

file(GLOB SRC_FILES
    "${PROJECT_SOURCE_DIR}/src/*.hpp"
    "${PROJECT_SOURCE_DIR}/src/*.cpp")

add_subdirectory(extern/pybind11)

# add_executable(USP test/test.cpp src/tools/VTK/point_cloud.cpp)
pybind11_add_module(straighjt_flow_2d bind.cpp)

