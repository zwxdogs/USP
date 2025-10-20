cmake_minimum_required(VERSION 3.10)

project(USP)
# set(CMAKE_BUILD_TYPE "Release")
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

find_package(Eigen3 REQUIRED)
find_package(VTK REQUIRED)

file(GLOB SRC_FILES
    "${PROJECT_SOURCE_DIR}/src/*.hpp"
    "${PROJECT_SOURCE_DIR}/src/*.cpp")

add_executable(USP test/test.cpp src/tools/VTK/point_cloud.cpp)

target_link_libraries(USP PRIVATE ${VTK_LIBRARIES})
