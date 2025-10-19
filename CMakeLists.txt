cmake_minimum_required(VERSION 3.10)

project(USP)

find_package(Eigen3 REQUIRED)

file(GLOB SRC_FILES
    "${PROJECT_SOURCE_DIR}/src/*.H"
    "${PROJECT_SOURCE_DIR}/src/*.C")

add_executable(USP test/test.C)
