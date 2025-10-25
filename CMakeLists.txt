cmake_minimum_required(VERSION 3.18...3.27)

project(USP)

# set(CMAKE_BUILD_TYPE "Debug")
set(DEV_MODULE Development.Module) 
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE) 
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo") 
endif()

set(Python_EXECUTABLE /home/zwxdogs/work/USP/.venv/bin/python) # absolute path.
find_package(Python 3.10 COMPONENTS Interpreter ${DEV_MODULE} REQUIRED)
execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT
)
find_package(nanobind CONFIG REQUIRED)
find_package(Eigen3 REQUIRED)

file(GLOB SRC_FILES
    "${PROJECT_SOURCE_DIR}/src/*.h"
    "${PROJECT_SOURCE_DIR}/src/*.cpp")

nanobind_add_module(linear_mesh src/tools/linear_mesh.cpp)
nanobind_add_module(SIR_calc src/SIR/SIR_calc.cpp)

