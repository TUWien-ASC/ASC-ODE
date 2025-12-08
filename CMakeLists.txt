cmake_minimum_required(VERSION 3.20)
project(ASC_ode LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include_directories(src nanoblas/src)
add_subdirectory(nanoblas)

add_executable(test_ode demos/test_ode.cpp)
target_link_libraries(test_ode PUBLIC nanoblas)

add_executable(demo_autodiff demos/demo_autodiff.cpp)
target_link_libraries(demo_autodiff PUBLIC nanoblas)

# Python bindings
find_package(pybind11 REQUIRED CONFIG)

# 🔥 FIX: Also link pybind11 for the C++ executables
target_link_libraries(test_ode PRIVATE pybind11::pybind11)
target_link_libraries(demo_autodiff PRIVATE pybind11::pybind11)

pybind11_add_module(asc_ode
    src/bindings.cpp
)

target_link_libraries(asc_ode PRIVATE nanoblas)

set_target_properties(asc_ode PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
)
