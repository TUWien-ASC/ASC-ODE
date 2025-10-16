cmake_minimum_required(VERSION 3.20)
project(ASC_ode)

set (CMAKE_CXX_STANDARD 20)


# use linalg from NGSolve, otherwise comment it out
# find_package(NGSolve CONFIG REQUIRED
#   HINTS /usr $ENV{NETGENDIR}/.. /tmp/inst/usr /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen
# )


include_directories(src nanoblas/src)

add_subdirectory (src)
add_subdirectory (nanoblas)


add_executable (test_ode demos/test_ode.cpp)
# target_link_libraries (test_ode PRIVATE LAPACK::LAPACK)
# target_link_libraries (test_ode PUBLIC ngbla)

# add_executable (test_newmark demos/test_newmark.cpp)
# target_link_libraries (test_newmark PUBLIC ngbla)

# add_executable (test_alpha demos/test_alpha.cpp)
# target_link_libraries (test_alpha PUBLIC ngbla)


# add_subdirectory (mass_spring)
