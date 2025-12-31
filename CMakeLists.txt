cmake_minimum_required(VERSION 3.14)
project(orbit_propagator 
    VERSION 1.0.0
    DESCRIPTION "C++ Orbit Propagator with J2 Perturbations"
    LANGUAGES CXX)

# C++17 standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
elseif(MSVC)
    add_compile_options(/W4)
endif()

# Header-only library interface
add_library(orbit_propagator INTERFACE)
target_include_directories(orbit_propagator INTERFACE 
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Examples
option(BUILD_EXAMPLES "Build example programs" ON)
if(BUILD_EXAMPLES)
    add_executable(iss_propagation examples/iss_propagation.cpp)
    target_link_libraries(iss_propagation PRIVATE orbit_propagator)

    add_executable(sunsync_orbit examples/sunsync_orbit.cpp)
    target_link_libraries(sunsync_orbit PRIVATE orbit_propagator)
endif()

# Installation
include(GNUInstallDirs)
install(TARGETS orbit_propagator
    EXPORT orbit_propagator_targets
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(FILES include/orbit_propagator.hpp 
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Package config
install(EXPORT orbit_propagator_targets
    FILE orbit_propagator-config.cmake
    NAMESPACE orbitprop::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/orbit_propagator
)
