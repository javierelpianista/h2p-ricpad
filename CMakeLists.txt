cmake_minimum_required(VERSION 3.9)

# Project name and description
project(
    h2p-ricpad VERSION 0.1
    DESCRIPTION "This program allows to solve the Schrödinger equation for the H2+ molecule in the Born-Oppenheimer approximation using the Riccati-Padé method."
    LANGUAGES CXX
    )

# Release build type by default
if ( NOT CMAKE_BUILD_TYPE )
    set(CMAKE_BUILD_TYPE "Release")
endif()

# Compile flags
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")

# C++ standard options
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# GMP, MPFR, MPC
find_package(GMP REQUIRED)
find_package(MPFR REQUIRED)
find_package(MPC REQUIRED)

# Eigen
find_package(Eigen3 REQUIRED)
if ( NOT EIGEN3_FOUND ) 
    message("Eigen3 not found. Maybe set EIGEN3_INCLUDE_DIR?")
endif()

# Boost
find_package(Boost REQUIRED)
message(STATUS "Boost version: ${Boost_VERSION}")

set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

add_subdirectory(ricpad)

# The H2+ program
add_executable(h2p src/h2p.cpp)
target_link_libraries(
    h2p PUBLIC ricpad Boost::boost Eigen3::Eigen mpfr mpc cln ginac )

add_executable(main src/main.cpp)
target_link_libraries(
    main PUBLIC ricpad Boost::boost Eigen3::Eigen mpfr mpc cln ginac )

# Installation
install(TARGETS h2p
    RUNTIME DESTINATION bin)
