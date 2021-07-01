cmake_minimum_required(VERSION 3.13)

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
set(CMAKE_CXX_FLAGS_RELEASE "-O2")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")

# C++ standard options
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# GMP, MPFR, MPC
find_library(gmp "gmp")
foreach(lib gmp mpfr ginac cln)
    find_library("${lib}" "${lib}")
    if ( ${${lib}} STREQUAL "${lib}-NOTFOUND" ) 
        message( FATAL_ERROR 
"Library ${lib} not found. One way to specify the path to libraries is by adding -DCMAKE_LIBRARY_PATH= <list-of-library-paths> to the cmake call, for example: 
cmake .. -DCMAKE_LIBRARY_PATH=\"/home/user/libs\"
If several paths are specified, they should be separated by ; .")
    endif()
endforeach()

# Eigen
find_package(Eigen3 REQUIRED)
if ( NOT EIGEN3_FOUND ) 
    message("Eigen3 not found. Maybe set EIGEN3_INCLUDE_DIR?")
endif()

# Boost
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

find_package(Boost 1.69 
    COMPONENTS program_options
    REQUIRED
    )

if ( Boost_FOUND ) 
    include_directories(${Boost_INCLUDE_DIRS})
endif()

add_subdirectory(ricpad)

add_executable(h2p-ricpad src/main.cpp)
target_link_libraries(
    h2p-ricpad PUBLIC 
    ricpad 
    Boost::boost
    Boost::program_options
    Eigen3::Eigen
    mpfr
    cln
    ginac 
    )

target_include_directories(
    h2p-ricpad PUBLIC src/include)
