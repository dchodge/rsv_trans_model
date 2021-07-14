cmake_minimum_required (VERSION 2.8)

project (rsv_model)
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(GCC_COVERAGE_COMPILE_FLAGS "-fsyntax-only")

set(BOOST_ROOT "/home/lshdh2/rsv_trans_model/lib/boost_1_64_0")
set(BOOST_INCLUDE_DIRS "${BOOST_ROOT}/include")
set(BOOST_LIBRARY_DIRS "${BOOST_ROOT}/lib")
set(Boost_NO_SYSTEM_PATHS on CACHE BOOL "Do not search system for Boost" )
set(APPEND CMAKE_PREFIX_PATH ${BOOST_ROOT})
message("BOOST_ROOT is ${BOOST_ROOT}")

include_directories(
        headers 
        ${BOOST_LIBRARY_DIRS}
        ${BOOST_INCLUDE_DIRS} 
)

add_executable(
        rsv_model
        src/pre.h
        src/epmgp.h
        src/cea.h
        src/interventions.h
        src/mcmc.h
        src/model.h             
        src/outcomes.h
        src/writing.h           
        src/main.cpp
)

target_link_libraries(
        rsv_model
)