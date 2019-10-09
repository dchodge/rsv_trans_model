cmake_minimum_required (VERSION 3.0)

project (rsv_model)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories(
	headers 
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
	boost_random 
)
