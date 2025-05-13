cmake_minimum_required(VERSION 2.7)
project(BIDON_MPI)

find_package(MPI REQUIRED)

set (CMAKE_CXX_STANDARD 11)
include_directories(${MPI_INCLUDE_PATH})

add_executable(fake_program fake_program.cpp)
target_link_libraries(fake_program ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(fake_program PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(fake_program PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
else()
    message(WARNING "Could not find OpenMP.")
endif()
