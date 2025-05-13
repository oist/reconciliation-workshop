cmake_minimum_required(VERSION 3.0.2 FATAL_ERROR)
project(cmake-example)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../../" "bin")

add_executable(newton)
target_sources(newton PRIVATE newton.c)
target_include_directories(newton PRIVATE "${CORAX_INCLUDE_DIRS}")
target_link_libraries(newton PRIVATE ${CORAX_LIB_TARGET})

