find_package(benchmark REQUIRED)

set(CORAX_BENCH_DATA_DIR ${CMAKE_CURRENT_SOURCE_DIR}/data)

# Register the corax-bench executable.
add_executable(corax-bench)
target_compile_definitions(corax-bench PRIVATE BENCH_DATA=${CORAX_BENCH_DATA_DIR})
target_include_directories(corax-bench PRIVATE ${GBENCH_INCLUDE_DIRS} ${CORAX_INCLUDE_DIRS})
target_link_libraries(corax-bench PRIVATE ${CORAX_LIB_TARGET} benchmark::benchmark)

add_subdirectory(src)
