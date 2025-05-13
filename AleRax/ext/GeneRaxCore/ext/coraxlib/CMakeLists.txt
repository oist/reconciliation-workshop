cmake_minimum_required(VERSION 3.11.0 FATAL_ERROR)
project(corax)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
include(simd)
include(warningFlags)

# GoogleTest needs at least C++14; we decided against C++17 for now, as it's not widely available on
# cluster systems.
set(CMAKE_CXX_STANDARD 14)
include(CTest)

# Enable ctest for the current subdirectory and below.
enable_testing()

# Options - user settable
option(CORAX_BUILD_DOCS "Build (Doxygen) documentation" Off)
option(CORAX_BUILD_TESTS "Build the automatic tests" Off)
option(CORAX_BUILD_BENCHMARKS "Build the benchmarks" Off)
option(CORAX_BUILD_DIFFICULTY_PREDICTION "Build difficulty prediction module" Off)
option(CORAX_ENABLE_SSE "Enable SSE SIMD kernels if supported by compiler, OS, and CPU" On)
option(CORAX_ENABLE_AVX "Enable AVX SIMD kernels if supported by compiler, OS, and CPU" On)
option(CORAX_ENABLE_AVX2 "Enable AVX2 SIMD kernels if supported by compiler, OS, and CPU" On)
option(CORAX_BUILD_SHARED "Build a shared library" Off)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
    "${PROJECT_SOURCE_DIR}/bin"
    CACHE PATH "Where put all LIBRARY target files when built."
)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    "${PROJECT_SOURCE_DIR}/bin"
    CACHE PATH "Where to put all ARCHIVE target files when built."
)

if (DEFINED BUILD_SHARED_LIBS AND NOT (CORAX_BUILD_SHARED))
    message(
            WARNING
            "BUILD_SHARED_LIBS is set, but CORAX_BUILD_SHARED is Off. Setting CORAX_BUILD_SHARED to On."
    )
    set(CORAX_BUILD_SHARED On)
endif ()

# Default the build type to Release
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
            "Release"
            CACHE STRING "Build type" FORCE
            )
    message(STATUS "CMAKE_BUILD_TYPE not set, defaulting to Release")
endif ()

# Check for available SIMD variants and disable not available ones.
check_sse_available()
if (NOT SSE_AVAILABLE)
    message(STATUS "SSE SIMD not supported -> disabling")
    set(CORAX_ENABLE_SSE Off)
endif ()

check_avx_available()
if (NOT AVX_AVAILABLE)
    message(STATUS "AVX SIMD not supported -> disabling")
    set(CORAX_ENABLE_AVX Off)
endif ()

check_avx2_available()
if (NOT AVX2_AVAILABLE)
    message(STATUS "AVX2 SIMD not supported -> disabling")
    set(CORAX_ENABLE_AVX2 Off)
endif ()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
if (CORAX_BUILD_SHARED)
    message(STATUS "Building coraxlib as a shared library.")
else ()
    message(STATUS "Building coraxlib as a static library.")
endif ()
message(STATUS "Building documentation: ${CORAX_BUILD_DOCS}")
message(STATUS "Building tests: ${CORAX_BUILD_TESTS}")
message(STATUS "Building benchmarks: ${CORAX_BUILD_BENCHMARKS}")
message(STATUS "Building difficutly prediction: ${CORAX_BUILD_DIFFICULTY_PREDICTION}")
message(STATUS "Enable SSE SIMD kernels: ${CORAX_ENABLE_SSE}")
message(STATUS "Enable AVX SIMD kernels: ${CORAX_ENABLE_AVX}")
message(STATUS "Enable AVX2 SIMD kernels: ${CORAX_ENABLE_AVX2}")

# Configuration which is not user settable
set(CORAX_INCLUDE_DIRS
        ${PROJECT_SOURCE_DIR}/src
        CACHE INTERNAL "${PROJECT_NAME}: Include dirs"
        )

# Build the difficulty prediction if requested.
if (CORAX_BUILD_DIFFICULTY_PREDICTION)
    add_subdirectory(lib/difficulty_prediction)
    list(APPEND CORAX_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/lib/difficulty_prediction/src)
endif ()

# Build the corax library
add_subdirectory(src/corax)

# Build the documentation if requested.
if (CORAX_BUILD_DOCS)
    add_subdirectory(docs)
endif ()

# Build the tests if requested.
if (CORAX_BUILD_TESTS)
    add_subdirectory(test)
    if(CORAX_BUILD_DIFFICULTY_PREDICTION)
        add_subdirectory(lib/difficulty_prediction/test)
    endif()
endif (CORAX_BUILD_TESTS)

# Build the benchmarks if requested.
if (CORAX_BUILD_BENCHMARKS)
    add_subdirectory(bench)
endif (CORAX_BUILD_BENCHMARKS)
