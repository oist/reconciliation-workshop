include(CheckIncludeFile)
include(targetSourcesLocal)

# Add the corax library
set(CORAX_LIB_TARGET
    "corax"
    CACHE INTERNAL "${PROJECT_NAME}: Library target to link against."
)

if (CORAX_BUILD_SHARED)
    add_library(${CORAX_LIB_TARGET} SHARED)
else ()
    add_library(${CORAX_LIB_TARGET} STATIC)
endif ()

# For other programs to include coraxlib, it has to have position independent code.
if (CORAX_BUILD_SHARED)
    set_property(TARGET ${CORAX_LIB_TARGET} PROPERTY POSITION_INDEPENDENT_CODE 1)
endif ()

# Add sources and recurse into subdirs
target_sources_local(
    ${CORAX_LIB_TARGET}
    PUBLIC corax_core.h
           corax_io.h
           corax_kernel.h
           corax_model.h
           corax_optimize.h
           corax_tree.h
           corax_util.h
           corax.hpp
           corax_topology.h
)

add_subdirectory(core)
add_subdirectory(optimize)
add_subdirectory(tree)
add_subdirectory(io)
add_subdirectory(kernel)
add_subdirectory(model)
add_subdirectory(util)
add_subdirectory(topology)

# Add include directories.
target_include_directories(${CORAX_LIB_TARGET} PUBLIC ${CORAX_INCLUDE_DIRS})

# Add warning flags
target_compile_options(${CORAX_LIB_TARGET} PRIVATE ${CORAX_WARNING_FLAGS})

option(CORAX_NONREV "Enable non-reversible model computation" OFF)

if(CORAX_NONREV)
  target_compile_definitions(${CORAX_LIB_TARGET} PUBLIC -DCORAX_NONREV)
  message("Non-rev enabled")

  find_package(PkgConfig REQUIRED)
  pkg_search_module(BLAS blas)
  pkg_search_module(CBLAS cblas)
  pkg_search_module(LAPACK lapack)
  pkg_search_module(LAPACKE lapacke)
  if((${BLAS_FOUND} AND ${CBLAS_FOUND}) AND (${LAPACK_FOUND} AND ${LAPACKE_FOUND}))
      target_include_directories(${CORAX_LIB_TARGET} PUBLIC 
        ${BLAS_INCLUDE_DIRS}
        ${CBLAS_INCLUDE_DIRS}
        ${LAPACK_INCLUDE_DIRS}
        ${LAPACKE_INCLUDE_DIRS}
        )
      target_link_libraries(${CORAX_LIB_TARGET} PUBLIC 
        ${BLAS_LDFLAGS}
        ${CBLAS_LDFLAGS}
        ${LAPACK_LDFLAGS}
        ${LAPACKE_LDFLAGS}
        )
    else()
        pkg_search_module(OpenBLAS REQUIRED openblas)
        target_include_directories(${CORAX_LIB_TARGET} PUBLIC ${OpenBLAS_INCLUDE_DIRS})
        target_link_libraries(${CORAX_LIB_TARGET} PUBLIC ${OpenBLAS_LDFLAGS})
        if(${LAPACKE_FOUND})
          target_link_libraries(${CORAX_LIB_TARGET} PUBLIC ${LAPACKE_LDFLAGS})
        endif()
    endif()
endif()

# TODO Decide if we want to use DEBUG, NDEBUG, or DEBUG_MODE and unify all usages.
if (CORAX_DEBUG)
    target_compile_definitions(${CORAX_LIB_TARGET} PRIVATE "DEBUG")
endif ()

# We want GNU instead of POSIX behaviour.
target_compile_definitions(${CORAX_LIB_TARGET} PRIVATE "_GNU_SOURCE")

# If the difficulty prediction is enabled, link it to the corax_obj.
if (CORAX_BUILD_DIFFICULTY_PREDICTION)
    target_link_libraries(${CORAX_LIB_TARGET} PRIVATE coraxlib_difficulty_prediction_lib)
endif ()

# Set the SIMD flags corresponding to the requested and supported variant.
if (CORAX_ENABLE_SSE)
    target_add_sse_definitions(${CORAX_LIB_TARGET})
endif ()

if (CORAX_ENABLE_AVX)
    target_add_avx_definitions(${CORAX_LIB_TARGET})
endif ()

if (CORAX_ENABLE_AVX2)
    target_add_avx2_definitions(${CORAX_LIB_TARGET})
endif ()

# Set the library name to "corax", regardless of if we build the static or shared version. We
# currently do not need this, as we call both of them "corax".
# set_target_properties(${CORAX_LIB_TARGET} PROPERTIES OUTPUT_NAME "corax")
