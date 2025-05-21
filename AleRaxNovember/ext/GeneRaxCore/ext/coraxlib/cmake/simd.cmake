# Flags to enable the different SIMD variants.
set(SSE_FLAGS "-msse3")
set(AVX_FLAGS "-mavx")
set(AVX2_FLAGS "-mfma -mavx2")

# Check if the compiler, OS, and architecture support the different SIMD variants. It seems like the
# only reasonable way to check if a SIMD variant is supported by the whole toolchain is to actually
# compile small test programs.
function (check_sse_available)
    set(TEST_CODE " #include <immintrin.h>
        int main() {__m128d a = _mm_setzero_pd();  return 1;}"
    )
    set(TEST_FILE ${CMAKE_CURRENT_BINARY_DIR}/test_sse.c)
    file(WRITE "${TEST_FILE}" "${TEST_CODE}")
    try_compile(
        SSE_COMPILED ${CMAKE_CURRENT_BINARY_DIR}
        ${TEST_FILE}
        COMPILE_DEFINITIONS ${SSE_FLAGS}
    )
    set(SSE_AVAILABLE
        ${SSE_COMPILED}
        PARENT_SCOPE
    )
endfunction ()

function (check_avx_available)
    set(TEST_CODE " #include <immintrin.h>
        int main() {__m256d a = _mm256_setzero_pd ();  return 1;}"
    )
    set(TEST_FILE ${CMAKE_CURRENT_BINARY_DIR}/test_avx.c)
    file(WRITE "${TEST_FILE}" "${TEST_CODE}")
    try_compile(
        AVX_COMPILED ${CMAKE_CURRENT_BINARY_DIR}
        ${TEST_FILE}
        COMPILE_DEFINITIONS ${AVX_FLAGS}
    )
    set(AVX_AVAILABLE
        ${AVX_COMPILED}
        PARENT_SCOPE
    )
endfunction ()

function (check_avx2_available)
    set(TEST_CODE " #include <immintrin.h>
        int main() {__m256i a, b; b =  _mm256_abs_epi16(a); return 1;}"
    )
    set(TEST_FILE ${CMAKE_CURRENT_BINARY_DIR}/test_avx2.c)
    file(WRITE "${TEST_FILE}" "${TEST_CODE}")
    try_compile(
        AVX2_COMPILED ${CMAKE_CURRENT_BINARY_DIR}
        ${TEST_FILE}
        COMPILE_DEFINITIONS ${AVX2_FLAGS}
    )
    set(AVX2_AVAILABLE
        ${AVX2_COMPILED}
        PARENT_SCOPE
    )
endfunction ()

# Add compile definitions for the different SIMD variants.
function(__corax_add_simd_definitions TARGET SIMD_VARIANT)
    target_compile_definitions(${TARGET} PRIVATE "-DHAVE_${SIMD_VARIANT}" "-DHAVE_X86INTRIN_H")
    target_compile_options(${TARGET} PRIVATE -march=native)
endfunction()

function (target_add_sse_definitions TARGET)
    __corax_add_simd_definitions(${TARGET} SSE3)
endfunction ()

function (target_add_avx_definitions TARGET)
    __corax_add_simd_definitions(${TARGET} AVX)
endfunction ()

function (target_add_avx2_definitions TARGET)
    __corax_add_simd_definitions(${TARGET} AVX2)
endfunction ()
