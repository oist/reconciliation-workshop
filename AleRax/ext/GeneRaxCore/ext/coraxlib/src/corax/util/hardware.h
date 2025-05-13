/**
 * @file util/hardware.h
 *
 * Contains functions used to detect CPU features, for example available SIMD
 * variants or a population count instruction.
 */
#ifndef CORAX_UTIL_HARDWARE_H_
#define CORAX_UTIL_HARDWARE_H_

#include "corax/core/common.h"

/**
 * @brief Describes the presence of CPU features (e.g. SIMD variants).
 */
typedef struct corax_hardware_s
{
  int is_initialized;
  /* cpu features */
  int altivec_present;
  int mmx_present;
  int sse_present;
  int sse2_present;
  int sse3_present;
  int ssse3_present;
  int sse41_present;
  int sse42_present;
  int popcnt_present;
  int avx_present;
  int avx2_present;

  /* TODO: add chip,core,mem info */
} corax_hardware_t;

/**
 * @brief Global variable holding information about the CPU (e.g. present SIMD
 * variants).
 */
CORAX_EXPORT extern __thread corax_hardware_t corax_hardware;

#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in hardware.c */

  /**
   * @brief Detect the available CPU features and save them to the global \c
   * corax_hardware variable.
   *
   * @return Always returns \c CORAX_SUCCESS.
   */
  CORAX_EXPORT int corax_hardware_probe(void);

  /**
   * @brief Detects and shows the available CPU features (e.g. SIMD variants).
   *
   * If \c corax_hardware is not initialized, calls corax_hardware_probe()
   * to detect the CPU features, therefore populates the global \c
   * corax_hardware.
   */
  CORAX_EXPORT void corax_hardware_dump(void);

  /**
   * @brief Ignores which hardware is actually available.
   *
   * Unconditionally sets all CPU features stored in the global \c corax_hardware to 1
   * ("available"), regardless of the CPU.
   */
  CORAX_EXPORT void corax_hardware_ignore(void);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CORAX_UTIL_HARDWARE_H_ */
