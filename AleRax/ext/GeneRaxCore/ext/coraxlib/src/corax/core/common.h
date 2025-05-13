#ifndef CORAX_CORE_COMMON_H_
#define CORAX_CORE_COMMON_H_

#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#if (!defined(__clang__) && defined(__GNUC__)                                  \
     && (__GNUC__ < 4 || (__GNUC__ == 4 && __GNUC_MINOR__ < 7)))
#if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 6))
#if (defined(HAVE_AVX2))
#error "GCC 4.6.x. Please run ./configure --disable-avx2"
#endif
#else
#if (defined(HAVE_AVX2) || defined(HAVE_AVX))
#error "GCC < 4.6. Please run ./configure --disable-avx --disable-avx2"
#endif
#endif
#endif

#ifdef HAVE_X86INTRIN_H
#include <x86intrin.h>
#endif

/* platform specific */

#if (!defined(__APPLE__) && !defined(__WIN32__) && !defined(__WIN64__))
#include <sys/sysinfo.h>
#endif

#if (defined(__WIN32__) || defined(__WIN64__))
#define CORAX_EXPORT __declspec(dllexport)
#else
#define CORAX_EXPORT
#endif

/* macros */

#define CORAX_MIN(a, b) ((a) < (b) ? (a) : (b))
#define CORAX_MAX(a, b) ((a) > (b) ? (a) : (b))
#define CORAX_SWAP(x, y)                                                       \
  do {                                                                         \
    __typeof__(x) _t = x;                                                      \
    x                = y;                                                      \
    y                = _t;                                                     \
  } while (0)
#define CORAX_HAS_CPU_FEATURE(x)                                               \
  ((corax_hardware.is_initialized || corax_hardware_probe())                   \
   && corax_hardware.x)
#define CORAX_UNUSED(expr)                                                     \
  do {                                                                         \
    (void)(expr);                                                              \
  } while (0)

#define CORAX_SUBST_RATE_COUNT(states) (states * (states - 1) / 2)
#define CORAX_SUBST_RATE_COUNT_NONREV(states) (states * (states - 1))

/** @defgroup corax_defines Constant Definitions
 * @{
 */
/* constants */
#define CORAX_FAILURE 0
#define CORAX_SUCCESS 1

#define CORAX_FALSE 0
#define CORAX_TRUE 1

/* branch linkage modes */
#define CORAX_BRLEN_LINKED 0
#define CORAX_BRLEN_SCALED 1
#define CORAX_BRLEN_UNLINKED 2

/* parallel reduction modes */
#define CORAX_REDUCE_SUM 0
#define CORAX_REDUCE_MAX 1
#define CORAX_REDUCE_MIN 2

/* memory block alignment of respective SIMD instructions */
#define CORAX_ALIGNMENT_CPU 8
#define CORAX_ALIGNMENT_SSE 16
#define CORAX_ALIGNMENT_AVX 32

#define CORAX_LINEALLOC 2048

#define CORAX_ASCII_SIZE 256

#define CORAX_ERRMSG_LEN 200

#define CORAX_SCALE_FACTOR                                                     \
  115792089237316195423570985008687907853269984665640564039457584007913129639936.0 /*  2**256 (exactly)  */
#define CORAX_SCALE_THRESHOLD (1.0 / CORAX_SCALE_FACTOR)
#define CORAX_SCALE_FACTOR_SQRT                                                \
  340282366920938463463374607431768211456.0 /* 2**128 */
#define CORAX_SCALE_THRESHOLD_SQRT (1.0 / CORAX_SCALE_FACTOR_SQRT)
#define CORAX_SCALE_BUFFER_NONE -1

/* in per-rate scaling mode, maximum difference between scalers
 * please see https://github.com/xflouris/libpll/issues/44  */
#define CORAX_SCALE_RATE_MAXDIFF 4

#define CORAX_MISC_EPSILON 1e-8
#define CORAX_ONE_EPSILON 1e-15
#define CORAX_ONE_MIN (1 - CORAX_ONE_EPSILON)
#define CORAX_ONE_MAX (1 + CORAX_ONE_EPSILON)
#define CORAX_EIGEN_MINFREQ 1e-6

#define CORAX_TREE_DEFAULT_BRANCH_LENGTH 0.1
/** @} */

// TODO: this must be adapted for MSVC
#define CORAX_POPCNT32 __builtin_popcount
#define CORAX_POPCNT64 __builtin_popcountll
#define CORAX_CTZ32 __builtin_ctz
#define CORAX_CTZ64 __builtin_ctzll

/* structures and data types */

#define CORAX_STATE_POPCNT CORAX_POPCNT64
#define CORAX_STATE_CTZ CORAX_CTZ64

/** @defgroup corax_errors Error Codes
 * Error codes for coraxlib.
 * @{
 */
/* error codes */
#define CORAX_ERROR_NOT_IMPLEMENTED 13
#define CORAX_ERROR_INVALID_RANGE 21
#define CORAX_ERROR_INVALID_NODE_TYPE 22
#define CORAX_ERROR_INVALID_INDEX 23
#define CORAX_ERROR_INVALID_PARAM 24
#define CORAX_ERROR_INVALID_TREE 25
#define CORAX_ERROR_INVALID_TREE_SIZE 26
#define CORAX_ERROR_INVALID_SPLIT 27
#define CORAX_ERROR_INVALID_THRESHOLD 28

#define CORAX_ERROR_FILE_OPEN 100
#define CORAX_ERROR_FILE_SEEK 101
#define CORAX_ERROR_FILE_EOF 102
#define CORAX_ERROR_FASTA_ILLEGALCHAR 201
#define CORAX_ERROR_FASTA_UNPRINTABLECHAR 202
#define CORAX_ERROR_FASTA_INVALIDHEADER 203
#define CORAX_ERROR_FASTA_NONALIGNED 204
#define CORAX_ERROR_PHYLIP_SYNTAX 231
#define CORAX_ERROR_PHYLIP_LONGSEQ 232
#define CORAX_ERROR_PHYLIP_NONALIGNED 233
#define CORAX_ERROR_PHYLIP_ILLEGALCHAR 234
#define CORAX_ERROR_PHYLIP_UNPRINTABLECHAR 235
#define CORAX_ERROR_NEWICK_SYNTAX 111
#define CORAX_ERROR_MEM_ALLOC 112
#define CORAX_ERROR_TIPDATA_ILLEGALSTATE 114
#define CORAX_ERROR_TIPDATA_ILLEGALFUNCTION 115
#define CORAX_ERROR_TREE_CONVERSION 116
#define CORAX_ERROR_INVAR_INCOMPAT 117
#define CORAX_ERROR_INVAR_PROPORTION 118
#define CORAX_ERROR_INVAR_PARAMINDEX 119
#define CORAX_ERROR_INVAR_NONEFOUND 120
#define CORAX_ERROR_AB_INVALIDMETHOD 121
#define CORAX_ERROR_AB_NOSUPPORT 122
#define CORAX_ERROR_SPR_TERMINALBRANCH 123
#define CORAX_ERROR_SPR_NOCHANGE 124
#define CORAX_ERROR_NNI_INVALIDMOVE 125
#define CORAX_ERROR_NNI_TERMINALBRANCH 126
#define CORAX_ERROR_STEPWISE_STRUCT 127
#define CORAX_ERROR_STEPWISE_TIPS 128
#define CORAX_ERROR_STEPWISE_UNSUPPORTED 129
#define CORAX_ERROR_EINVAL 130
#define CORAX_ERROR_MSA_EMPTY 131
#define CORAX_ERROR_MSA_MAP_INVALID 132
/** @} */


typedef unsigned long long corax_state_t;
typedef int                corax_bool_t;

CORAX_EXPORT extern __thread int              corax_errno;
CORAX_EXPORT extern __thread char             corax_errmsg[CORAX_ERRMSG_LEN];


#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in common.c */

  void corax_set_error(int _errno, const char *errmsg_fmt, ...);
  void corax_reset_error();

  CORAX_EXPORT void *corax_aligned_alloc(size_t size, size_t alignment);
  CORAX_EXPORT void corax_aligned_free(void *ptr);

#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CORAX_CORE_COMMON_H_ */
