#ifndef CORAX_KERNEL_AVX_CORE_PMATRIX_H_
#define CORAX_KERNEL_AVX_CORE_PMATRIX_H_

#include "corax/core/common.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in core_pmatrix_avx.c */

  CORAX_EXPORT int
  corax_core_update_pmatrix_4x4_avx(double **           pmatrix,
                                    unsigned int        rate_cats,
                                    const double *      rates,
                                    const double *      branch_lengths,
                                    const unsigned int *matrix_indices,
                                    const unsigned int *params_indices,
                                    const double *      prop_invar,
                                    double *const *     eigenvals,
                                    double *const *     eigenvecs,
                                    double *const *     inv_eigenvecs,
                                    unsigned int        count);

  CORAX_EXPORT int
  corax_core_update_pmatrix_20x20_avx(double **           pmatrix,
                                      unsigned int        rate_cats,
                                      const double *      rates,
                                      const double *      branch_lengths,
                                      const unsigned int *matrix_indices,
                                      const unsigned int *params_indices,
                                      const double *      prop_invar,
                                      double *const *     eigenvals,
                                      double *const *     eigenvecs,
                                      double *const *     inv_eigenvecs,
                                      unsigned int        count);

#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CORAX_KERNEL_AVX_CORE_PMATRIX_H_ */
