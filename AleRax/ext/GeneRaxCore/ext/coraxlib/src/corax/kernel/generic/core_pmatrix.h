#ifndef CORAX_KERNEL_GENERIC_CORE_PMATRIX_H_
#define CORAX_KERNEL_GENERIC_CORE_PMATRIX_H_

#include "corax/core/common.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in core_pmatrix.c */

  CORAX_EXPORT int corax_core_update_pmatrix(double            **pmatrix,
                                             unsigned int        states,
                                             unsigned int        rate_cats,
                                             const double       *rates,
                                             const double       *branch_lengths,
                                             const unsigned int *matrix_indices,
                                             const unsigned int *params_indices,
                                             const double       *prop_invar,
                                             double *const      *eigenvals,
                                             double *const      *eigenvecs,
                                             double *const      *inv_eigenvecs,
                                             unsigned int        count,
                                             unsigned int        attrib);
#ifdef CORAX_NONREV
  CORAX_EXPORT int
  corax_core_update_pmatrix_nonrev(double            **pmatrix,
                                   size_t              states,
                                   size_t              rate_cats,
                                   const double       *rates,
                                   const double       *branch_lengths,
                                   const unsigned int *matrix_indices,
                                   const unsigned int *params_indices,
                                   const double       *prop_invar,
                                   double *const      *params,
                                   unsigned int        count,
                                   unsigned int        attrib);
#endif // CORAX_NONREV

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CORAX_KERNEL_GENERIC_CORE_PMATRIX_H_ */
