#ifndef CORAX_KERNEL_SSE_CORE_DERIVATIVES_H_
#define CORAX_KERNEL_SSE_CORE_DERIVATIVES_H_

#include "corax/core/common.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in core_derivatives_sse.c */

  CORAX_EXPORT int
  corax_core_update_sumtable_ii_sse(unsigned int        states,
                                    unsigned int        sites,
                                    unsigned int        rate_cats,
                                    const double *      parent_clv,
                                    const double *      child_clv,
                                    const unsigned int *parent_scaler,
                                    const unsigned int *child_scaler,
                                    double *const *     eigenvecs,
                                    double *const *     inv_eigenvecs,
                                    double *const *     freqs,
                                    double *            sumtable,
                                    unsigned int        attrib);

  CORAX_EXPORT int
  corax_core_update_sumtable_ti_sse(unsigned int         states,
                                    unsigned int         sites,
                                    unsigned int         rate_cats,
                                    const double *       parent_clv,
                                    const unsigned char *left_tipchars,
                                    const unsigned int * parent_scaler,
                                    double *const *      eigenvecs,
                                    double *const *      inv_eigenvecs,
                                    double *const *      freqs,
                                    const corax_state_t *tipmap,
                                    double *             sumtable,
                                    unsigned int         attrib);

  CORAX_EXPORT int corax_core_update_sumtable_repeats_generic_sse(
      unsigned int        states,
      unsigned int        sites,
      unsigned int        parent_sites,
      unsigned int        rate_cats,
      const double *      clvp,
      const double *      clvc,
      const unsigned int *parent_scaler,
      const unsigned int *child_scaler,
      double *const *     eigenvecs,
      double *const *     inv_eigenvecs,
      double *const *     freqs,
      double *            sumtable,
      const unsigned int *parent_site_id,
      const unsigned int *child_site_id,
      double *            bclv_buffer,
      unsigned int        inv,
      unsigned int        attrib);

#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CORAX_KERNEL_SSE_CORE_DERIVATIVES_H_ */
