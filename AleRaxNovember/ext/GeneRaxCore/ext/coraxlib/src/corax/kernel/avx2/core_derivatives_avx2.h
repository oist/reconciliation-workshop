#ifndef CORAX_KERNEL_AVX2_CORE_DERIVATIVES_H_
#define CORAX_KERNEL_AVX2_CORE_DERIVATIVES_H_

#include "corax/core/common.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in core_derivatives_avx2.c */

  CORAX_EXPORT int
  corax_core_update_sumtable_ii_avx2(unsigned int        states,
                                     unsigned int        sites,
                                     unsigned int        rate_cats,
                                     const double *      clvp,
                                     const double *      clvc,
                                     const unsigned int *parent_scaler,
                                     const unsigned int *child_scaler,
                                     double *const *     eigenvecs,
                                     double *const *     inv_eigenvecs,
                                     double *const *     freqs,
                                     double *            sumtable,
                                     unsigned int        attrib);

  CORAX_EXPORT int
  corax_core_update_sumtable_ti_avx2(unsigned int         states,
                                     unsigned int         sites,
                                     unsigned int         rate_cats,
                                     const double *       parent_clv,
                                     const unsigned char *left_tipchars,
                                     const unsigned int * parent_scaler,
                                     double *const *      eigenvecs,
                                     double *const *      inv_eigenvecs,
                                     double *const *      freqs,
                                     const corax_state_t *tipmap,
                                     unsigned int         tipmap_size,
                                     double *             sumtable,
                                     unsigned int         attrib);

  CORAX_EXPORT
  int corax_core_likelihood_derivatives_avx2(
      unsigned int        states,
      unsigned int        states_padded,
      unsigned int        rate_cats,
      unsigned int        ef_sites,
      const unsigned int *pattern_weights,
      const double *      rate_weights,
      const int *         invariant,
      const double *      prop_invar,
      double *const *     freqs,
      const double *      sumtable,
      const double *      diagptable,
      double *            d_f,
      double *            dd_f);

  CORAX_EXPORT int corax_core_update_sumtable_repeats_generic_avx2(
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


#endif /* CORAX_KERNEL_AVX2_CORE_DERIVATIVES_H_ */
