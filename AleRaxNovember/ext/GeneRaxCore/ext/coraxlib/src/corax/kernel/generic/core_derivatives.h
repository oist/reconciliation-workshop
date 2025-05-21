#ifndef CORAX_KERNEL_GENERIC_CORE_DERIVATIVES_H_
#define CORAX_KERNEL_GENERIC_CORE_DERIVATIVES_H_

#include "corax/core/common.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in core_derivatives.c */

  CORAX_EXPORT int
  corax_core_update_sumtable_repeats(unsigned int        states,
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

  CORAX_EXPORT int
  corax_core_update_sumtable_repeats_generic(unsigned int        states,
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
  CORAX_EXPORT int
  corax_core_update_sumtable_ti_4x4(unsigned int         sites,
                                    unsigned int         rate_cats,
                                    const double *       parent_clv,
                                    const unsigned char *left_tipchars,
                                    const unsigned int * parent_scaler,
                                    double *const *      eigenvecs,
                                    double *const *      inv_eigenvecs,
                                    double *const *      freqs,
                                    double *             sumtable,
                                    unsigned int         attrib);

  CORAX_EXPORT int
  corax_core_update_sumtable_ii(unsigned int        states,
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
  corax_core_update_sumtable_ti(unsigned int         states,
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

  CORAX_EXPORT int
  corax_core_likelihood_derivatives(unsigned int        states,
                                    unsigned int        sites,
                                    unsigned int        rate_cats,
                                    const double *      rate_weights,
                                    const unsigned int *parent_scaler,
                                    const unsigned int *child_scaler,
                                    unsigned int        parent_ids,
                                    unsigned int        child_ids,
                                    const int *         invariant,
                                    const unsigned int *pattern_weights,
                                    double              branch_length,
                                    const double *      prop_invar,
                                    double *const *     freqs,
                                    const double *      rates,
                                    double *const *     eigenvals,
                                    const double *      sumtable,
                                    double *            d_f,
                                    double *            dd_f,
                                    unsigned int        attrib);

  CORAX_EXPORT int
  corax_core_update_sumtable_repeats_avx(unsigned int        states,
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


#endif /* CORAX_KERNEL_GENERIC_CORE_DERIVATIVES_H_ */
