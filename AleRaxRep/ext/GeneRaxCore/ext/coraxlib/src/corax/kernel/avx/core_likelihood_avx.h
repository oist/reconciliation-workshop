#ifndef CORAX_KERNEL_AVX_CORE_LIKELIHOOD_H_
#define CORAX_KERNEL_AVX_CORE_LIKELIHOOD_H_

#include "corax/core/common.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in core_likelihood_avx.c */

  CORAX_EXPORT double
  corax_core_edge_loglikelihood_ii_avx(unsigned int         states,
                                       unsigned int         sites,
                                       unsigned int         rate_cats,
                                       const double *       parent_clv,
                                       const unsigned int * parent_scaler,
                                       const double *       child_clv,
                                       const unsigned int * child_scaler,
                                       const double *       pmatrix,
                                       const double *const *frequencies,
                                       const double *       rate_weights,
                                       const unsigned int * pattern_weights,
                                       const double *       invar_proportion,
                                       const int *          invar_indices,
                                       const unsigned int * freqs_indices,
                                       double *             persite_lnl,
                                       unsigned int         attrib);

  CORAX_EXPORT double
  corax_core_edge_loglikelihood_ii_4x4_avx(unsigned int         sites,
                                           unsigned int         rate_cats,
                                           const double *       parent_clv,
                                           const unsigned int * parent_scaler,
                                           const double *       child_clv,
                                           const unsigned int * child_scaler,
                                           const double *       pmatrix,
                                           const double *const *frequencies,
                                           const double *       rate_weights,
                                           const unsigned int * pattern_weights,
                                           const double *      invar_proportion,
                                           const int *         invar_indices,
                                           const unsigned int *freqs_indices,
                                           double *            persite_lnl,
                                           unsigned int        attrib);

  CORAX_EXPORT double
  corax_core_edge_loglikelihood_ti_4x4_avx(unsigned int         sites,
                                           unsigned int         rate_cats,
                                           const double *       parent_clv,
                                           const unsigned int * parent_scaler,
                                           const unsigned char *tipchars,
                                           const double *       pmatrix,
                                           const double *const *frequencies,
                                           const double *       rate_weights,
                                           const unsigned int * pattern_weights,
                                           const double *      invar_proportion,
                                           const int *         invar_indices,
                                           const unsigned int *freqs_indices,
                                           double *            persite_lnl,
                                           unsigned int        attrib);

  CORAX_EXPORT double corax_core_edge_loglikelihood_ti_20x20_avx(
      unsigned int         sites,
      unsigned int         rate_cats,
      const double *       parent_clv,
      const unsigned int * parent_scaler,
      const unsigned char *tipchars,
      const corax_state_t *tipmap,
      unsigned int         tipmap_size,
      const double *       pmatrix,
      const double *const *frequencies,
      const double *       rate_weights,
      const unsigned int * pattern_weights,
      const double *       invar_proportion,
      const int *          invar_indices,
      const unsigned int * freqs_indices,
      double *             persite_lnl,
      unsigned int         attrib);

  CORAX_EXPORT double
  corax_core_edge_loglikelihood_ti_avx(unsigned int         states,
                                       unsigned int         sites,
                                       unsigned int         rate_cats,
                                       const double *       parent_clv,
                                       const unsigned int * parent_scaler,
                                       const unsigned char *tipchars,
                                       const corax_state_t *tipmap,
                                       const double *       pmatrix,
                                       const double *const *frequencies,
                                       const double *       rate_weights,
                                       const unsigned int * pattern_weights,
                                       const double *       invar_proportion,
                                       const int *          invar_indices,
                                       const unsigned int * freqs_indices,
                                       double *             persite_lnl,
                                       unsigned int         attrib);

  CORAX_EXPORT double
  corax_core_root_loglikelihood_4x4_avx(unsigned int         sites,
                                        unsigned int         rate_cats,
                                        const double *       clv,
                                        const unsigned int * scaler,
                                        const double *const *frequencies,
                                        const double *       rate_weights,
                                        const unsigned int * pattern_weights,
                                        const double *       invar_proportion,
                                        const int *          invar_indices,
                                        const unsigned int * freqs_indices,
                                        double *             persite_lnl);

  CORAX_EXPORT double
  corax_core_root_loglikelihood_avx(unsigned int         states,
                                    unsigned int         sites,
                                    unsigned int         rate_cats,
                                    const double *       clv,
                                    const unsigned int * scaler,
                                    const double *const *frequencies,
                                    const double *       rate_weights,
                                    const unsigned int * pattern_weights,
                                    const double *       invar_proportion,
                                    const int *          invar_indices,
                                    const unsigned int * freqs_indices,
                                    double *             persite_lnl);

  CORAX_EXPORT double
  corax_core_root_loglikelihood_repeats_avx(unsigned int         states,
                                            unsigned int         sites,
                                            unsigned int         rate_cats,
                                            const double *       clv,
                                            const unsigned int * site_id,
                                            const unsigned int * scaler,
                                            const double *const *frequencies,
                                            const double *       rate_weights,
                                            const unsigned int *pattern_weights,
                                            const double *invar_proportion,
                                            const int *   invar_indices,
                                            const unsigned int *freqs_indices,
                                            double *            persite_lnl);

  CORAX_EXPORT double corax_core_edge_loglikelihood_repeats_generic_avx(
      unsigned int         states,
      unsigned int         sites,
      const unsigned int   child_sites,
      unsigned int         rate_cats,
      const double *       parent_clv,
      const unsigned int * parent_scaler,
      const double *       child_clv,
      const unsigned int * child_scaler,
      const double *       pmatrix,
      const double *const *frequencies,
      const double *       rate_weights,
      const unsigned int * pattern_weights,
      const double *       invar_proportion,
      const int *          invar_indices,
      const unsigned int * freqs_indices,
      double *             persite_lnl,
      const unsigned int * parent_site_id,
      const unsigned int * child_site_id,
      double *             bclv,
      unsigned int         attrib);

  CORAX_EXPORT double corax_core_edge_loglikelihood_repeats_4x4_avx(
      unsigned int         states,
      unsigned int         sites,
      const unsigned int   child_sites,
      unsigned int         rate_cats,
      const double *       parent_clv,
      const unsigned int * parent_scaler,
      const double *       child_clv,
      const unsigned int * child_scaler,
      const double *       pmatrix,
      const double *const *frequencies,
      const double *       rate_weights,
      const unsigned int * pattern_weights,
      const double *       invar_proportion,
      const int *          invar_indices,
      const unsigned int * freqs_indices,
      double *             persite_lnl,
      const unsigned int * parent_site_id,
      const unsigned int * child_site_id,
      double *             bclv,
      unsigned int         attrib);

  CORAX_EXPORT double corax_core_edge_loglikelihood_repeatsbclv_4x4_avx(
      unsigned int         states,
      unsigned int         sites,
      const unsigned int   child_sites,
      unsigned int         rate_cats,
      const double *       parent_clv,
      const unsigned int * parent_scaler,
      const double *       child_clv,
      const unsigned int * child_scaler,
      const double *       pmatrix,
      const double *const *frequencies,
      const double *       rate_weights,
      const unsigned int * pattern_weights,
      const double *       invar_proportion,
      const int *          invar_indices,
      const unsigned int * freqs_indices,
      double *             persite_lnl,
      const unsigned int * parent_site_id,
      const unsigned int * child_site_id,
      double *             bclv,
      unsigned int         attrib);

#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CORAX_KERNEL_AVX_CORE_LIKELIHOOD_H_ */
