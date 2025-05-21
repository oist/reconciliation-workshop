#ifndef CORAX_KERNEL_AVX2_CORE_CLV_H_
#define CORAX_KERNEL_AVX2_CORE_CLV_H_

#include "corax/core/common.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in core_clvs_avx2.c */

  CORAX_EXPORT void
  corax_core_update_clv_ti_avx2(unsigned int         states,
                                unsigned int         sites,
                                unsigned int         rate_cats,
                                double *             parent_clv,
                                unsigned int *       parent_scaler,
                                const unsigned char *left_tipchars,
                                const double *       right_clv,
                                const double *       left_matrix,
                                const double *       right_matrix,
                                const unsigned int * right_scaler,
                                const corax_state_t *tipmap,
                                unsigned int         tipmap_size,
                                unsigned int         attrib);

  CORAX_EXPORT
  void corax_core_update_clv_ti_20x20_avx2(unsigned int         sites,
                                           unsigned int         rate_cats,
                                           double *             parent_clv,
                                           unsigned int *       parent_scaler,
                                           const unsigned char *left_tipchar,
                                           const double *       right_clv,
                                           const double *       left_matrix,
                                           const double *       right_matrix,
                                           const unsigned int * right_scaler,
                                           const corax_state_t *tipmap,
                                           unsigned int         tipmap_size,
                                           unsigned int         attrib);

  CORAX_EXPORT void
  corax_core_update_clv_ii_avx2(unsigned int        states,
                                unsigned int        sites,
                                unsigned int        rate_cats,
                                double *            parent_clv,
                                unsigned int *      parent_scaler,
                                const double *      left_clv,
                                const double *      right_clv,
                                const double *      left_matrix,
                                const double *      right_matrix,
                                const unsigned int *left_scaler,
                                const unsigned int *right_scaler,
                                unsigned int        attrib);

  CORAX_EXPORT void
  corax_core_update_clv_repeats_generic_avx2(unsigned int        states,
                                             unsigned int        parent_sites,
                                             unsigned int        left_sites,
                                             unsigned int        right_sites,
                                             unsigned int        rate_cats,
                                             double *            parent_clv,
                                             unsigned int *      parent_scaler,
                                             const double *      left_clv,
                                             const double *      right_clv,
                                             const double *      left_matrix,
                                             const double *      right_matrix,
                                             const unsigned int *left_scaler,
                                             const unsigned int *right_scaler,
                                             const unsigned int *parent_id_site,
                                             const unsigned int *left_site_id,
                                             const unsigned int *right_site_id,
                                             double *            bclv_buffer,
                                             unsigned int        attrib);

#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CORAX_KERNEL_AVX2_CORE_CLV_H_ */
