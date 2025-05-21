#ifndef CORAX_KERNEL_GENERIC_CORE_CLV_H_
#define CORAX_KERNEL_GENERIC_CORE_CLV_H_

#include "corax/core/common.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in core_clvs.c */

  CORAX_EXPORT void corax_core_create_lookup(unsigned int         states,
                                             unsigned int         rate_cats,
                                             double *             lookup,
                                             const double *       left_matrix,
                                             const double *       right_matrix,
                                             const corax_state_t *tipmap,
                                             unsigned int         tipmap_size,
                                             unsigned int         attrib);

  CORAX_EXPORT void
  corax_core_update_clv_tt(unsigned int         states,
                           unsigned int         sites,
                           unsigned int         rate_cats,
                           double *             parent_clv,
                           unsigned int *       parent_scaler,
                           const unsigned char *left_tipchars,
                           const unsigned char *right_tipchars,
                           unsigned int         tipmap_size,
                           const double *       lookup,
                           unsigned int         attrib);

  CORAX_EXPORT void corax_core_update_clv_ti(unsigned int         states,
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

  CORAX_EXPORT void corax_core_update_clv_ii(unsigned int        states,
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
  corax_core_update_clv_repeats(unsigned int        states,
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

  CORAX_EXPORT void
  corax_core_update_clv_repeats_generic(unsigned int        states,
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

  CORAX_EXPORT void
  corax_core_update_clv_repeatsbclv_generic(unsigned int        states,
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

  CORAX_EXPORT void corax_core_create_lookup_4x4(unsigned int  rate_cats,
                                                 double *      lookup,
                                                 const double *left_matrix,
                                                 const double *right_matrix);

  CORAX_EXPORT void
  corax_core_update_clv_tt_4x4(unsigned int         sites,
                               unsigned int         rate_cats,
                               double *             parent_clv,
                               unsigned int *       parent_scaler,
                               const unsigned char *left_tipchars,
                               const unsigned char *right_tipchars,
                               const double *       lookup,
                               unsigned int         attrib);

  CORAX_EXPORT void
  corax_core_update_clv_ti_4x4(unsigned int         sites,
                               unsigned int         rate_cats,
                               double *             parent_clv,
                               unsigned int *       parent_scaler,
                               const unsigned char *left_tipchars,
                               const double *       right_clv,
                               const double *       left_matrix,
                               const double *       right_matrix,
                               const unsigned int * right_scaler,
                               unsigned int         attrib);

#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CORAX_KERNEL_GENERIC_CORE_CLV_H_ */
