#ifndef CORAX_KERNEL_AVX_CORE_FAST_PARSIMONY_H_
#define CORAX_KERNEL_AVX_CORE_FAST_PARSIMONY_H_

#include "corax/core/parsimony.h"
#include "corax/core/partition.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in fast_parsimony_avx.c */

  CORAX_EXPORT void
  corax_fastparsimony_update_vector_4x4_avx(corax_parsimony_t *parsimony,
                                            const corax_pars_buildop_t *op);

  CORAX_EXPORT unsigned int
  corax_fastparsimony_edge_score_4x4_avx(const corax_parsimony_t *parsimony,
                                         unsigned int node1_score_index,
                                         unsigned int node2_score_index);

  CORAX_EXPORT void
  corax_fastparsimony_update_vector_avx(corax_parsimony_t *         parsimony,
                                        const corax_pars_buildop_t *op);

  CORAX_EXPORT unsigned int
  corax_fastparsimony_edge_score_avx(const corax_parsimony_t *parsimony,
                                     unsigned int             node1_score_index,
                                     unsigned int node2_score_index);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CORAX_KERNEL_AVX_CORE_FAST_PARSIMONY_H_ */
