#ifndef CORAX_KERNEL_GENERIC_CORE_FAST_PARSIMONY_H_
#define CORAX_KERNEL_GENERIC_CORE_FAST_PARSIMONY_H_

#include "corax/core/parsimony.h"
#include "corax/core/partition.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in fast_parsimony.c */

  CORAX_EXPORT corax_parsimony_t *
               corax_fastparsimony_init(const corax_partition_t *partition);

  CORAX_EXPORT void
  corax_fastparsimony_update_vectors(corax_parsimony_t *         parsimony,
                                     const corax_pars_buildop_t *ops,
                                     unsigned int                count);

  CORAX_EXPORT unsigned int
  corax_fastparsimony_root_score(const corax_parsimony_t *parsimony,
                                 unsigned int             root_index);

  CORAX_EXPORT unsigned int
  corax_fastparsimony_edge_score(const corax_parsimony_t *parsimony,
                                 unsigned int             node1_score_index,
                                 unsigned int             node2_score_index);

  CORAX_EXPORT void
  corax_fastparsimony_update_vector_4x4(corax_parsimony_t *         parsimony,
                                        const corax_pars_buildop_t *op);

  CORAX_EXPORT unsigned int
  corax_fastparsimony_edge_score_4x4(const corax_parsimony_t *parsimony,
                                     unsigned int             node1_score_index,
                                     unsigned int node2_score_index);

  CORAX_EXPORT void
  corax_fastparsimony_update_vector(corax_parsimony_t *         parsimony,
                                    const corax_pars_buildop_t *op);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CORAX_KERNEL_GENERIC_CORE_FAST_PARSIMONY_H_ */
