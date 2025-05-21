#ifndef CORAX_KERNEL_GENERIC_WEIGHTED_PARSIMONY_H_
#define CORAX_KERNEL_GENERIC_WEIGHTED_PARSIMONY_H_

#include "corax/core/parsimony.h"
#include "corax/core/partition.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in weighted_parsimony.c */

  CORAX_EXPORT double
  corax_parsimony_build(corax_parsimony_t *         pars,
                        const corax_pars_buildop_t *operations,
                        unsigned int                count);

  CORAX_EXPORT void
  corax_parsimony_reconstruct(corax_parsimony_t *       pars,
                              const corax_state_t *     map,
                              const corax_pars_recop_t *operations,
                              unsigned int              count);

  CORAX_EXPORT double corax_parsimony_score(corax_parsimony_t *pars,
                                            unsigned int score_buffer_index);


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CORAX_KERNEL_GENERIC_WEIGHTED_PARSIMONY_H_ */
