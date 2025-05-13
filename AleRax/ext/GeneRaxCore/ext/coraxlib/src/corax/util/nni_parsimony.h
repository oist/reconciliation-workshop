#ifndef CORAX_UTIL_NNI_PARSIMONY_H_
#define CORAX_UTIL_NNI_PARSIMONY_H_

#include "corax/corax.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /* function in nni_parsimony.c */
  CORAX_EXPORT unsigned int
                corax_algo_nni_parsimony_local(corax_unode_t* q,
                                              corax_parsimony_t **list,
                                              unsigned int pars_score,
                                              unsigned int partition_count,
                                              unsigned int tips_count);

  CORAX_EXPORT unsigned int
               corax_algo_nni_round_parsimony(corax_utree_t * tree,
                                            corax_parsimony_t **list,
                                            unsigned int partition_count,
                                            unsigned int * score);

#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CORAX_UTIL_NNI_PARSIMONY_H_ */
