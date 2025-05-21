#ifndef CORAX_UTIL_STEPWISE_H_
#define CORAX_UTIL_STEPWISE_H_

#include "corax/corax_core.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in stepwise.c */

  CORAX_EXPORT corax_utree_t *
               corax_fastparsimony_stepwise(corax_parsimony_t **list,
                                            const char *const * labels,
                                            unsigned int *      score,
                                            unsigned int        count,
                                            unsigned int        seed);

  CORAX_EXPORT int 
                corax_fastparsimony_stepwise_spr_round(corax_utree_t * tree,
                                                   corax_parsimony_t ** pars_list,
                                                   unsigned int pars_count,
                                                   const unsigned int * tip_msa_idmap,
                                                   unsigned int seed,
                                                   const int * clv_index_map,
                                                   unsigned int * cost);
    
    CORAX_EXPORT int 
                  corax_fastparsimony_stepwise_extend(corax_utree_t * tree,
                                                 corax_parsimony_t ** pars_list,
                                                 unsigned int pars_count,
                                                 char * const * labels,
                                                 const unsigned int * tip_msa_idmap,
                                                 unsigned int seed,
                                                 unsigned int * cost);

#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CORAX_UTIL_STEPWISE_H_ */
