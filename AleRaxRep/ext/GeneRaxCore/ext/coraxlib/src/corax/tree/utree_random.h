#ifndef CORAX_TREE_UTREE_RANDOM_H_
#define CORAX_TREE_UTREE_RANDOM_H_

#include "corax/tree/utree.h"

#ifdef __cplusplus
extern "C"
{
#endif

  CORAX_EXPORT corax_utree_t *corax_utree_random_create(unsigned int taxa_count,
                                                        const char *const *names,
                                                        unsigned int random_seed);

  CORAX_EXPORT int corax_utree_random_extend(corax_utree_t *    tree,
                                             unsigned int       ext_taxa_count,
                                             const char *const *ext_names,
                                             unsigned int       random_seed);

  CORAX_EXPORT corax_utree_t *
               corax_utree_random_resolve_multi(const corax_utree_t *multi_tree,
                                                unsigned int         random_seed,
                                                int *                clv_index_map);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CORAX_TREE_UTREE_RANDOM_H_ */
