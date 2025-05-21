#ifndef CORAX_TREE_UTREE_COMPARE_H_
#define CORAX_TREE_UTREE_COMPARE_H_

#include "utree_split.h"

#define CORAX_UTREE_WEIGHT_EPSILON 1e-12

typedef struct split_system_t
{
  corax_split_t *splits;
  double *       support;
  unsigned int   split_count;
  double         max_support;
} corax_split_system_t;

typedef struct consensus_data_t
{
  corax_split_t split;
  unsigned int  bit_count;
  double        support;
} corax_consensus_data_t;

typedef struct consensus_utree_t
{
  corax_unode_t *         tree;
  corax_consensus_data_t *branch_data;
  unsigned int            tip_count;
  unsigned int            branch_count;
} corax_consensus_utree_t;

#ifdef __cplusplus
extern "C"
{
#endif

  /* check that node ids and tip labels agree in both trees */
  CORAX_EXPORT int corax_utree_consistency_check(const corax_utree_t *t1,
                                                 const corax_utree_t *t2);

  /* if 2 different trees are parsed from newick node ids might have been set
     in a different order, so this function sets node ids in t2 such that
     node ids and tip labels agree in both trees */
  CORAX_EXPORT int corax_utree_consistency_set(corax_utree_t *t1,
                                               corax_utree_t *t2);

  CORAX_EXPORT unsigned int corax_utree_rf_distance(const corax_unode_t *t1,
                                                    const corax_unode_t *t2,
                                                    unsigned int tip_count);

  CORAX_EXPORT corax_consensus_utree_t *
               corax_utree_from_splits(const corax_split_system_t *split_system,
                                       unsigned int                tip_count,
                                       const char *const *const    tip_labels);

  CORAX_EXPORT corax_split_system_t *corax_utree_split_consensus(
      bitv_hashtable_t *splits_hash, unsigned int tip_count, double threshold);

  CORAX_EXPORT corax_consensus_utree_t *
               corax_utree_weight_consensus(const corax_utree_t *const *trees,
                                            const double *const         weights,
                                            double                      threshold,
                                            unsigned int                tree_count);

  // TODO: implement Newick->splits parser
  #if 0
  CORAX_EXPORT corax_consensus_utree_t * corax_utree_consensus(
                                                      const char * trees_filename,
                                                      double threshold,
                                                      unsigned int * tree_count);
  #endif

  CORAX_EXPORT void
  corax_utree_split_system_destroy(corax_split_system_t *split_system);

  CORAX_EXPORT void corax_utree_consensus_destroy(corax_consensus_utree_t *tree);

#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CORAX_TREE_UTREE_COMPARE_H_ */
