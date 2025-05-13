#ifndef CORAX_TREE_UTREE_SPLIT_H_
#define CORAX_TREE_UTREE_SPLIT_H_

#include "corax/tree/utree.h"

typedef unsigned int        corax_split_base_t;
typedef corax_split_base_t *corax_split_t;

/* set of splits with equal dimensions, e.g. extracted from a tree */
typedef struct
{
  unsigned int tip_count;    /* number of taxa */
  unsigned int split_size;   /* size of a split storage element in bits */
  unsigned int split_len;    /* number of storage elements per split */
  unsigned int split_count;  /* number of splits currently set */
  corax_split_t * splits;      /* array of pointers to splits */
  int *id_to_split;          /* map between node/subnode ids and splits */
} corax_split_set_t;

typedef unsigned int hash_key_t;

typedef struct bitv_hash_entry
{
  hash_key_t    key;
  corax_split_t bit_vector;
  unsigned int *tree_vector;
  unsigned int  tip_count;
  double        support;
  unsigned int  bip_number;

  struct bitv_hash_entry *next;
} bitv_hash_entry_t;

typedef struct
{
  unsigned int        table_size;
  bitv_hash_entry_t **table;
  unsigned int        entry_count;
  unsigned int        bit_count; /* number of bits per entry */
  unsigned int        bitv_len;  /* bitv length */
} bitv_hashtable_t;

typedef struct string_hash_entry
{
  hash_key_t                key;
  int                       node_number;
  char *                    word;
  struct string_hash_entry *next;
} string_hash_entry_t;

typedef struct
{
  char **               labels;
  unsigned int          table_size;
  string_hash_entry_t **table;
  unsigned int          entry_count;
} string_hashtable_t;

#ifdef __cplusplus
extern "C"
{
#endif

  CORAX_EXPORT corax_split_t *
               corax_utree_split_create(const corax_unode_t *tree,
                                        unsigned int         tip_count,
                                        corax_unode_t **     split_to_node_map);

  CORAX_EXPORT corax_split_t
  corax_utree_split_from_tips(const unsigned int *subtree_tip_ids,
                              unsigned int        subtree_size,
                              unsigned int        tip_count);

  CORAX_EXPORT void corax_utree_split_normalize_and_sort(corax_split_t *s,
                                                         unsigned int   tip_count,
                                                         unsigned int   n_splits,
                                                         int keep_first);

  CORAX_EXPORT void corax_utree_split_show(const corax_split_t split,
                                           unsigned int        tip_count);

  CORAX_EXPORT void corax_utree_split_destroy(corax_split_t *split_list);

  CORAX_EXPORT unsigned int corax_utree_split_lightside(const corax_split_t split,
                                                        unsigned int tip_count);

  CORAX_EXPORT unsigned int corax_utree_split_hamming_distance(
      const corax_split_t s1, const corax_split_t s2, unsigned int tip_count);

  CORAX_EXPORT int corax_utree_split_compatible(const corax_split_t s1,
                                                const corax_split_t s2,
                                                unsigned int        split_len,
                                                unsigned int        tip_count);

  CORAX_EXPORT int corax_utree_split_find(const corax_split_t *split_list,
                                          const corax_split_t  split,
                                          unsigned int         tip_count);

  CORAX_EXPORT unsigned int corax_utree_split_rf_distance(const corax_split_t *s1,
                                                          const corax_split_t *s2,
                                                          unsigned int tip_count);

  /**
   * Extract all non-directed non-trivial splits, i.e. one split per inner tree edge.
   *
   * @return pointer to a newly created pll_split_set_t struct on success
   *         NULL on error
   */
  CORAX_EXPORT corax_split_set_t * corax_utree_splitset_create(const corax_utree_t * tree);

  /**
   * Extract all directed splits, i.e. two splits per every tree edge (inner+outer).
   *
   * @return pointer to a newly created pll_split_set_t struct on success
   *         NULL on error
   */
  CORAX_EXPORT corax_split_set_t * corax_utree_splitset_create_all(const corax_utree_t * tree);

  /**
   * Extract all directed splits from tree and store them in a pre-allocated split set.
   *
   * @param  split_set  existing split set with compatible dimensions
   *                    (use pllmod_utree_splitset_create_all() to initialize)
   * @param  tree       topology to extract splits from
   *
   * @return CORAX_SUCCESS if extraction and update was successful
   *         CORAX_FAILURE on error
   */
  CORAX_EXPORT int corax_utree_splitset_update_all(corax_split_set_t * split_set,
                                                   const corax_utree_t * tree);

  /**
   * Deallocate split set memory
   */
  CORAX_EXPORT void corax_utree_splitset_destroy(corax_split_set_t * split_set);

  
  /**
   * Check if a given topology is compatible with a set of constraint splits.
   *
   * @param  cons_splits  splits from a constraint tree
   * @param  tree         tree to be checked
   *
   * @return CORAX_SUCCESS if topology is compatible
   *         CORAX_FAILURE otherwise
   */
  CORAX_EXPORT int corax_utree_constraint_check_splits_tree(corax_split_set_t * cons_splits,
                                                            const corax_utree_t * tree);
  
  /**
   * Check if a given topology (query tree) is compatible with a topological constraint,
   * while both are specified as a set of splits. Essentially, it checks whether
   * every split in the constraint tree is present in the query tree.
   *
   * @param  cons_splits  splits from a constraint tree
   * @param  tree_splits  splits from a tree to be checked
   *
   * @return CORAX_SUCCESS if topology is compatible
   *         CORAX_FAILURE otherwise
   */
  CORAX_EXPORT int corax_utree_constraint_check_splits(corax_split_set_t * cons_splits, 
                                                       corax_split_set_t * tree_splits);
  
  /**
   * Check if a given topology is compatible with a topological constraint.
   *
   * @param  cons_tree  a constraint tree
   * @param  tree       tree to be checked
   *
   * @return CORAX_SUCCESS if topology is compatible
   *         CORAX_FAILURE otherwise
   */
  CORAX_EXPORT int corax_utree_constraint_check_tree(const corax_utree_t * cons_tree,
                                                     const corax_utree_t * tree);

  /**
   * Check if an NNI is compatible with a given topological constraint
   *
   * @param  cons_splits  splits from a constraint tree
   * @param  tree_splits  splits from the ORIGINAL tree BEFORE PRUNING
   * @param  edge         NNI edge
   * @param  nni_type     CORAX_UTREE_MOVE_NNI_LEFT or CORAX_UTREE_MOVE_NNI_RIGHT
   *
   * @return CORAX_SUCCESS if NNI is compatible with constraint
   *         CORAX_FAILURE otherwise
   */
  CORAX_EXPORT int corax_utree_constraint_check_nni(corax_split_set_t * cons_splits,
                                                    corax_split_set_t * tree_splits,
                                                    corax_unode_t * edge,
                                                    int nni_type);

  /**
   * Check if an SPR is compatible with a given topological constraint
   *
   * @param  cons_splits  splits from a constraint tree
   * @param  tree_splits  splits from the ORIGINAL tree BEFORE PRUNING
   * @param  p_edge  pruned subtree
   * @param  r_edge  re-instertion edge
   *
   * @return CORAX_SUCCESS if SPR is compatible with constraint
   *         CORAX_FAILURE otherwise
   */
  CORAX_EXPORT int corax_utree_constraint_check_spr(corax_split_set_t * cons_splits,
                                                    corax_split_set_t * tree_splits,
                                                    corax_unode_t * p_edge,
                                                    corax_unode_t * r_edge);

  /**
   * Check if constraint is relevant for a given subtree. For instance, a subtree comprising
   * only "free" taxa (e.g. those absent from the constraint tree) is not affected,
   * and hence constraint check is not required before regrafting this subtree.
   *
   * @param  cons_splits  splits from a constraint tree
   * @param  tree_splits  splits from the ORIGINAL tree BEFORE PRUNING
   * @param  p_edge       pruned subtree
   *
   * @return CORAX_SUCCESS subtree is affected by the constraint
   *         CORAX_FAILURE otherwise
   */
  CORAX_EXPORT int corax_utree_constraint_subtree_affected(const corax_split_set_t * cons_splits,
                                                           const corax_split_set_t * tree_splits,
                                                           corax_unode_t * p_edge);


  // TODO: implement Newick->splits parser
  #if 0
  CORAX_EXPORT corax_split_t * corax_utree_split_newick_string(char * s,
                                                         unsigned int tip_count,
                                                         string_hashtable_t * names_hash);
  #endif

  /* split hashtable */

  CORAX_EXPORT
  bitv_hashtable_t *corax_utree_split_hashtable_create(unsigned int tip_count,
                                                       unsigned int slot_count);

  CORAX_EXPORT bitv_hash_entry_t *corax_utree_split_hashtable_insert_single(
      bitv_hashtable_t *splits_hash, const corax_split_t split, double support);

  CORAX_EXPORT bitv_hashtable_t *
               corax_utree_split_hashtable_insert(bitv_hashtable_t *splits_hash,
                                                  corax_split_t *   splits,
                                                  unsigned int      tip_count,
                                                  unsigned int      split_count,
                                                  const double *    support,
                                                  int               update_only);

  CORAX_EXPORT bitv_hash_entry_t *
               corax_utree_split_hashtable_lookup(bitv_hashtable_t *  splits_hash,
                                                  const corax_split_t split,
                                                  unsigned int        tip_count);

  CORAX_EXPORT
  void corax_utree_split_hashtable_destroy(bitv_hashtable_t *hash);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CORAX_TREE_UTREE_SPLIT_H_ */
