/**
 * @file utree.h
 *
 * @brief This header file contains functions related to the corax_utree_t
 * struct
 *
 * @author whoever
 */
#ifndef CORAX_TREE_UTREE_H_
#define CORAX_TREE_UTREE_H_

#include "corax/corax_core.h"

#define CORAX_UTREE_IS_TIP(node) (node->next == NULL)

/**
 * A structure that is a fundamental element of `corax_utree_t`. It contains a
 * next and back pointer. For more information, please see docs/corax_utree_t.md
 *
 */
typedef struct corax_unode_s
{
  /**
   * Label for the tree. Optional. If not present, then should be set to
   * `nullptr`
   */
  char *label;

  /**
   * Length of the edge, which is represented by the back pointer
   */
  double length;

  /**
   * Index of this node in the `nodes` buffer of `corax_utree_t`. Each
   * "super"-node shares and index. I.E. the index is on the "tree node" level,
   * not on the corax_unode_t level.
   */
  unsigned int node_index;

  /**
   * Index into the CLV buffer when computing a likelihood. For more
   * information, please see the documentation on `corax_partition_t`.
   */
  unsigned int clv_index;

  /**
   * Index into the scalar array to represent the CLV scaler. Please see the
   * documentation on `corax_partition_t` for more information.
   */
  int scaler_index;

  /**
   * Index into the array of probability matrices. These probability matrices
   * will be computed based on the branch length `length`. For more information
   * please see the documentation on `corax_partition_t`.
   */
  unsigned int pmatrix_index;

  /**
   * See the explaination in the concepts section of `docs/corax_utree_t.md`
   */
  struct corax_unode_s *next;

  /**
   * See the explaination in the concepts section of `docs/corax_utree_t.md`
   */
  struct corax_unode_s *back;

  /**
   * An extra pointer to store "user data". In praactice, this section can be
   * used for any task, but exsiting functions might use it, so be careful.
   */
  void *data;
} corax_unode_t;

/** @defgroup corax_utree_t corax_utree_t
 * Module for the `corax_utree_t` struct and associated functions
 */

/**
 * The data structure is made up of two different structs. The first,
 * corax_utree_t wraps the tree. In general, when a tree is used for a function,
 * it requires a corax_utree_t. Some important things to know about this
 * structure: the first inner_count nodes in the nodes array are assumed to be
 * "inner nodes". This means that they have a non-null next pointer. Several
 * functions that use corax_utree_ts don't check for this, so they may fail when
 * this assumption is violated. To avoid this, use the corax_utree_wraptree
 * function discussed below to create a corax_utree_t.
 *
 * @ingroup corax_utree_t
 */
typedef struct corax_utree_s
{
  /**
   * Number of tips in the tree
   */
  unsigned int tip_count;

  /**
   * Number of inner nodes. Not the number of `corax_unode_t` that make up the
   * tree, but the number of conceptual nodes on the phylogenetic tree.
   */
  unsigned int inner_count;

  /**
   * The number of edges in the tree
   */
  unsigned int edge_count;

  /**
   * Flag indicating if the tree is binary
   */
  int binary;

  /**
   * An array of `corax_unode_t` pointers
   */
  corax_unode_t **nodes;

  /**
   * A pointer to the virtual root. By convention, this is always an inner node.
   * All tree manipulation functions such as `corax_utree_wraptree` follow this
   * convention.
   */
  corax_unode_t *vroot;
} corax_utree_t;

#ifdef __cplusplus
extern "C"
{
#endif

  /**
   * Deallocate the memory associated with a utree. `cb_destroy` is used to delete
   * the user data allocated in data.
   *
   * @ingroup corax_utree_t
   */
  CORAX_EXPORT void corax_utree_destroy(corax_utree_t *tree,
                                        void (*cb_destroy)(void *));

  CORAX_EXPORT void corax_utree_reset_template_indices(corax_unode_t *node,
                                                       unsigned int   tip_count);

  CORAX_EXPORT void corax_utree_graph_destroy(corax_unode_t *root,
                                              void (*cb_destroy)(void *));

  /**
   * Takes a tree, represented by a node, and optionally a tip count. Will produce
   * a corax_utree_t that contains that tree. The pointer to the original node is
   * not invalidated.
   *
   * @param root Pointer to the virtual root. Should be an "inner node".
   *
   * @param tip_count Number of tips in contained in the tree represented by
   * `root`
   *
   * @ingroup corax_utree_t
   */
  CORAX_EXPORT corax_utree_t *corax_utree_wraptree(corax_unode_t *root,
                                                   unsigned int   tip_count);

  CORAX_EXPORT corax_utree_t *corax_utree_wraptree_multi(
      corax_unode_t *root, unsigned int tip_count, unsigned int inner_count);

  CORAX_EXPORT corax_unode_t *corax_utree_create_node(unsigned int clv_index,
                                                      int          scaler_index,
                                                      char *       label,
                                                      void *       data);

  CORAX_EXPORT int corax_unode_is_rooted(const corax_unode_t *root);

  CORAX_EXPORT int corax_utree_is_rooted(const corax_utree_t *tree);


  /**
   * Given the `corax_unode_t**` from a traversal using `corax_utree_traverse`,
   * this will create a list of `corax_operation_t`.
   *
   * @param trav_buffer
   *
   * @param trav_buffer_size
   *
   * @param[out] branches A buffer to store the branch length parameters used in
   * the operations. Optional.
   *
   * @param[out] pmatrix_indices A buffer to store the pmatrix indices used in the
   * operations. Optional.
   *
   * @param[out] ops Buffer to store the created ops. For a full traversal, the
   * allocated size should be equal to the number of the number of branches in the
   * tree. On an unrooted tree, this is `n-3`.
   *
   * @param[out] matrix_count Out parameter indicating the number matrices
   * required to perform the operations
   *
   * @param[out] ops_count Out parameter indicating the actual number of
   * operations.
   *
   * @ingroup corax_operation_t
   */
  CORAX_EXPORT void
  corax_utree_create_operations(const corax_unode_t *const *trav_buffer,
                                unsigned int                trav_buffer_size,
                                double *                    branches,
                                unsigned int *              pmatrix_indices,
                                corax_operation_t *         ops,
                                unsigned int *              matrix_count,
                                unsigned int *              ops_count);

  CORAX_EXPORT int corax_utree_check_integrity(const corax_utree_t *root);

  CORAX_EXPORT corax_unode_t *corax_utree_graph_clone(const corax_unode_t *root);

  /**
   * Clone a tree. This is a semi-deep copy. The fields `label` and the pointers
   * `next` and `back` are deep copied, but the data field is shallowly copied.
   *
   * @param root The tree to clone.
   */
  CORAX_EXPORT corax_utree_t *corax_utree_clone(const corax_utree_t *root);

  CORAX_EXPORT int corax_utree_set_clv_minimal(corax_unode_t *root,
                                               unsigned int   tip_count);


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CORAX_TREE_UTREE_H_ */
