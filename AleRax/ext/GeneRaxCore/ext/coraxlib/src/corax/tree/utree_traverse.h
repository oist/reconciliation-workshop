#ifndef CORAX_TREE_UTREE_TRAVERSE_H_
#define CORAX_TREE_UTREE_TRAVERSE_H_

#include "corax/corax_core.h"
#include "corax/tree/utree.h"

#define CORAX_TREE_TRAVERSE_POSTORDER 1
#define CORAX_TREE_TRAVERSE_PREORDER 2

#define CORAX_TREE_TRAVERSE_FULL              1
#define CORAX_TREE_TRAVERSE_PARTIAL           2
#define CORAX_TREE_TRAVERSE_NONE              3

#ifdef __cplusplus
extern "C"
{
#endif

  /**
   * Creates a list of nodes from a traversal, starting at root. The order of the
   * traversal can be controlled with the traversal argument, which accepts either
   * `CORAX_TREE_TRAVERSE_POSTORDER `or `CORAX_TREE_TRAVERSE_PREORDER`. The
   * callback function controls which nodes are traversed by returning true for a
   * node node which should be added to `outbuffer`. If false is returned instead,
   * then the traversal is halted for that subtree, and the node which returned it
   * is not added to `outbuffer`. By doing this the traversal of the tree can be
   * halted early, which is useful for partial traversals. `trav_size `is an out
   * parameter.  Returns `CORAX_SUCCESS` on a traversal without errors, and
   * `CORAX_FAILURE` if there was an error.
   *
   * If a full traversal is desired, the callback should return true for all
   * inputs, for example:
   *
   * ```
   * int trav_cb(corax_unode_t *){
   *   return CORAX_SUCCESS;
   * }
   * ```
   *
   * The intended use of this function is to build a traversal buffer, which is to
   * be used in later computations. For example, the `outbuffer `parameter is used
   * to build the operations for likelihood computations (see
   * `corax_utree_create_operations`). The function will fail if `nullptr `is
   * passed in for `outbuffer`. If the behavior of this function is needed, but
   * the buffer is not needed, please use `corax_utree_traverse_apply`.
   *
   * A full traversal of an unrooted tree will take `2n-3` elements. This means
   * that all the buffers here should be that long in the case of a full
   * traversal.
   *
   * @param root Node which the traversal will start from.
   *
   * @param traversal Type of traversal. Can be `CORAX_TREE_TRAVERSE_POSTORDER` or
   * `CORAX_TREE_TRAVERSE_PREORDER`
   *
   * @param cbtrav Callback function which can be used to control the traversal.
   * Please see the general function documentation for more detail.
   *
   * @param[out] outbuffer Buffer which will contain the traversed nodes in order.
   * Buffer should be allocated.
   *
   * @param[out] trav_size Pointer a buffer which will contain the final traversal
   * size.
   *
   * @return `CORAX_SUCCESS` on a traversal without errors. `CORAX_FAILURE`
   * otherwise.
   *
   * @ingroup corax_utree_t
   */
  CORAX_EXPORT int corax_utree_traverse(corax_unode_t *root,
                                        int            traversal,
                                        int (*cbtrav)(corax_unode_t *),
                                        corax_unode_t **outbuffer,
                                        unsigned int *  trav_size);

  CORAX_EXPORT int corax_utree_traverse_const(const corax_unode_t *root,
                                        int            traversal,
                                        int (*cbtrav)(const corax_unode_t *),
                                        corax_unode_t const ** outbuffer,
                                        unsigned int *  trav_size);

  CORAX_EXPORT int corax_utree_traverse_subtree(corax_unode_t *root,
                                                int            traversal,
                                                int (*cbtrav)(corax_unode_t *),
                                                corax_unode_t **outbuffer,
                                                unsigned int *  trav_size);
  CORAX_EXPORT int corax_utree_every(corax_utree_t *tree,
                                     int (*cb)(const corax_utree_t *,
                                               const corax_unode_t *));

  CORAX_EXPORT int corax_utree_every_const(const corax_utree_t *tree,
                                           int (*cb)(const corax_utree_t *tree,
                                                     const corax_unode_t *));

  CORAX_EXPORT int
  corax_utree_traverse_apply(corax_unode_t *root,
                             int (*cb_pre_trav)(corax_unode_t *, void *),
                             int (*cb_in_trav)(corax_unode_t *, void *),
                             int (*cb_post_trav)(corax_unode_t *, void *),
                             void *data);

  CORAX_EXPORT int corax_utree_nodes_at_node_dist(corax_unode_t * node,
                                                  corax_unode_t **outbuffer,
                                                  unsigned int *  node_count,
                                                  unsigned int    min_distance,
                                                  unsigned int    max_distance);

  CORAX_EXPORT int corax_utree_nodes_at_edge_dist(corax_unode_t * edge,
                                                  corax_unode_t **outbuffer,
                                                  unsigned int *  node_count,
                                                  unsigned int    min_distance,
                                                  unsigned int    max_distance);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CORAX_TREE_UTREE_TRAVERSE_H_ */
