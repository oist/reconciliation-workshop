#ifndef CORAX_TREE_UTREE_OPS_H_
#define CORAX_TREE_UTREE_OPS_H_

#include "corax/tree/utree.h"

#define CORAX_TREE_ERROR_POLYPHYL_OUTGROUP 3970 // B + {10...}

#ifdef __cplusplus
extern "C"
{
#endif

  CORAX_EXPORT void corax_utree_set_length(corax_unode_t *edge, double length);

  CORAX_EXPORT void corax_utree_set_length_recursive(corax_utree_t *tree,
                                                     double         length,
                                                     int            missing_only);

  CORAX_EXPORT void corax_utree_scale_branches(corax_utree_t *tree,
                                               double branch_length_scaler);

  CORAX_EXPORT void corax_utree_scale_branches_all(corax_unode_t *root,
                                                   double branch_length_scaler);

  CORAX_EXPORT void
  corax_utree_scale_subtree_branches(corax_unode_t *root,
                                     double         branch_length_scaler);

  CORAX_EXPORT int corax_utree_collapse_branches(corax_utree_t *tree,
                                                 double         min_brlen);

  CORAX_EXPORT corax_unode_t *corax_utree_unroot_inplace(corax_unode_t *root);

  CORAX_EXPORT int corax_utree_root_inplace(corax_utree_t *tree);

  CORAX_EXPORT int corax_utree_outgroup_root(corax_utree_t *tree,
                                             unsigned int * outgroup_tip_ids,
                                             unsigned int   outgroup_size,
                                             int            add_root_node);

  CORAX_EXPORT int corax_utree_draw_support(const corax_utree_t *ref_tree,
                                            const double *       support,
                                            corax_unode_t **     node_map,
                                            char *(*cb_serialize)(double));

  CORAX_EXPORT corax_unode_t *corax_utree_serialize(corax_unode_t *tree,
                                                    unsigned int   tip_count);

  CORAX_EXPORT corax_utree_t *corax_utree_expand(corax_unode_t *serialized_tree,
                                                 unsigned int   tip_count);

#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CORAX_TREE_UTREE_OPS_H_ */
