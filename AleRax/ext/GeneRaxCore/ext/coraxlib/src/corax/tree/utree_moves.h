#ifndef CORAX_TREE_UTREE_MOVES_H_
#define CORAX_TREE_UTREE_MOVES_H_

#include "corax/tree/utree.h"

/* structures for handling topological rearrangement move rollbacks */

#define CORAX_TREE_REARRANGE_SPR 0
#define CORAX_TREE_REARRANGE_NNI 1
#define CORAX_TREE_REARRANGE_TBR 2

#define CORAX_UTREE_MOVE_NNI_LEFT 1
#define CORAX_UTREE_MOVE_NNI_RIGHT 2

/* error codes (for this module, 3000-4000) ; B = 2^10+2^11*/
/* TBR errors (B + {2^2,2^1,2^0}) */
#define CORAX_TREE_ERROR_TBR_LEAF_BISECTION 3073   // B + {001}
#define CORAX_TREE_ERROR_TBR_OVERLAPPED_NODES 3074 // B + {010}
#define CORAX_TREE_ERROR_TBR_SAME_SUBTREE 3075     // B + {011}
#define CORAX_TREE_ERROR_TBR_MASK 3079             // B + {111}

/* NNI errors (B + {2^4,2^3}) */
#define CORAX_TREE_ERROR_NNI_INVALID_MOVE 3080 // B + {01...}
#define CORAX_TREE_ERROR_NNI_LEAF 3081         // B + {01...}
#define CORAX_TREE_ERROR_NNI_MASK 3096         // B + {11...}

/* SPR errors (B + {2^6,2^5}) */
#define CORAX_TREE_ERROR_SPR_INVALID_NODE 3104 // B + {01...}
#define CORAX_TREE_ERROR_SPR_MASK 3168         // B + {11...}

typedef struct corax_utree_edge
{
  corax_unode_t *parent;
  corax_unode_t *child;
  double         length;
} corax_utree_edge_t;

typedef struct
{
  int    rearrange_type;
  int    rooted;
  double likelihood;

  union
  {
    struct
    {
      corax_unode_t *prune_edge;
      corax_unode_t *regraft_edge;
      double         prune_bl; //! length of the pruned branch
      double prune_left_bl;    //! length of the removed branch when pruning
      double prune_right_bl;   //! length of the removed branch when pruning
      double regraft_bl;       //! length of the splitted branch when regrafting
    } SPR;
    struct
    {
      corax_unode_t *edge;
      double         left_left_bl;
      double         left_right_bl;
      double         right_left_bl;
      double         right_right_bl;
      double         edge_bl;
      int            type;
    } NNI;
    struct
    {
      corax_unode_t     *bisect_edge;
      corax_utree_edge_t reconn_edge;
      double             bisect_left_bl;
      double             bisect_right_bl;
      double             reconn_parent_left_bl;
      double             reconn_parent_right_bl;
      double             reconn_child_left_bl;
      double             reconn_child_right_bl;
    } TBR;
  };
} corax_tree_rollback_t;


#ifdef __cplusplus
extern "C"
{
#endif

  CORAX_EXPORT int corax_utree_connect_nodes(corax_unode_t *parent,
                                             corax_unode_t *child,
                                             double         length);

  CORAX_EXPORT int corax_utree_bisect(corax_unode_t * edge,
                                      corax_unode_t **parent_subtree,
                                      corax_unode_t **child_subtree);

  CORAX_EXPORT corax_utree_edge_t
  corax_utree_reconnect(corax_utree_edge_t *edge, corax_unode_t *pruned_edge);

  CORAX_EXPORT corax_unode_t *corax_utree_prune(corax_unode_t *edge);

  CORAX_EXPORT int corax_utree_regraft(corax_unode_t *edge, corax_unode_t *tree);

  CORAX_EXPORT int corax_utree_interchange(corax_unode_t *edge1,
                                           corax_unode_t *edge2);

  CORAX_EXPORT int corax_utree_tbr(corax_unode_t *        b_edge,
                                   corax_utree_edge_t *   r_edge,
                                   corax_tree_rollback_t *rollback_info);

  CORAX_EXPORT int corax_utree_spr(corax_unode_t *        p_edge,
                                   corax_unode_t *        r_edge,
                                   corax_tree_rollback_t *rollback_info);

CORAX_EXPORT int corax_utree_nni(corax_unode_t         *edge,
                                 int                    type,
                                 corax_tree_rollback_t *rollback_info);
/**
 * Performs an NNI move on the tree. There are two options for moves here:
 *
 * - CORAX_NNI_NEXT,
 * - CORAX_NNI_NEXTNETX.
 *
 * This controls which of the two possible topology modifications that can be
 * done. If `CORAX_NNI_NEXT` is given, then the `edge->back->next` subtree is
 * swapped with `edge->next`. This means that `edge->back->next` becomes sister
 * to `edge`. If `CORAX_NNI_NEXTNEXT` is given, then `edge->back->next->next` is
 * used instead.
 *
 * The CLV, scaler and pmatrix indices are updated.
 *
 * @param[in] edge Node whos back pointer will be the "pivot" for the NNI move.
 *
 * @param[in] type One of either: `CORAX_NNI_NEXT` or `CORAX_NNI_NEXTNEXT`.
 *
 * @param[out] rollback_info Rollback information for undoing this move.
 *                           If it is NULL, rollback information is ignored.
 *
 * @return CORAX_SUCCESS if the move was applied correctly,
 *         CORAX_FAILURE otherwise (check corax_errmsg for details)
 *
 * @ingroup corax_utree_t
 */
  CORAX_EXPORT int corax_utree_nni(corax_unode_t *        edge,
                                   int                    type,
                                   corax_tree_rollback_t *rollback_info);

  CORAX_EXPORT int corax_tree_rollback(corax_tree_rollback_t *rollback_info);

#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CORAX_TREE_UTREE_MOVES_H_ */
