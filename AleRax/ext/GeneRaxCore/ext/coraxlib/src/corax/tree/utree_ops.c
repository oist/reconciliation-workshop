#include "utree_ops.h"
#include "corax/corax.h"

static corax_unode_t *unode_prev(corax_unode_t *node)
{
  if (node->next)
  {
    corax_unode_t *prev = node;
    while (prev->next != node) prev = prev->next;
    return prev;
  }
  else
    return NULL;
}

/* This function removes a branch between lnode and lnode->back by
 * "dissolving" a roundabout ("inner node triplet") that contains lnode->back,
 * and merging its remainders into the "left" roundabout as show below:
 *
 *     *-l2-*          *-r2-*                  *-l2-*
      /      \        /      \                /      \
 * --l1  x1  l3------r1  x2  r3--   ---->  --l1  x1  r2--
 *    \      /        \      /                \      /
 *     *----*          *----*                  *-r3-*
 *
 *  where l3 = lnode, r1 = lnode->back
 */
static int remove_branch(corax_unode_t *lnode)
{
  corax_unode_t *rnode = lnode->back;

  /* can only remove a branch between two inner nodes */
  if (!lnode->next || !rnode->next) return CORAX_FAILURE;

  corax_unode_t *lnode_prev = unode_prev(lnode);
  corax_unode_t *lnode_next = lnode->next;
  corax_unode_t *rnode_prev = unode_prev(rnode);
  corax_unode_t *rnode_next = rnode->next;

  /* merge remaining subnodes of left and right nodes */
  lnode_prev->next = rnode_next;
  rnode_prev->next = lnode_next;

  /* update clv_index and scaler_index in right node remainder */
  while (rnode_next != lnode_next)
  {
    rnode_next->clv_index    = lnode_prev->clv_index;
    rnode_next->scaler_index = lnode_prev->scaler_index;
    rnode_next->label        = lnode_prev->label;
    rnode_next               = rnode_next->next;
  }

  /* destroy both subnodes adjacent to the deleted branch */
  free(lnode);
  free(rnode->label);
  free(rnode);

  return CORAX_SUCCESS;
}

static char *default_support_fmt(double support)
{
  char *sup_str;
  int   size_alloced = asprintf(&sup_str, "%lf", support);

  return size_alloced >= 0 ? sup_str : NULL;
}

CORAX_EXPORT void corax_utree_set_length(corax_unode_t *edge, double length)
{
  edge->length = edge->back->length = length;
}

CORAX_EXPORT void corax_utree_set_length_recursive(corax_utree_t *tree,
                                                   double         length,
                                                   int            missing_only)
{
  /* set branch lengths */
  unsigned int i;
  unsigned int tip_count   = tree->tip_count;
  unsigned int inner_count = tree->inner_count;
  for (i = 0; i < tip_count + inner_count; ++i)
  {
    corax_unode_t *node = tree->nodes[i];
    if (!node->length || !missing_only) corax_utree_set_length(node, length);
    if (node->next)
    {
      if (!node->next->length || !missing_only)
        corax_utree_set_length(node->next, length);
      if (!node->next->next->length || !missing_only)
        corax_utree_set_length(node->next->next, length);
    }
  }
}

CORAX_EXPORT void corax_utree_scale_branches(corax_utree_t *tree,
                                             double branch_length_scaler)
{
  /* scale branch lengths */
  unsigned int    i;
  unsigned int    tip_count   = tree->tip_count;
  unsigned int    inner_count = tree->inner_count;
  corax_unode_t **nodes       = tree->nodes;
  for (i = 0; i < tip_count; ++i) { nodes[i]->length *= branch_length_scaler; }
  for (i = tip_count; i < tip_count + inner_count; ++i)
  {
    nodes[i]->length *= branch_length_scaler;
    nodes[i]->next->length *= branch_length_scaler;
    nodes[i]->next->next->length *= branch_length_scaler;
  }
}

CORAX_EXPORT void corax_utree_scale_branches_all(corax_unode_t *root,
                                                 double branch_length_scaler)
{
  double root_length;

  /* scale all branches in a tree */
  corax_utree_scale_subtree_branches(root, branch_length_scaler);
  root_length = root->length;
  corax_utree_scale_subtree_branches(root->back, branch_length_scaler);

  /* undo duplicated scaling */
  root->length = root->back->length = root_length;
}

CORAX_EXPORT void
corax_utree_scale_subtree_branches(corax_unode_t *root,
                                   double         branch_length_scaler)
{
  /* scale all branches in a subtree rooted at node */
  root->length *= branch_length_scaler;
  root->back->length *= branch_length_scaler;

  if (root->next)
  {
    corax_utree_scale_subtree_branches(root->next->back, branch_length_scaler);
    corax_utree_scale_subtree_branches(root->next->next->back,
                                       branch_length_scaler);
  }
}

CORAX_EXPORT int corax_utree_collapse_branches(corax_utree_t *tree,
                                               double         min_brlen)
{
  if (!tree || !tree->vroot)
  {
    corax_set_error(CORAX_ERROR_INVALID_PARAM, "Empty tree specified!");
    return CORAX_FAILURE;
  }

  double        brlen_cutoff  = min_brlen + CORAX_ONE_EPSILON;
  unsigned int  tip_count     = tree->tip_count;
  unsigned int  inner_count   = tree->inner_count;
  unsigned int  node_count    = inner_count + tip_count;
  unsigned int  removed_count = 0;
  unsigned int *clv2pos_map =
      (unsigned int *)calloc(inner_count, sizeof(unsigned int));

  if (!clv2pos_map)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Cannot allocate memory for clv2pos map!");
    return CORAX_FAILURE;
  }

  /* to avoid making assumptions about node ordering in tree->nodes,
   * we build this map indexed based on clv_index */
  for (unsigned int i = tip_count; i < node_count; ++i)
  {
    unsigned int inner_clv_idx = tree->nodes[i]->clv_index - tip_count;
    clv2pos_map[inner_clv_idx] = i;
  }

  for (unsigned int i = tip_count; i < node_count; ++i)
  {
    corax_unode_t *node = tree->nodes[i];

    /* this node has been removed in a previous iteration -> skip */
    if (!node) continue;

    assert(!CORAX_UTREE_IS_TIP(node));

    corax_unode_t *start_node = NULL;
    do {
      corax_unode_t *anode = node->back;
      if (CORAX_UTREE_IS_TIP(anode) || node->length > brlen_cutoff)
      {
        if (!start_node) start_node = node;
        node = node->next;
      }
      else
      {
        /* remove branch and merge adjacent inner nodes */
        corax_unode_t *prev = unode_prev(node);
        if (tree->vroot == node || tree->vroot == anode) tree->vroot = prev;

        /* find out position of to-be-removed node in the tree->nodes array,
         * and earmark it for deletion by setting respective entry to NULL */
        unsigned int anode_pos = clv2pos_map[anode->clv_index - tip_count];
        assert(anode_pos >= tip_count && anode_pos < node_count);
        tree->nodes[anode_pos] = NULL;
        tree->nodes[i]         = prev;

        remove_branch(node);
        removed_count++;

        node = prev->next;
      }
    } while (node && node != start_node);
  }

  if (removed_count > 0)
  {
    /* compress tree->nodes array by excluding removed inner nodes */
    unsigned int idx            = tip_count;
    unsigned int new_node_count = node_count - removed_count;
    for (unsigned int i = tip_count; i < node_count; ++i)
    {
      corax_unode_t *node = tree->nodes[i];
      if (node) tree->nodes[idx++] = node;
    }
    assert(idx == new_node_count);

    /* update corax_utree_t metadata */
    tree->inner_count -= removed_count;
    tree->edge_count -= removed_count;
    tree->binary = 0;
    tree->nodes  = (corax_unode_t **)realloc(
        tree->nodes, new_node_count * sizeof(corax_unode_t *));
  }

  free(clv2pos_map);

  return CORAX_SUCCESS;
}

CORAX_EXPORT corax_unode_t *corax_utree_unroot_inplace(corax_unode_t *root)
{
  /* check for a bifurcation at the root */
  if (corax_unode_is_rooted(root))
  {
    if (root->next == root)
    {
      corax_set_error(CORAX_ERROR_NEWICK_SYNTAX,
                      "Unifurcation detected at root");
      return CORAX_FAILURE;
    }
    corax_unode_t *left  = root->back;
    corax_unode_t *right = root->next->back;

    if (root->label) free(root->label);
    free(root->next);
    free(root);

    double new_length = left->length + right->length;
    left->back        = right;
    right->back       = left;
    left->length = right->length = new_length;
    left->pmatrix_index          = right->pmatrix_index =
        CORAX_MIN(left->pmatrix_index, right->pmatrix_index);

    return left->next ? left : right;
  }
  else
    return root;
}

CORAX_EXPORT int corax_utree_root_inplace(corax_utree_t *tree)
{
  if (!tree)
  {
    corax_set_error(CORAX_ERROR_INVALID_PARAM, "Empty tree specified!");
    return CORAX_FAILURE;
  }

  /* check if tree is already rooted */
  if (tree->vroot->next && tree->vroot->next->next == tree->vroot)
    return CORAX_SUCCESS;

  corax_unode_t *root       = tree->vroot;
  corax_unode_t *root_back  = root->back;
  corax_unode_t *root_left  = (corax_unode_t *)calloc(1, sizeof(corax_unode_t));
  corax_unode_t *root_right = (corax_unode_t *)calloc(1, sizeof(corax_unode_t));
  root_left->next           = root_right;
  root_right->next          = root_left;
  double       root_brlen   = root->length / 2.;
  unsigned int last_clv_index     = 0;
  int          last_scaler_index  = 0;
  unsigned int last_node_index    = 0;
  unsigned int last_pmatrix_index = 0;
  unsigned int node_count         = tree->inner_count + tree->tip_count;

  for (unsigned int i = 0; i < node_count; ++i)
  {
    const corax_unode_t *node = tree->nodes[i];
    last_clv_index            = CORAX_MAX(last_clv_index, node->clv_index);
    last_scaler_index = CORAX_MAX(last_scaler_index, node->scaler_index);
    do {
      last_node_index    = CORAX_MAX(last_node_index, node->node_index);
      last_pmatrix_index = CORAX_MAX(last_pmatrix_index, node->pmatrix_index);
      node               = node->next;
    } while (node && node != tree->nodes[i]);
  }

  root_left->clv_index = root_right->clv_index = ++last_clv_index;
  root_left->scaler_index = root_right->scaler_index = ++last_scaler_index;
  root_left->node_index                              = ++last_node_index;
  root_right->node_index                             = ++last_node_index;
  root_right->pmatrix_index                          = ++last_pmatrix_index;

  corax_utree_connect_nodes(root, root_left, root_brlen);
  corax_utree_connect_nodes(root_right, root_back, root_brlen);

  tree->vroot = root_left;
  tree->inner_count++;
  tree->edge_count++;
  node_count++;

  tree->nodes                 = (corax_unode_t **)realloc(tree->nodes,
                                          node_count * sizeof(corax_unode_t *));
  tree->nodes[node_count - 1] = root_left;

  return CORAX_SUCCESS;
}

CORAX_EXPORT int corax_utree_outgroup_root(corax_utree_t *tree,
                                           unsigned int * outgroup_tip_ids,
                                           unsigned int   outgroup_size,
                                           int            add_root_node)
{
  corax_unode_t **split_to_node_map = NULL;
  corax_split_t * tree_splits       = NULL;
  corax_unode_t * new_root          = NULL;
  unsigned int    tip_count;
  unsigned int    split_count;

  if (!tree || !outgroup_tip_ids || !outgroup_size)
  {
    corax_set_error(CORAX_ERROR_INVALID_PARAM,
                    "Empty tree and/or outgroup specified!");
    return CORAX_FAILURE;
  }

  if (outgroup_size == 1)
  {
    // special case single-taxon outgroup: just find a tip by node_index
    for (unsigned int i = 0; i < tree->tip_count; ++i)
    {
      const corax_unode_t *node = tree->nodes[i];
      if (node->node_index == outgroup_tip_ids[0])
      {
        new_root = node->back;
        break;
      }
    }
  }
  else
  {
    tip_count   = tree->tip_count;
    split_count = tip_count - 3;

    split_to_node_map =
        (corax_unode_t **)calloc(split_count, sizeof(corax_unode_t *));

    if (!split_to_node_map)
    {
      corax_set_error(CORAX_ERROR_MEM_ALLOC,
                      "Cannot allocate memory for split->node map!");
      return CORAX_FAILURE;
    }

    tree_splits = corax_utree_split_create(
        tree->vroot, tree->tip_count, split_to_node_map);

    if (!tree_splits)
    {
      assert(corax_errno);
      free(split_to_node_map);
      return CORAX_FAILURE;
    }

    // create outgroup split
    corax_split_t outgroup_split =
        corax_utree_split_from_tips(outgroup_tip_ids, outgroup_size, tip_count);

    // check if this split is in the tree
    int root_idx =
        corax_utree_split_find(tree_splits, outgroup_split, tip_count);
    if (root_idx >= 0) new_root = split_to_node_map[root_idx];

    corax_utree_split_destroy(tree_splits);
    free(split_to_node_map);
    free(outgroup_split);
  }

  // set tree->vroot to the outgroup split node
  if (new_root)
  {
    tree->vroot = new_root;
    if (add_root_node)
      return corax_utree_root_inplace(tree);
    else
      return CORAX_SUCCESS;
  }
  else
  {
    corax_set_error(CORAX_TREE_ERROR_POLYPHYL_OUTGROUP,
                    "Outgroup is not monophyletic!");
    return CORAX_FAILURE;
  }
}

CORAX_EXPORT int corax_utree_draw_support(const corax_utree_t *ref_tree,
                                          const double *       support,
                                          corax_unode_t **     node_map,
                                          char *(*cb_serialize)(double))
{
  if (!ref_tree || !support)
  {
    corax_set_error(CORAX_ERROR_INVALID_PARAM, "Parameter is NULL!\n");
    return CORAX_FAILURE;
  }

  unsigned int split_count = ref_tree->edge_count - ref_tree->tip_count;
  for (size_t i = 0; i < split_count; ++i)
  {
    corax_unode_t *node =
        node_map ? node_map[i] : ref_tree->nodes[ref_tree->tip_count + i];

    /* this has to be an inner node! */
    assert(node->next);

    if (node->label) free(node->label);

    node->label = node->next->label = node->next->next->label =
        cb_serialize ? cb_serialize(support[i])
                     : default_support_fmt(support[i]);
  }

  return CORAX_SUCCESS;
}

/* auxiliary structure for the callback function below */
struct serial_tree_s
{
  corax_unode_t *serialized_tree;
  unsigned int   node_count;
  unsigned int   max_nodes;
};

/* callback function to fill the serialized tree */
static int cb_serialize(corax_unode_t *tree, void *data)
{
  struct serial_tree_s *list            = (struct serial_tree_s *)data;
  corax_unode_t *       serialized_tree = list->serialized_tree;
  unsigned int          cur_pos         = list->node_count;

  assert(cur_pos < list->max_nodes);

  memcpy(&(serialized_tree[cur_pos]), tree, sizeof(corax_unode_t));
  serialized_tree[cur_pos].data  = 0;
  serialized_tree[cur_pos].label = 0;
  if (!CORAX_UTREE_IS_TIP(tree))
  {
    /* set to arbitrary non-junk value */
    serialized_tree[cur_pos].next = (corax_unode_t *)1;
  }

  ++list->node_count;
  return 1;
}

// TODO: serialize/expand using a compressed format instead of corax_unode_t
CORAX_EXPORT corax_unode_t *corax_utree_serialize(corax_unode_t *tree,
                                                  unsigned int   tip_count)
{
  unsigned int         node_count;
  corax_unode_t *      serialized_tree;
  struct serial_tree_s data;

  node_count = 2 * tip_count - 2;

  /* allocate the serialized structure */
  serialized_tree = (corax_unode_t *)malloc(node_count * sizeof(corax_unode_t));

  /* fill data for callback function */
  data.serialized_tree = serialized_tree;
  data.node_count      = 0;
  data.max_nodes       = node_count;

  /* if tree is a tip, move to its back position */
  if (CORAX_UTREE_IS_TIP(tree)) tree = tree->back;

  /* apply callback function to serialize */
  corax_utree_traverse_apply(tree, NULL, NULL, cb_serialize, &data);

  if (data.node_count != data.max_nodes)
  {
    /* if the number of serialized nodes is not correct, return error */
    corax_set_error(CORAX_ERROR_INVALID_TREE,
                    "tree structure ot tip_count are invalid");
    free(serialized_tree);
    serialized_tree = NULL;
  }

  return serialized_tree;
}

CORAX_EXPORT corax_utree_t *corax_utree_expand(corax_unode_t *serialized_tree,
                                               unsigned int   tip_count)
{
  unsigned int    i, node_count, next_node_index;
  corax_unode_t **tree_stack;
  corax_unode_t * tree;
  unsigned int    tree_stack_top;

  corax_reset_error();

  node_count = 2 * tip_count - 2;

  /* allocate stack for at most 'n_tips' nodes */
  tree_stack = (corax_unode_t **)malloc(tip_count * sizeof(corax_unode_t *));
  tree_stack_top = 0;

  next_node_index = tip_count;

  /* read nodes */
  for (i = 0; i < node_count; ++i)
  {
    corax_unode_t *t   = 0;                  /* new node */
    corax_unode_t  t_s = serialized_tree[i]; /* serialized node */
    if (t_s.next)
    {
      /* build inner node and connect */
      corax_unode_t *t_cr, *t_r, *t_cl, *t_l;
      t   = corax_utree_create_node(t_s.clv_index,
                                  t_s.scaler_index,
                                  0,  /* label */
                                  0); /* data */
      t_l = t->next;
      t_r = t->next->next;

      t->node_index   = next_node_index++;
      t_r->node_index = next_node_index++;
      t_l->node_index = next_node_index++;

      /* pop and connect */
      t_cr       = tree_stack[--tree_stack_top];
      t_r->back  = t_cr;
      t_cr->back = t_r;
      t_cl       = tree_stack[--tree_stack_top];
      t_l->back  = t_cl;
      t_cl->back = t_l;

      /* set branch attributes */
      t->length          = t_s.length;
      t->pmatrix_index   = t_s.pmatrix_index;
      t_r->pmatrix_index = t_cr->pmatrix_index;
      t_r->length        = t_cr->length;
      t_l->pmatrix_index = t_cl->pmatrix_index;
      t_l->length        = t_cl->length;
    }
    else
    {
      t = (corax_unode_t *)calloc(1, sizeof(corax_unode_t));
      memcpy(t, &t_s, sizeof(corax_unode_t));
      assert(t->node_index < tip_count);
    }

    /* push */
    tree_stack[tree_stack_top++] = t;
  }

  /* root vertices must be in the stack */
  assert(tree_stack_top == 2);
  assert(next_node_index == (4 * tip_count - 6));

  tree             = tree_stack[--tree_stack_top];
  tree->back       = tree_stack[--tree_stack_top];
  tree->back->back = tree;

  if (tree->pmatrix_index != tree->back->pmatrix_index)
  {
    /* if pmatrix indices differ, connecting branch must be a tip */
    if (tree->back->next)
    {
      corax_set_error(CORAX_ERROR_INVALID_TREE,
                      "pmatrix indices do not match in serialized tree");
    }
    tree->pmatrix_index = tree->back->pmatrix_index;
  }

  if (tree->length != tree->back->length)
  {
    corax_set_error(CORAX_ERROR_INVALID_TREE,
                    "branch lengths do not matchin serialized tree");
  }

  if (corax_errno)
  {
    corax_utree_graph_destroy(tree, NULL);
    tree = 0;
  }

  free(tree_stack);

  return corax_utree_wraptree(tree, tip_count);
}
