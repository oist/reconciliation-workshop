#include "corax/corax.h"
#include "hashtable.h"

static void dealloc_graph_recursive(corax_unode_t *node)
{
  if (node->label) free(node->label);

  if (!node->next)
  {
    free(node);
    return;
  }

  corax_unode_t *sibling = node->next;
  while (node != sibling)
  {
    corax_unode_t *cur_node = sibling;
    dealloc_graph_recursive(cur_node->back);
    sibling = sibling->next;
    free(cur_node);
  }

  free(node);
}

static int
is_subsplit(corax_split_t child, corax_split_t parent, unsigned int split_len)
{
  unsigned int i;
  for (i = 0; i < split_len; ++i)
  {
    if ((child[i] & parent[i]) != child[i]) return 0;
  }
  return 1;
}

static void reverse_split(corax_split_t split, unsigned int tip_count)
{
  unsigned int split_size   = sizeof(corax_split_base_t) * 8;
  unsigned int split_offset = tip_count % split_size;
  unsigned int split_len    = tip_count / split_size + (split_offset > 0);
  unsigned int i;

  for (i = 0; i < split_len; ++i) split[i] = ~split[i];

  // there are no bits to mask and the mask will overflow
  if (!split_offset) {
    return;
  }

  // set all bits that do not represent any tip to 0
  corax_split_base_t mask = (1u << split_offset) - 1;
  split[split_len - 1] &= mask;
}

static corax_split_t clone_split(const corax_split_t from,
                                 unsigned int        split_len)
{
  corax_split_t to =
      (corax_split_t)calloc(split_len, sizeof(corax_split_base_t));
  memcpy(to, from, sizeof(corax_split_base_t) * split_len);

  return to;
}

/* reverse sort splits by weight */
static int sort_by_weight(const void *a, const void *b)
{
  double ca, cb;

  ca = (*((bitv_hash_entry_t **)a))->support;
  cb = (*((bitv_hash_entry_t **)b))->support;

  if (ca == cb)
    return 0;
  else
    return ((ca < cb) ? 1 : -1);
}

static void mre(bitv_hashtable_t *    h,
                corax_split_system_t *consensus,
                unsigned int          split_len,
                unsigned int          max_splits)
{
  bitv_hash_entry_t **split_list;

  unsigned int i = 0, j = 0;

  /* queue all splits */

  split_list = (bitv_hash_entry_t **)malloc(sizeof(bitv_hash_entry_t *)
                                            * h->entry_count);

  j = 0;
  for (i = 0; i < h->table_size; i++) /* copy hashtable h to list sbw */
  {
    bitv_hash_entry_t *e = h->table[i];
    while (e != NULL)
    {
      split_list[j] = e;
      ++j;
      e = e->next;
    }
  }
  assert(h->entry_count == j);

  /* sort by weight descending */
  qsort(
      split_list, h->entry_count, sizeof(bitv_hash_entry_t *), sort_by_weight);

  for (i = 0; (i < h->entry_count) && (consensus->split_count < max_splits);
       i++)
  {
    int                compatible      = 1;
    bitv_hash_entry_t *split_candidate = split_list[i];
    for (j = consensus->split_count; j > 0; --j)
    {
      corax_split_t split_consolidated = consensus->splits[j - 1];
      if (!corax_utree_split_compatible(split_candidate->bit_vector,
                                        split_consolidated,
                                        split_len,
                                        h->bit_count))
      {
        compatible = 0;
        break;
      }
    }

    if (compatible)
    {
      consensus->splits[consensus->split_count] =
          clone_split(split_candidate->bit_vector, split_len);
      consensus->support[consensus->split_count] = split_candidate->support;
      ++(consensus->split_count);
    }
  }

  free(split_list);

  return;
}

static int get_split_id(corax_split_t split, unsigned int split_len)
{
  unsigned int i, base_id, ctz;
  unsigned int taxa_per_split = 8 * sizeof(corax_split_base_t);
  int          id             = -1;
  //  unsigned int n_bits = setbit_count(split, split_len);
  unsigned int n_bits =
      bitv_popcount(split, taxa_per_split * split_len, split_len);

  if (n_bits != 1)
  {
    corax_set_error(CORAX_ERROR_INVALID_SPLIT, "Invalid trivial split");
    return -1;
  }

  for (i = 0; i < split_len; ++i)
  {
    if (split[i])
    {
      base_id = i * taxa_per_split;
      ctz     = __builtin_ctz(split[i]);
      assert(ctz < taxa_per_split);
      id = base_id + ctz;
      break;
    }
  }

  /* if the assertion below fails, there is an error either in this algorithm,
     or in setbit_count. */
  assert(id != -1);

  return id;
}

static void build_tips_recurse(corax_unode_t *    tree,
                               const char *const *tip_labels,
                               unsigned int       split_len)
{
  corax_consensus_data_t *data = (corax_consensus_data_t *)tree->data;
  corax_unode_t *         next_root;

  next_root = tree->next;
  if (next_root == tree)
  {
    assert(data);
    /* create tips */
    int tip_id = get_split_id(data->split, split_len);

    assert(tip_id != -1);

    if (tip_labels)
    {
      tree->label = (char *)malloc(strlen(tip_labels[tip_id]) + 1);
      strcpy(tree->label, tip_labels[tip_id]);
    }
    else
    {
      tree->label = NULL;
    }
    tree->node_index = (unsigned int)tip_id;
    tree->next       = NULL;
  }
  else
  {
    while (next_root != tree)
    {
      /* recurse next branch */
      build_tips_recurse(next_root->back, tip_labels, split_len);
      next_root = next_root->next;
    }
  }
}

static corax_unode_t *find_splitnode_recurse(corax_split_t  split,
                                             corax_unode_t *root,
                                             unsigned int   split_len)
{
  corax_unode_t *         next_root, *ret_node;
  corax_consensus_data_t *data = (corax_consensus_data_t *)root->data;

  if (is_subsplit(split, data->split, split_len))
  {
    /* check children */
    next_root = root->next;
    while (next_root != root)
    {
      if ((ret_node = find_splitnode_recurse(split, next_root->back, split_len))
          != NULL)
      {
        return ret_node;
      }
      next_root = next_root->next;
    }
    return root;
  }
  else
  {
    return NULL;
  }
}

static void connect_consensus_node(corax_unode_t *parent,
                                   corax_unode_t *child,
                                   unsigned int   split_len,
                                   int            auto_rearrange)
{
  corax_consensus_data_t *data_p, *data_c, *data_aux;
  corax_unode_t *         new_node, *aux_node, *aux_node2;
  data_p = (corax_consensus_data_t *)parent->data;
  data_c = (corax_consensus_data_t *)child->data;

  assert(is_subsplit(data_c->split, data_p->split, split_len));
  if (child->back)
  {
    new_node = child->back;

    /* disconnect from parent */
    aux_node = new_node;
    assert(new_node->next != new_node);
    while (aux_node->next != child->back) aux_node = aux_node->next;

    aux_node->next = new_node->next;
    new_node->next = 0;
  }
  else
  {
    new_node = (corax_unode_t *)malloc(sizeof(corax_unode_t));
  }

  if (auto_rearrange)
  {
    aux_node = parent->next;
    while (aux_node != parent)
    {
      aux_node2 = aux_node->next;
      data_aux  = (corax_consensus_data_t *)aux_node->back->data;
      if (is_subsplit(data_aux->split, data_c->split, split_len))
      {
        connect_consensus_node(child, aux_node->back, split_len, 0);
      }
      aux_node = aux_node2;
    }
  }

  /* connect new node */
  new_node->data = data_p;
  new_node->next = parent->next;
  parent->next   = new_node;

  /* connect child */
  new_node->back = child;
  child->back    = new_node;
}

static corax_unode_t *create_consensus_node(corax_unode_t *parent,
                                            corax_split_t  split,
                                            double         support,
                                            unsigned int   split_len)
{
  corax_unode_t *new_node      = (corax_unode_t *)malloc(sizeof(corax_unode_t));
  corax_consensus_data_t *data = 0;

  if (support > 0)
  {
    data = (corax_consensus_data_t *)malloc(sizeof(corax_consensus_data_t));
    data->split   = split;
    data->support = support;
  }

  new_node->data  = data;
  new_node->label = NULL;
  new_node->back  = NULL;

  /* self link */
  new_node->next = new_node;

  if (parent) { connect_consensus_node(parent, new_node, split_len, 1); }

  return new_node;
}

static void recursive_assign_indices(corax_unode_t *node,
                                     unsigned int * inner_clv_index,
                                     int *          inner_scaler_index,
                                     unsigned int * inner_node_index)
{
  if (!node) return;

  if (!node->next)
  {
    node->clv_index     = node->node_index;
    node->pmatrix_index = node->node_index;
    node->scaler_index  = CORAX_SCALE_BUFFER_NONE;
    return;
  }

  corax_unode_t *sibling = node->next;
  while (sibling != node)
  {
    recursive_assign_indices(
        sibling->back, inner_clv_index, inner_scaler_index, inner_node_index);
    sibling = sibling->next;
  }

  node->node_index             = *inner_node_index;
  node->next->node_index       = *inner_node_index + 1;
  node->next->next->node_index = *inner_node_index + 2;

  node->clv_index             = *inner_clv_index;
  node->next->clv_index       = *inner_clv_index;
  node->next->next->clv_index = *inner_clv_index;

  node->pmatrix_index             = *inner_clv_index;
  node->next->pmatrix_index       = node->next->back->pmatrix_index;
  node->next->next->pmatrix_index = node->next->next->back->pmatrix_index;

  node->scaler_index             = *inner_scaler_index;
  node->next->scaler_index       = *inner_scaler_index;
  node->next->next->scaler_index = *inner_scaler_index;

  *inner_clv_index    = *inner_clv_index + 1;
  *inner_scaler_index = *inner_scaler_index + 1;
  *inner_node_index   = *inner_node_index + 3;
}

static void reset_template_indices(corax_unode_t *node, unsigned int tip_count)
{
  unsigned int inner_clv_index    = tip_count;
  unsigned int inner_node_index   = tip_count;
  int          inner_scaler_index = 0;

  if (CORAX_UTREE_IS_TIP(node))
  {
    node = node->back;
    assert(!CORAX_UTREE_IS_TIP(node));
  }

  if (node->back)
    recursive_assign_indices(
        node->back, &inner_clv_index, &inner_scaler_index, &inner_node_index);

  corax_unode_t *sibling = node->next;
  while (sibling != node)
  {
    recursive_assign_indices(sibling->back,
                             &inner_clv_index,
                             &inner_scaler_index,
                             &inner_node_index);
    sibling = sibling->next;
  }

  node->node_index   = inner_node_index++;
  node->clv_index    = inner_clv_index;
  node->scaler_index = inner_scaler_index;

  sibling = node->next;
  while (sibling != node)
  {
    sibling->node_index   = inner_node_index;
    sibling->clv_index    = inner_clv_index;
    sibling->scaler_index = inner_scaler_index;
    inner_node_index++;

    if (sibling->back) sibling->pmatrix_index = sibling->back->pmatrix_index;

    sibling = sibling->next;
  }
}

static void consensus_data_destroy(void *data, int destroy_split)
{
  corax_consensus_data_t *cdata = (corax_consensus_data_t *)data;
  if (destroy_split) free(cdata->split);
  free(cdata);
}

static void fill_consensus_recurse(corax_consensus_utree_t *consensus_tree,
                                   corax_unode_t *          node,
                                   unsigned int *           cur_branch)
{
  assert(consensus_tree);
  assert(cur_branch);
  assert(node);

  corax_unode_t *child;
  unsigned int   max_degree = consensus_tree->tip_count;

  if (CORAX_UTREE_IS_TIP(node))
  {
    /* free tip data pointer */
    consensus_data_destroy(node->data, 1);

    /* unlink connected data pointer */
    node->back->data = 0;
    return;
  }

  child = node->next;
  while (child != node)
  {
    /* prevent infinite loop */
    --max_degree;
    assert(max_degree);

    fill_consensus_recurse(consensus_tree, child->back, cur_branch);
    child = child->next;
  }
  if (node->data)
  {
    memcpy(consensus_tree->branch_data + *cur_branch,
           node->data,
           sizeof(corax_consensus_data_t));
    consensus_data_destroy(node->data, 0);

    node->back->data = node->data = consensus_tree->branch_data + *cur_branch;
    ++(*cur_branch);
  }
}

static void fill_consensus(corax_consensus_utree_t *consensus_tree)
{
  assert(consensus_tree);
  assert(consensus_tree->tree);

  corax_unode_t *root       = consensus_tree->tree;
  unsigned int   cur_branch = 0;

  free(root->data);
  root->data = 0;

  unsigned int max_degree = consensus_tree->tip_count;
  fill_consensus_recurse(consensus_tree, root->back, &cur_branch);
  corax_unode_t *child = root->next;
  while (child != root)
  {
    /* prevent infinite loop */
    --max_degree;
    assert(max_degree);

    fill_consensus_recurse(consensus_tree, child->back, &cur_branch);
    child = child->next;
  }
}

/**
 * Check whether tip node indices in 2 trees are consistent to each other.
 *
 * Data pointers must contain the node index at first position
 * @param  t1 first tree
 * @param  t2 second tree
 *
 * @return CORAX_SUCCESS if node indices are consistent,
 *         CORAX_FAILURE otherwise (check corax_errmsg for details)
 */
CORAX_EXPORT int corax_utree_consistency_check(const corax_utree_t *t1,
                                               const corax_utree_t *t2)
{
  unsigned int    i;
  unsigned int    node_id;
  int             retval = CORAX_SUCCESS;
  corax_unode_t **tipnodes;
  char **         tipnames;
  unsigned int    tip_count;

  if (!t1 || !t2 || t1->tip_count != t2->tip_count)
  {
    corax_set_error(CORAX_ERROR_INVALID_TREE,
                    "Trees do not have the same number of tips\n");
    return CORAX_FAILURE;
  }

  tipnodes  = t1->nodes;
  tip_count = t1->tip_count;

  tipnames = (char **)malloc(tip_count * sizeof(char *));
  if (!tipnames)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Cannot allocate memory for tipnodes and tipnames\n");
    return CORAX_FAILURE;
  }

  /* fill names table */
  for (i = 0; i < tip_count; ++i)
  {
    node_id           = tipnodes[i]->node_index;
    tipnames[node_id] = tipnodes[i]->label;
  }

  /* check names consistency */
  tipnodes = t2->nodes;
  for (i = 0; i < tip_count; ++i)
  {
    node_id = tipnodes[i]->node_index;
    if (strcmp(tipnames[node_id], tipnodes[i]->label))
    {
      retval = CORAX_FAILURE;
      break;
    }
  }

  free(tipnames);
  return retval;
}

/**
 * Set t2 tip node indices consistent with t1.
 *
 * Data pointers must contain the node index at first position
 * @param  t1 reference tree
 * @param  t2 second tree
 *
 * @return CORAX_SUCCESS if the operation was applied correctly,
 *         CORAX_FAILURE otherwise (check corax_errmsg for details)
 */
CORAX_EXPORT int corax_utree_consistency_set(corax_utree_t *t1,
                                             corax_utree_t *t2)
{
  unsigned int    i, j;
  unsigned int    node_id;
  int             retval = CORAX_SUCCESS, checkval;
  corax_unode_t **tipnodes;
  char **         tipnames;
  unsigned int    tip_count;

  if (!t1 || !t2 || t1->tip_count != t2->tip_count)
  {
    corax_set_error(CORAX_ERROR_INVALID_TREE,
                    "Trees do not have the same number of tips\n");
    return CORAX_FAILURE;
  }

  tipnodes  = t1->nodes;
  tip_count = t1->tip_count;

  tipnames = (char **)malloc(tip_count * sizeof(char *));
  if (!tipnames)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Cannot allocate memory for tipnames\n");
    return CORAX_FAILURE;
  }

  /* fill names table */
  for (i = 0; i < tip_count; ++i)
  {
    node_id           = tipnodes[i]->node_index;
    tipnames[node_id] = tipnodes[i]->label;
  }

  /* set names consistency */
  tipnodes = t2->nodes;
  for (i = 0; i < tip_count; ++i)
  {
    // node_id = get_utree_node_id(tipnodes[i]);
    corax_unode_t *tipnode = tipnodes[i];
    checkval               = 0;
    for (j = 0; j < tip_count; ++j)
    {
      if (!strcmp(tipnames[j], tipnode->label))
      {
        checkval            = 1;
        tipnode->node_index = j;
        break;
      }
    }
    if (!checkval)
    {
      retval = CORAX_FAILURE;
      break;
    }
  }

  free(tipnames);
  return retval;
}

CORAX_EXPORT unsigned int corax_utree_rf_distance(const corax_unode_t *t1,
                                                  const corax_unode_t *t2,
                                                  unsigned int tip_count)
{
  unsigned int rf_distance;

  /* reset corax_error */
  corax_errno = 0;

  /* split both trees */
  corax_split_t *s1 = corax_utree_split_create(t1, tip_count, NULL);
  corax_split_t *s2 = corax_utree_split_create(t2, tip_count, NULL);

  /* compute distance */
  rf_distance = corax_utree_split_rf_distance(s1, s2, tip_count);

  /* clean up */
  corax_utree_split_destroy(s1);
  corax_utree_split_destroy(s2);

  assert(rf_distance <= 2 * (tip_count - 3));
  return rf_distance;
}

/* Consensus */

CORAX_EXPORT corax_consensus_utree_t *
             corax_utree_from_splits(const corax_split_system_t *split_system,
                                     unsigned int                tip_count,
                                     const char *const *const    tip_labels)
{
  const corax_split_t *    splits     = split_system->splits;
  unsigned int             split_size = sizeof(corax_split_base_t) * 8;
  unsigned int             split_len  = bitv_length(tip_count);
  unsigned int             i;
  corax_unode_t *          tree, *next_parent;
  corax_consensus_utree_t *return_tree;
  corax_split_t            rootsplit1, rootsplit2, next_split;
  double *                 support_values;
  corax_split_t *          all_splits;

  return_tree =
      (corax_consensus_utree_t *)malloc(sizeof(corax_consensus_utree_t));

  if (!return_tree)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Cannot allocate memory for consensus tree!");
    return NULL;
  }

  return_tree->tip_count    = tip_count;
  return_tree->branch_count = split_system->split_count;
  return_tree->branch_data  = (corax_consensus_data_t *)malloc(
      return_tree->branch_count * sizeof(corax_consensus_data_t));

  if (!return_tree->branch_data)
  {
    free(return_tree);
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Cannot allocate memory for consensus tree!");
    return NULL;
  }

  if (split_system->split_count == 0)
  {
    // build star tree
    rootsplit1 = (corax_split_t)calloc(split_len, sizeof(corax_split_base_t));
    rootsplit1[0] = 1;
    rootsplit2    = clone_split(rootsplit1, split_len);
    reverse_split(rootsplit1, tip_count);

    tree             = create_consensus_node(NULL, rootsplit1, 1.0, split_len);
    tree->back       = create_consensus_node(NULL, rootsplit2, 1.0, split_len);
    tree->back->back = tree;

    for (i = 1; i < tip_count; ++i)
    {
      unsigned int split_id = i / split_size;
      next_split = (corax_split_t)calloc(split_len, sizeof(corax_split_base_t));
      next_split[split_id] = (1 << (i % split_size));
      corax_unode_t *next_n =
          create_consensus_node(tree, next_split, 1.0, split_len);
      assert(next_n);
    }

    return_tree->tree = tree;
  }
  else
  {
    all_splits = (corax_split_t *)malloc((split_system->split_count + tip_count)
                                         * sizeof(corax_split_t));
    support_values = (double *)malloc((split_system->split_count + tip_count)
                                      * sizeof(double));
    memcpy(
        all_splits, splits, split_system->split_count * sizeof(corax_split_t));

    /* fill trivial splits */
    for (i = 0; i < tip_count; ++i)
    {
      support_values[split_system->split_count + i] = 1.0;
      all_splits[split_system->split_count + i] =
          (corax_split_t)calloc(split_len, sizeof(corax_split_base_t));
      {
        unsigned int split_id = i / split_size;
        all_splits[split_system->split_count + i][split_id] =
            (1 << (i % split_size));
      }
    }

    /* compute support for other splits */
    if (split_system->support)
      for (i = 0; i < split_system->split_count; ++i)
        support_values[i] =
            1.0 * split_system->support[i] / split_system->max_support;
    else
      for (i = 0; i < split_system->split_count; ++i)
        support_values[i] = 1.0 * split_system->max_support;

    /* create initial tree with 2 connected nodes of degree 1 */
    rootsplit1 = clone_split(all_splits[0], split_len);
    rootsplit2 = clone_split(all_splits[0], split_len);
    reverse_split(rootsplit2, tip_count);

    /* build tree out of the first split */
    tree =
        create_consensus_node(NULL, rootsplit1, support_values[0], split_len);
    tree->back =
        create_consensus_node(NULL, rootsplit2, support_values[0], split_len);
    tree->back->back = tree;

    return_tree->tree = tree;

    /* add splits individually */
    for (i = 1; i < (split_system->split_count + tip_count); ++i)
    {
      next_split = clone_split(all_splits[i], split_len);

      /* select branch */
      if (is_subsplit(next_split, rootsplit1, split_len))
        next_parent = tree;
      else if (is_subsplit(next_split, rootsplit2, split_len))
        next_parent = tree->back;
      else
      {
        reverse_split(next_split, tip_count);
        if (is_subsplit(next_split, rootsplit1, split_len))
          next_parent = tree;
        else if (is_subsplit(next_split, rootsplit2, split_len))
          next_parent = tree->back;
        else
        {
          corax_set_error(CORAX_ERROR_INVALID_SPLIT, "Splits are incompatible");
          free(return_tree->branch_data);
          free(return_tree);
          return_tree = NULL;
          break;
        }
      }

      /* select node */
      next_parent = find_splitnode_recurse(next_split, next_parent, split_len);
      assert(next_parent);

      /* create new node for the split*/
      create_consensus_node(
          next_parent, next_split, support_values[i], split_len);
    }

    /* clean */
    for (i = 0; i < tip_count; ++i)
      free(all_splits[split_system->split_count + i]);
    free(all_splits);
    free(support_values);
    free(rootsplit1);
  }

  build_tips_recurse(tree, tip_labels, split_len);
  if (tree->back) build_tips_recurse(tree->back, tip_labels, split_len);

  reset_template_indices(tree, tip_count);

  if (return_tree) fill_consensus(return_tree);

  /* return_tree == tree if success, or null if the algorithm failed */
  return return_tree;
}

CORAX_EXPORT corax_split_system_t *corax_utree_split_consensus(
    bitv_hashtable_t *splits_hash, unsigned int tip_count, double threshold)
{
  unsigned int          i;
  unsigned int          max_splits = tip_count - 3;
  unsigned int          split_len  = bitv_length(tip_count);
  corax_split_system_t *split_system;
  double                min_support = threshold;
  double                max_support = 1.0;
  double                thr_support = CORAX_MAX(min_support, .5 + CORAX_UTREE_WEIGHT_EPSILON);

  thr_support = CORAX_MIN(thr_support, max_support - CORAX_UTREE_WEIGHT_EPSILON);

  split_system =
      (corax_split_system_t *)calloc(1, sizeof(corax_split_system_t));

  if (!split_system)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Cannot allocate memory for split system!");
    return NULL;
  }

  split_system->splits =
      (corax_split_t *)calloc(max_splits, sizeof(corax_split_t));
  split_system->support     = (double *)calloc(max_splits, sizeof(double));
  split_system->split_count = 0;
  split_system->max_support = max_support;

  if (!split_system->splits || !split_system->support)
  {
    free(split_system);
    if (split_system->splits) free(split_system->splits);
    if (split_system->support) free(split_system->support);
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Cannot allocate memory for split system!");
    return NULL;
  }

  for (i = 0; i < splits_hash->table_size; ++i)
  {
    bitv_hash_entry_t * e     = splits_hash->table[i];
    bitv_hash_entry_t **e_ptr = &(splits_hash->table[i]);
    while (e != NULL)
    {
      int delete_split = 0;
      if (e->support > thr_support)
      {
        assert(split_system->split_count < max_splits);
        split_system->support[split_system->split_count] = e->support;
        split_system->splits[split_system->split_count] =
            clone_split(e->bit_vector, split_len);
        split_system->split_count++;
        delete_split = 1;
      }
      else
        delete_split = (min_support > 0.) && (e->support <= min_support);

      if (delete_split)
      {
        /* remove entry */
        hash_remove(splits_hash, e_ptr, e);

        e = *e_ptr;
      }
      else
      {
        e_ptr = &(e->next);
        e     = e->next;
      }
    }
  }

  if (min_support < .5 && splits_hash->entry_count > 0)
  {
    mre(splits_hash, split_system, split_len, max_splits);
  }

  return split_system;
}

/**
 * Build a consensus tree out of a list of trees and weights
 *
 * @param  treess           set of trees
 * @param  threshold        consensus threshold in [0,1].
 *                          1.0 -> strict
 *                          0.5 -> majority rule
 *                          0.0 -> extended majority rule
 * @param  weights          list of tree weights
 * @param  tree_count       number of trees
 * @return                  consensus unrooted tree structure
 */
CORAX_EXPORT corax_consensus_utree_t *
             corax_utree_weight_consensus(const corax_utree_t *const *trees,
                                          const double *              weights,
                                          double                      threshold,
                                          unsigned int                tree_count)
{

  const corax_utree_t *reference_tree =
      trees[0]; /* reference tree for consistency */
  corax_consensus_utree_t *consensus_tree = NULL; /* final consensus tree */
  corax_unode_t **         tipnodes;              /* tips from reference tree */
  bitv_hashtable_t *       splits_hash      = NULL;
  string_hashtable_t *     string_hashtable = NULL;
  corax_split_t *          tree_splits;
  unsigned int             i, j, tip_count = reference_tree->tip_count,
                     n_splits = tip_count - 3;

  /* validate threshold */
  if (threshold > 1 || threshold < 0)
  {
    corax_set_error(
        CORAX_ERROR_INVALID_THRESHOLD,
        "Invalid consensus threshold (%f). Should be in range [0.0,1.0]",
        threshold);
    return NULL;
  }

  /* validate weights */
  double sum_weights = 0.0;
  for (i = 0; i < tree_count; ++i)
  {
    if (weights[i] < 0)
    {
      sum_weights = 0;
      break;
    }
    sum_weights += weights[i];
  }
  if (fabs(1.0 - sum_weights) > CORAX_UTREE_WEIGHT_EPSILON)
  {
    corax_set_error(
        CORAX_ERROR_INVALID_TREE, "Invalid tree weights", threshold);
    return NULL;
  }

  /* store taxa names */
  tipnodes = reference_tree->nodes;

  string_hashtable = string_hash_init(10 * tip_count, tip_count);

  for (i = 0; i < tip_count; ++i)
  {
    string_hash_insert(
        tipnodes[i]->label, string_hashtable, (int)tipnodes[i]->node_index);
  }

  /* create hashtable */
  splits_hash = hash_init(tip_count * 10, tip_count);

  corax_errno = 0;

  for (i = 0; i < tree_count; ++i)
  {
    const corax_utree_t *current_tree = trees[i];
    if (current_tree->tip_count != tip_count)
    {
      /*TODO: error */
      return NULL;
    }

    tree_splits = corax_utree_split_create(
        current_tree
            ->nodes[current_tree->tip_count + current_tree->inner_count - 1],
        tip_count,
        NULL);

    /* insert normalized splits */
    for (j = 0; j < n_splits; ++j)
    {
      bitv_normalize(tree_splits[j], tip_count);

      hash_insert(
          tree_splits[j], splits_hash, j, HASH_KEY_UNDEF, weights[i], 0);
    }
    corax_utree_split_destroy(tree_splits);
  }

  if (corax_errno)
  {
    /* cleanup and spread error */
    char *aux_errmsg = (char *)malloc(strlen(corax_errmsg) + 1);
    strcpy(aux_errmsg, corax_errmsg);
    corax_set_error(corax_errno, "%s [tree #%u]", aux_errmsg, i);
    free(aux_errmsg);
    string_hash_destroy(string_hashtable);
    hash_destroy(splits_hash);
    return NULL;
  }

  /* build final split system */
  corax_split_system_t *split_system =
      corax_utree_split_consensus(splits_hash, tip_count, threshold);

  /* buld tree from splits */
  consensus_tree = corax_utree_from_splits(
      split_system, tip_count, (const char **)string_hashtable->labels);
  /*
   * The cast here is to silence a warning about
   *   char** -> const char *const *const.
   * This should be a safe conversion, but gcc and clange will warn about it
   * anyways. The problem is that pointers to pointers can be weird in the case
   * of concurrent code. However, this part should be fine, because it had no
   * const before.
   */

  /* cleanup */
  string_hash_destroy(string_hashtable);
  hash_destroy(splits_hash);
  corax_utree_split_system_destroy(split_system);
  return consensus_tree;
}

/**
 * Build a consensus tree out of a set of trees in a file in NEWICK format
 *
 * @param  trees_filename   trees filename
 * @param  threshold        consensus threshold in [0,1].
 *                          1.0 -> strict
 *                          0.5 -> majority rule
 *                          0.0 -> extended majority rule
 * @param[out] tree_count   number of trees parsed
 * @return                  consensus unrooted tree structure
 */
#if 0

static int cb_set_indices(corax_unode_t * node, void *data)
{
  string_hashtable_t * string_hashtable = (string_hashtable_t *) data;

  if (!node->next)
  {
    int index = string_hash_lookup(node->label, string_hashtable);
    if (index == -1)
      return 0;
    node->node_index = node->clv_index = (unsigned int) index;
  }

  return 1;
}

#define FCHUNK_LEN 2000

static corax_split_t * read_splits(FILE * file,
                                 int * n_chunks,
                                 unsigned int * tip_count,
                                 string_hashtable_t * string_hashtable)
{
  int line_size = *n_chunks * FCHUNK_LEN;
  char * tree_str = (char *) malloc((size_t)line_size);
  char * tree_str_ptr = tree_str;
  char * read;
  corax_split_t * splits;

  while ((read = fgets(tree_str_ptr, line_size, file)) &&
         tree_str_ptr[strlen(tree_str_ptr) - 1] != '\n')
  {
    /* increase line size */
    ++*n_chunks;
    tree_str = (char *) realloc(tree_str, (size_t)*n_chunks * FCHUNK_LEN);
    tree_str_ptr = tree_str + strlen(tree_str);
    line_size = FCHUNK_LEN;
  }

  if (read == NULL)
  {
    splits = NULL;
  }
  else
  {
    splits = corax_utree_split_newick_string(tree_str,
                                             *tip_count,
                                             string_hashtable);
  }
  free(tree_str);

  return splits;
}

static corax_utree_t * read_tree(FILE * file,
                               int * n_chunks,
                               unsigned int * tip_count,
                               string_hashtable_t * string_hashtable)
{
  int line_size = *n_chunks * FCHUNK_LEN;
  char * tree_str = (char *) malloc((size_t)line_size);
  char * tree_str_ptr = tree_str;
  char * read;
  corax_utree_t * utree;
  corax_unode_t * root;

  while ((read = fgets(tree_str_ptr, line_size, file)) &&
         tree_str_ptr[strlen(tree_str_ptr) - 1] != '\n')
  {
    /* increase line size */
    ++*n_chunks;
    tree_str = (char *) realloc(tree_str, (size_t)*n_chunks * FCHUNK_LEN);
    tree_str_ptr = tree_str + strlen(tree_str);
    line_size = FCHUNK_LEN;
  }

  if (read == NULL)
  {
    utree = NULL;
  }
  else
  {
    utree = corax_utree_parse_newick_string(tree_str);
    if (!utree)
    {
      return NULL;
    }
    *tip_count = utree->tip_count;
    root = utree->nodes[utree->tip_count + utree->inner_count - 1];

    if (root && string_hashtable)
    {
      if (!corax_utree_traverse_apply(root,
                                       NULL,
                                       NULL,
                                       &cb_set_indices,
                                       (void *)string_hashtable))
      {
        corax_set_error(
          CORAX_ERROR_INVALID_TREE,
          "Cannot match labels");
        corax_utree_destroy(utree, NULL);
        utree = NULL;
      }
    }
  }
  free(tree_str);

  return utree;
}

#undef FCHUNK_LEN

static FILE *get_number_of_trees(unsigned int *tree_count,
                                 const char *filename)
{
  FILE
    *f = fopen(filename, "r");

  if (!f)
    return NULL;

  unsigned int trees = 0;
  int ch;

  while((ch = fgetc(f)) != EOF)
    if(ch == ';')
      trees++;

  *tree_count = trees;

  rewind(f);

  return f;
}

CORAX_EXPORT corax_consensus_utree_t * corax_utree_consensus(
                                                const char * trees_filename,
                                                double threshold,
                                                unsigned int * _tree_count)
{
  FILE * trees_file;
  corax_utree_t * reference_tree = NULL; /* reference tree for consistency */
  corax_consensus_utree_t * consensus_tree = NULL; /* final consensus tree */
  corax_unode_t ** tipnodes;             /* tips from reference tree */
  bitv_hashtable_t * splits_hash = NULL;
  string_hashtable_t * string_hashtable = NULL;
  int  n_chunks = 1; /* chunks for reading newick trees */
  corax_split_t * tree_splits;
  unsigned int i,
               tip_count,
               n_splits,
               tree_count,         /* number of trees */
               current_tree_index; /* for error management */
  double individual_support;

  /* validate threshold */
  if (threshold > 1 || threshold < 0)
  {
    corax_set_error(
      CORAX_ERROR_INVALID_THRESHOLD,
      "Invalid consensus threshold (%f). Should be in range [0.0,1.0]",
      threshold);
    return NULL;
  }

  /* open file */
  if (!(trees_file = get_number_of_trees(&tree_count, trees_filename)))
  {
    corax_set_error(CORAX_ERROR_FILE_OPEN, "Cannot open trees file" );
    return NULL;
  }

  individual_support = 1.0/tree_count;

  if (_tree_count)
    *_tree_count = tree_count;

  /* read first tree */
  reference_tree = read_tree(trees_file, &n_chunks, &tip_count, NULL);

  if(!reference_tree)
  {
    assert(corax_errno);
    fclose(trees_file);
    return NULL;
  }

  /* store taxa names */
  tipnodes = reference_tree->nodes;

  string_hashtable = string_hash_init(10 * tip_count, tip_count);

  for (i=0; i<tip_count; ++i)
  {
    string_hash_insert(tipnodes[i]->label,
                       string_hashtable,
                       (int) tipnodes[i]->node_index);
  }

  /* create hashtable */
  splits_hash = hash_init(tip_count * 10, tip_count);

  corax_errno = 0;
  current_tree_index = 1;
  n_splits = tip_count - 3;

  tree_splits = corax_utree_split_create(reference_tree->nodes[
                                              reference_tree->tip_count +
                                              reference_tree->inner_count - 1],
                                          tip_count,
                                          NULL);

  while (tree_splits)
  {
    ++current_tree_index;

    /* insert normalized splits */
    for (i=0; i<n_splits; ++i)
    {
      bitv_normalize(tree_splits[i], tip_count);

      hash_insert(tree_splits[i],
                  splits_hash,
                  i,
                  HASH_KEY_UNDEF,
                  individual_support,
                  0);
    }
    corax_utree_split_destroy(tree_splits);

    /* parse next tree */
    tree_splits = read_splits(trees_file, &n_chunks, &tip_count, string_hashtable);
  }
  fclose(trees_file);

  if (corax_errno)
  {
    /* cleanup and spread error */
    char * aux_errmsg = (char *) malloc(strlen(corax_errmsg) + 1);
    strcpy(aux_errmsg, corax_errmsg);
    corax_set_error(corax_errno, "%s [tree #%u]", aux_errmsg, current_tree_index);
    free(aux_errmsg);
    string_hash_destroy(string_hashtable);
    hash_destroy(splits_hash);
    corax_utree_destroy(reference_tree, NULL);
    return NULL;
  }

  /* build final split system */
  corax_split_system_t * split_system = corax_utree_split_consensus(splits_hash,
                                                                   tip_count,
                                                                   threshold);

  /* buld tree from splits */
  consensus_tree = corax_utree_from_splits(split_system,
                                      tip_count,
                                      string_hashtable->labels);

  /* cleanup */
  string_hash_destroy(string_hashtable);
  hash_destroy(splits_hash);
  corax_utree_destroy(reference_tree, NULL);
  corax_utree_split_system_destroy(split_system);

  return consensus_tree;
}
#endif

CORAX_EXPORT void
corax_utree_split_system_destroy(corax_split_system_t *split_system)
{
  unsigned int i;

  for (i = 0; i < split_system->split_count; ++i) free(split_system->splits[i]);
  free(split_system->splits);
  free(split_system->support);
  free(split_system);
}

CORAX_EXPORT void corax_utree_consensus_destroy(corax_consensus_utree_t *tree)
{
  unsigned int i;

  /* dealloc branch data */
  for (i = 0; i < tree->branch_count; ++i) free(tree->branch_data[i].split);
  free(tree->branch_data);

  /* dealloc tree structure */
  dealloc_graph_recursive(tree->tree->back);
  dealloc_graph_recursive(tree->tree);

  /* dealloc tree */
  free(tree);
}
