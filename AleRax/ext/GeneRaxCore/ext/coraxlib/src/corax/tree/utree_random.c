#include "corax/corax.h"

static void shuffle_tree_nodes(const corax_utree_t *tree, unsigned int seed)
{
  unsigned int        node_count = tree->tip_count + tree->inner_count;
  corax_random_state *rstate     = corax_random_create(seed);
  corax_unode_t **    subnodes =
      (corax_unode_t **)calloc(tree->tip_count, sizeof(corax_unode_t *));

  for (unsigned int i = tree->tip_count; i < node_count; ++i)
  {
    corax_unode_t *node   = tree->nodes[i];
    unsigned int   degree = 0;
    do {
      subnodes[degree] = node;
      degree++;
      node = node->next;
    } while (node != tree->nodes[i]);

    // Fisher–Yates shuffle
    for (unsigned int j = degree - 1; j > 0; --j)
    {
      unsigned int r = corax_random_getint(rstate, j + 1);
      CORAX_SWAP(subnodes[j], subnodes[r]);
    }

    // re-connect corax_unodes in the new, shuffled order
    tree->nodes[i] = node = subnodes[0];
    for (unsigned int j = 1; j < degree; ++j)
    {
      node->next = subnodes[j];
      node       = node->next;
    }

    // close roundabout
    node->next = tree->nodes[i];
  }

  corax_random_destroy(rstate);
  free(subnodes);
}

static void split_multi_node(corax_utree_t *tree,
                             corax_unode_t *first,
                             corax_unode_t *last,
                             unsigned int   degree)
{
  assert(last->next == first);
  if (degree > 3)
  {
    assert(first->next && first->next->next
           && first->next->next->next != first);

    // got a multifurcating node, split it in two
    unsigned int new_pmatrix_id = tree->edge_count;
    unsigned int new_node_id    = tree->edge_count * 2;
    unsigned int new_clv_id     = tree->tip_count + tree->inner_count;
    unsigned int new_scaler_id  = tree->inner_count;

    corax_unode_t *second = first->next;

    corax_unode_t *old_link = (corax_unode_t *)calloc(1, sizeof(corax_unode_t));
    corax_unode_t *new_link = (corax_unode_t *)calloc(1, sizeof(corax_unode_t));

    old_link->data = second->data;
    new_link->data = second->next->data;

    // close 'new' roundabout
    new_link->next = second->next;
    last->next     = new_link;

    // close 'old' roundabout
    old_link->next = first;
    second->next   = old_link;

    old_link->clv_index     = second->clv_index;
    old_link->scaler_index  = second->scaler_index;
    old_link->pmatrix_index = new_pmatrix_id;
    old_link->node_index    = new_node_id;

    assert(new_link->next && new_link->next->next);

    new_link->pmatrix_index = new_pmatrix_id;
    new_link->node_index    = new_node_id + 1;

    new_link->clv_index                    = new_link->next->clv_index =
        new_link->next->next->clv_index    = new_clv_id;
    new_link->scaler_index                 = new_link->next->scaler_index =
        new_link->next->next->scaler_index = (int)new_scaler_id;

    // set backpointers old<->new
    corax_utree_connect_nodes(
        old_link, new_link, CORAX_TREE_DEFAULT_BRANCH_LENGTH);

    tree->nodes[tree->inner_count + tree->tip_count] = new_link;

    tree->edge_count++;
    tree->inner_count++;

    // split new node if needed
    split_multi_node(tree, new_link, last, degree - 1);
  }
}

static int utree_insert_tips_random(corax_unode_t **nodes,
                                    unsigned int    taxa_count,
                                    unsigned int    start_tip,
                                    unsigned int    random_seed)
{
  unsigned int i;
  unsigned int start_inner_count     = start_tip - 2;
  unsigned int start_branches        = 2 * start_tip - 3;
  unsigned int max_branches          = 2 * taxa_count - 3;
  unsigned int placed_branches_count = 0;
  unsigned int last_branch_id        = 0;

  corax_unode_t **    branches = NULL;
  corax_random_state *rstate   = NULL;

  branches = (corax_unode_t **)calloc(max_branches, sizeof(corax_unode_t *));

  if (!branches)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Cannot allocate memory for branches!");
    return CORAX_FAILURE;
  }

  rstate = corax_random_create(random_seed);

  if (!rstate)
  {
    free(branches);
    return CORAX_FAILURE;
  }

  // check pmatrix indices on tip branches
  for (i = 0; i < taxa_count; ++i)
    last_branch_id = CORAX_MAX(last_branch_id, nodes[i]->pmatrix_index);

  for (i = taxa_count; i < taxa_count + start_inner_count; ++i)
  {
    corax_unode_t *snode = nodes[i];
    do {
      if (snode->clv_index > snode->back->clv_index)
      {
        branches[placed_branches_count++] = snode;
        last_branch_id = CORAX_MAX(last_branch_id, snode->pmatrix_index);
      }
      snode = snode->next;
    } while (snode != nodes[i]);
  }
  assert(placed_branches_count == start_branches);

  for (i = start_tip; i < taxa_count; ++i)
  {
    /* take tips iteratively */
    corax_unode_t *next_tip   = nodes[i];
    corax_unode_t *next_inner = nodes[taxa_count + i - 2];

    /* select random branch from the tree */
    int rand_branch_id = corax_random_getint(rstate, placed_branches_count);
    corax_unode_t *next_branch = branches[rand_branch_id];

    /* connect tip to selected branch */
    corax_utree_connect_nodes(
        next_branch->back, next_inner, CORAX_TREE_DEFAULT_BRANCH_LENGTH);
    corax_utree_connect_nodes(
        next_branch, next_inner->next, CORAX_TREE_DEFAULT_BRANCH_LENGTH);
    corax_utree_connect_nodes(
        next_tip, next_inner->next->next, CORAX_TREE_DEFAULT_BRANCH_LENGTH);

    if (CORAX_UTREE_IS_TIP(next_inner->back))
    {
      next_inner->next->pmatrix_index = next_inner->next->back->pmatrix_index =
          ++last_branch_id;
    }
    else
    {
      next_inner->pmatrix_index = next_inner->back->pmatrix_index =
          ++last_branch_id;
    }

    /* store branches */
    branches[placed_branches_count++] = next_inner;
    branches[placed_branches_count++] = next_inner->next->next;
  }
  assert(placed_branches_count == max_branches);

  /* clean */
  free(branches);
  corax_random_destroy(rstate);

  return CORAX_SUCCESS;
}

/**
 * Extend a tree by inserting new taxa to randomly chosen branches
 */
CORAX_EXPORT int corax_utree_random_extend(corax_utree_t *    tree,
                                           unsigned int       ext_taxa_count,
                                           const char *const *ext_names,
                                           unsigned int       random_seed)
{
  unsigned int old_taxa_count  = tree->tip_count;
  unsigned int old_inner_count = tree->inner_count;
  unsigned int old_node_count  = old_taxa_count + old_inner_count;
  unsigned int new_taxa_count  = old_taxa_count + ext_taxa_count;
  unsigned int new_inner_count = old_inner_count + ext_taxa_count;
  unsigned int new_node_count  = new_taxa_count + new_inner_count;

  unsigned int last_clv_id     = 0;
  unsigned int last_pmatrix_id = 0;
  unsigned int last_node_id    = 0;
  int          next_scaler_id  = 0;

  unsigned int i;
  int          retval;

  corax_unode_t **old_nodes = tree->nodes;
  corax_unode_t **new_nodes =
      (corax_unode_t **)calloc(new_node_count, sizeof(corax_unode_t *));

  if (!new_nodes)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC, "Cannot allocate memory for nodes!");
    return CORAX_FAILURE;
  }

  // 1:1 mapping for old tips
  for (i = 0; i < old_taxa_count; ++i) new_nodes[i] = old_nodes[i];

  // copy old inner nodes and adjust clvs
  for (i = old_taxa_count; i < old_node_count; ++i)
  {
    unsigned int new_idx = i + ext_taxa_count;
    new_nodes[new_idx]   = old_nodes[i];
    corax_unode_t *snode = new_nodes[new_idx];
    assert(snode->next);
    do {
      snode->clv_index += ext_taxa_count;
      snode->node_index += ext_taxa_count;
      last_clv_id     = CORAX_MAX(last_clv_id, snode->clv_index);
      last_node_id    = CORAX_MAX(last_node_id, snode->node_index);
      next_scaler_id  = CORAX_MAX(next_scaler_id, snode->scaler_index);
      last_pmatrix_id = CORAX_MAX(last_pmatrix_id, snode->pmatrix_index);
      snode           = snode->next;
    } while (snode != new_nodes[new_idx]);
  }

  // create new tip nodes
  for (i = old_taxa_count; i < new_taxa_count; ++i)
  {
    corax_unode_t *node = (corax_unode_t *)calloc(1, sizeof(corax_unode_t));
    node->clv_index     = i;
    node->node_index    = i;
    node->scaler_index  = CORAX_SCALE_BUFFER_NONE;
    node->pmatrix_index = ++last_pmatrix_id; // ????

    node->label = ext_names ? strdup(ext_names[i - old_taxa_count]) : NULL;

    new_nodes[i] = node;
  }

  // create new inner nodes
  for (i = old_node_count + ext_taxa_count; i < new_node_count; ++i)
  {
    corax_unode_t *node =
        corax_utree_create_node(++last_clv_id, ++next_scaler_id, NULL, NULL);

    node->node_index             = ++last_node_id;
    node->next->node_index       = ++last_node_id;
    node->next->next->node_index = ++last_node_id;

    new_nodes[i] = node;
  }

  retval = utree_insert_tips_random(
      new_nodes, new_taxa_count, old_taxa_count, random_seed);

  if (retval)
  {
    free(tree->nodes);
    tree->nodes       = new_nodes;
    tree->tip_count   = new_taxa_count;
    tree->inner_count = new_inner_count;
    tree->edge_count += 2 * ext_taxa_count;
    return CORAX_SUCCESS;
  }
  else
  {
    free(new_nodes);
    return CORAX_FAILURE;
  }
}

/**
 * Creates a random topology with default branch lengths
 */
CORAX_EXPORT corax_utree_t *corax_utree_random_create(unsigned int taxa_count,
                                                      const char *const *names,
                                                      unsigned int random_seed)
{
  /*
   * The algorithm works as follows:
   *    1. Build a minimal 3-tip tree
   *    2. Select a branch at random
   *    3. Connect next tip to that branch
   *    4. Repeat 2 and 3 until no tips left
   */
  unsigned int i;
  unsigned int tip_node_count   = taxa_count;
  unsigned int inner_node_count = taxa_count - 2;
  unsigned int node_count       = tip_node_count + inner_node_count;

  corax_unode_t **nodes =
      (corax_unode_t **)calloc(node_count, sizeof(corax_unode_t *));

  corax_unode_t *tree_root;

  corax_utree_t *wrapped_tree;

  unsigned int node_id = 0;

  /* allocate tips */
  for (i = 0; i < taxa_count; ++i)
  {
    nodes[i]                = (corax_unode_t *)calloc(1, sizeof(corax_unode_t));
    nodes[i]->clv_index     = i;
    nodes[i]->scaler_index  = CORAX_SCALE_BUFFER_NONE;
    nodes[i]->pmatrix_index = i;
    nodes[i]->node_index    = node_id++;

    if (names)
    {
      nodes[i]->label = (char *)malloc(strlen(names[i]) + 1);
      strcpy(nodes[i]->label, names[i]);
    }
    else
    {
      nodes[i]->label = NULL;
    }
  }

  /* allocate inner */
  for (i = taxa_count; i < node_count; ++i)
  {
    nodes[i] = corax_utree_create_node(i, (int)i, NULL, NULL);
    nodes[i]->scaler_index -= taxa_count;
    nodes[i]->next->scaler_index -= taxa_count;
    nodes[i]->next->next->scaler_index -= taxa_count;

    nodes[i]->node_index             = node_id++;
    nodes[i]->next->node_index       = node_id++;
    nodes[i]->next->next->node_index = node_id++;
  }
  assert(node_id == tip_node_count + inner_node_count * 3);

  /* set an inner node as return value */
  tree_root = nodes[taxa_count];

  /* build minimal tree with 3 tips and 1 inner node */
  corax_utree_connect_nodes(
      nodes[0], nodes[taxa_count], CORAX_TREE_DEFAULT_BRANCH_LENGTH);
  corax_utree_connect_nodes(
      nodes[1], nodes[taxa_count]->next, CORAX_TREE_DEFAULT_BRANCH_LENGTH);
  corax_utree_connect_nodes(nodes[2],
                            nodes[taxa_count]->next->next,
                            CORAX_TREE_DEFAULT_BRANCH_LENGTH);

  /* insert remaining taxa_count-3 tips into the tree */
  utree_insert_tips_random(nodes, taxa_count, 3, random_seed);

  /* clean */
  free(nodes);

  wrapped_tree = corax_utree_wraptree(tree_root, tip_node_count);
  return (wrapped_tree);
}

CORAX_EXPORT
corax_utree_t *corax_utree_random_resolve_multi(const corax_utree_t *multi_tree,
                                                unsigned int random_seed,
                                                int *        clv_index_map)
{
  if (!multi_tree)
  {
    corax_set_error(CORAX_ERROR_INVALID_PARAM, "Parameter multi_tree is NULL.");
    return NULL;
  }

  if (multi_tree->vroot->next
      && multi_tree->vroot->next->next == multi_tree->vroot)
  {
    corax_set_error(
        CORAX_ERROR_INVALID_TREE,
        "Unrooted tree is expected but a rooted tree was provided.");
    return NULL;
  }

  corax_utree_t *bin_tree = corax_utree_clone(multi_tree);

  unsigned int tip_count        = bin_tree->tip_count;
  unsigned int multi_node_count = bin_tree->tip_count + bin_tree->inner_count;
  unsigned int bin_node_count   = 2 * tip_count - 2;

  // 1:1 CLV index mapping for existing nodes
  if (clv_index_map)
  {
    for (unsigned int i = 0; i < multi_node_count; ++i)
    {
      const unsigned int clv_id = bin_tree->nodes[i]->clv_index;
      clv_index_map[clv_id]     = (int)clv_id;
    }
  }

  if (bin_tree->binary) return bin_tree;

  if (random_seed) shuffle_tree_nodes(bin_tree, random_seed);

  bin_tree->nodes = (corax_unode_t **)realloc(
      bin_tree->nodes, bin_node_count * sizeof(corax_unode_t *));

  // iterate over inner nodes, resolve multifurcations, and map new->old CLV
  // indices
  unsigned int old_inner_count = bin_tree->inner_count;
  for (unsigned int i = tip_count; i < multi_node_count; ++i)
  {
    corax_unode_t *start  = bin_tree->nodes[i];
    corax_unode_t *end    = NULL;
    corax_unode_t *snode  = start;
    unsigned int   degree = 0;
    do {
      end   = snode;
      snode = snode->next;
      degree++;
    } while (snode && snode != start);

    split_multi_node(bin_tree, start, end, degree);

    assert(bin_tree->inner_count == old_inner_count + degree - 3);
    if (clv_index_map)
    {
      for (unsigned int j = old_inner_count; j < bin_tree->inner_count; ++j)
      {
        unsigned int new_clv_id   = bin_tree->nodes[tip_count + j]->clv_index;
        clv_index_map[new_clv_id] = (int)start->clv_index;
      }
    }
    old_inner_count = bin_tree->inner_count;
  }

  assert(bin_tree->inner_count == bin_tree->tip_count - 2);

  bin_tree->binary = 1;

  /* re-assign node indices such that:
   * (1) all 3 corax_unode's of an inner node have consecutive indices: (x, x+1,
   * x+2) (2) for any two random multifurcation resolutions R1 and R2 holds (x,
   * x+1, x+2) in R1 iff (x, x+1, x+2) in R2
   */
  unsigned int max_node_index = tip_count;
  for (unsigned int i = tip_count; i < bin_node_count; ++i)
  {
    corax_unode_t *node          = bin_tree->nodes[i];
    node->node_index             = max_node_index++;
    node->next->node_index       = max_node_index++;
    node->next->next->node_index = max_node_index++;
  }
  assert(max_node_index == bin_tree->tip_count + 3 * bin_tree->inner_count);

  return bin_tree;
}
