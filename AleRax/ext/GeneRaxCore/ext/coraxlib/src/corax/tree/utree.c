/*
    Copyright (C) 2015-2021 Tomas Flouri, Alexey Kozlov

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
    Heidelberg Institute for Theoretical Studies,
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "utree.h"
#include "utree_traverse.h"

static void dealloc_data(corax_unode_t *node, void (*cb_destroy)(void *))
{
  if (node->data)
  {
    if (cb_destroy) cb_destroy(node->data);
  }
}

static void dealloc_graph_recursive(corax_unode_t *node,
                                    void (*cb_destroy)(void *),
                                    int level)
{
  if (!node->next)
  {
    /* tip node */
    dealloc_data(node, cb_destroy);
    free(node->label);
    free(node);
  }
  else
  {
    /* inner node */
    if (node->label) free(node->label);

    corax_unode_t *snode = node;
    do {
      if (node != snode || level == 0)
        dealloc_graph_recursive(snode->back, cb_destroy, level + 1);
      corax_unode_t *next = snode->next;
      dealloc_data(snode, cb_destroy);
      free(snode);
      snode = next;
    } while (snode && snode != node);
  }
}

CORAX_EXPORT void
corax_utree_create_operations(const corax_unode_t *const *trav_buffer,
                              unsigned int                trav_buffer_size,
                              double *                    branches,
                              unsigned int *              pmatrix_indices,
                              corax_operation_t *         ops,
                              unsigned int *              matrix_count,
                              unsigned int *              ops_count)
{
  const corax_unode_t *node;
  unsigned int         i;

  *ops_count = 0;
  if (matrix_count) *matrix_count = 0;

  for (i = 0; i < trav_buffer_size; ++i)
  {
    node = trav_buffer[i];

    /* if the current node is the second end-point of the edge
    shared with the root node, then do not add the edge to the
    list as it will be added in the end (avoid duplicate edges
    in the list) */
    if (node != trav_buffer[trav_buffer_size - 1]->back)
    {
      if (branches) *branches++ = node->length;
      if (pmatrix_indices) *pmatrix_indices++ = node->pmatrix_index;
      if (matrix_count) *matrix_count = *matrix_count + 1;
    }

    if (node->next)
    {
      ops[*ops_count].parent_clv_index    = node->clv_index;
      ops[*ops_count].parent_scaler_index = node->scaler_index;

      ops[*ops_count].child1_clv_index    = node->next->back->clv_index;
      ops[*ops_count].child1_scaler_index = node->next->back->scaler_index;
      ops[*ops_count].child1_matrix_index = node->next->back->pmatrix_index;

      ops[*ops_count].child2_clv_index = node->next->next->back->clv_index;
      ops[*ops_count].child2_scaler_index =
          node->next->next->back->scaler_index;
      ops[*ops_count].child2_matrix_index =
          node->next->next->back->pmatrix_index;

      *ops_count = *ops_count + 1;
    }
  }
}

/* a callback function for checking tree integrity */
static int cb_check_integrity_mult(const corax_utree_t *tree,
                                   const corax_unode_t *node)
{
  unsigned int clv_index     = node->clv_index;
  int          scaler_index  = node->scaler_index;
  unsigned int pmatrix_index = node->pmatrix_index;
  char *       label         = node->label;
  double       length        = node->length;
  unsigned int subnodes      = 1;

  /* edge attributes */
  if (node->back->length != length)
  {
    corax_set_error(0,
                    "Inconsistent branch lengths: %lf != %lf",
                    length,
                    node->back->length);
    return CORAX_FAILURE;
  }

  if (node->back->pmatrix_index != pmatrix_index)
  {
    corax_set_error(0,
                    "Inconsistent pmatrix indices: %u != %u",
                    pmatrix_index,
                    node->back->pmatrix_index);
    return CORAX_FAILURE;
  }

  if (node->next)
  {
    /* node attributes */
    corax_unode_t *snode = node->next;
    do {
      subnodes++;

      if (tree->binary && subnodes > 3)
      {
        corax_set_error(0,
                        "Multifurcation found in a binary tree "
                        "at node with clv_index = %u",
                        snode->clv_index);
        return CORAX_FAILURE;
      }

      if (subnodes > tree->tip_count)
      {
        corax_set_error(0,
                        "Multifurcation exceeding the tree size found "
                        "at node with clv_index = %u",
                        snode->clv_index);
        return CORAX_FAILURE;
      }

      if (snode->clv_index != clv_index)
      {
        corax_set_error(0,
                        "Inconsistent CLV indices: %u != %u",
                        clv_index,
                        snode->clv_index);
        return CORAX_FAILURE;
      }
      if (snode->scaler_index != scaler_index)
      {
        corax_set_error(0,
                        "Inconsistent scaler indices: %d != %d",
                        scaler_index,
                        snode->scaler_index);
        return CORAX_FAILURE;
      }
      if (snode->label != label)
      {
        corax_set_error(
            0, "Inconsistent node labels: '%s' != '%s'", label, snode->label);
        return CORAX_FAILURE;
      }
      if (!snode->next)
      {
        corax_set_error(0,
                        "Open roundabout (node->next is NULL) "
                        "at node with clv_index = %u",
                        snode->clv_index);
        return CORAX_FAILURE;
      }
      snode = snode->next;
    } while (snode != node);
  }

  return 1;
}

CORAX_EXPORT int corax_utree_check_integrity(const corax_utree_t *tree)
{
  return corax_utree_every_const(tree, cb_check_integrity_mult);
}

/* TODO: Memory allocation checks were not implemented in this function!!! */
static corax_unode_t *clone_node(const corax_unode_t *node)
{
  corax_unode_t *new_node = (corax_unode_t *)malloc(sizeof(corax_unode_t));
  memcpy(new_node, node, sizeof(corax_unode_t));

  if (node->label)
  {
    new_node->label = (char *)malloc(strlen(node->label) + 1);
    strcpy(new_node->label, node->label);
  }

  if (node->next)
  {
    corax_unode_t *snode     = node->next;
    corax_unode_t *new_snode = new_node;
    do {
      new_snode->next = (corax_unode_t *)malloc(sizeof(corax_unode_t));
      memcpy(new_snode->next, snode, sizeof(corax_unode_t));
      new_snode->next->label = new_node->label;
      snode                  = snode->next;
      new_snode              = new_snode->next;
    } while (snode != node);

    new_snode->next = new_node;
  }

  return new_node;
}

static void utree_recurse_clone(corax_unode_t *      new_root,
                                const corax_unode_t *root)
{
  const corax_unode_t *node = root->back;
  if (node)
  {
    new_root->back       = clone_node(node);
    new_root->back->back = new_root;

    if (node->next)
    {
      corax_unode_t *snode     = node->next;
      corax_unode_t *new_snode = new_root->back->next;
      do {
        utree_recurse_clone(new_snode, snode);
        snode     = snode->next;
        new_snode = new_snode->next;
      } while (snode && snode != node);
    }
  }
}

CORAX_EXPORT corax_unode_t *corax_utree_graph_clone(const corax_unode_t *root)
{
  corax_unode_t *new_root = clone_node(root);

  const corax_unode_t *snode     = root;
  corax_unode_t *      new_snode = new_root;
  do {
    utree_recurse_clone(new_snode, snode);
    snode     = snode->next;
    new_snode = new_snode->next;
  } while (snode && snode != root);

  return new_root;
}

CORAX_EXPORT corax_utree_t *corax_utree_clone(const corax_utree_t *tree)
{
  /* choose the last inner node as the starting point of the clone. It does not
    really matter which node to choose, but since the newick parser places the
    root node at the end of the list, we use the same notation here */
  corax_unode_t *root = corax_utree_graph_clone(tree->vroot);

  if (tree->binary)
    return corax_utree_wraptree(root, tree->tip_count);
  else
    return corax_utree_wraptree_multi(root, tree->tip_count, tree->inner_count);
}

CORAX_EXPORT void corax_utree_graph_destroy(corax_unode_t *root,
                                            void (*cb_destroy)(void *))
{
  if (!root) return;

  dealloc_graph_recursive(root, cb_destroy, 0);
}

CORAX_EXPORT void corax_utree_destroy(corax_utree_t *tree,
                                      void (*cb_destroy)(void *))
{
  unsigned int i;

  /* deallocate tip nodes */
  for (i = 0; i < tree->tip_count; ++i)
  {
    dealloc_data(tree->nodes[i], cb_destroy);
    if (tree->nodes[i]->label) free(tree->nodes[i]->label);
    free(tree->nodes[i]);
  }

  /* deallocate inner nodes */
  for (i = tree->tip_count; i < tree->tip_count + tree->inner_count; ++i)
  {
    corax_unode_t *first = tree->nodes[i];

    assert(first);

    if (first->label) free(first->label);

    corax_unode_t *node = first;
    do {
      corax_unode_t *next = node->next;
      dealloc_data(node, cb_destroy);
      free(node);
      node = next;
    } while (node && node != first);
  }

  /* deallocate tree structure */
  free(tree->nodes);
  free(tree);
}

static void recursive_assign_indices(corax_unode_t *node,
                                     unsigned int * tip_clv_index,
                                     unsigned int * inner_clv_index,
                                     int *          inner_scaler_index,
                                     unsigned int * inner_node_index,
                                     unsigned int   level)
{
  if (!node->next)
  {
    /* tip node */
    node->node_index    = *tip_clv_index;
    node->clv_index     = *tip_clv_index;
    node->pmatrix_index = *tip_clv_index;
    node->scaler_index  = CORAX_SCALE_BUFFER_NONE;
    *tip_clv_index      = *tip_clv_index + 1;
  }
  else
  {
    /* inner node */
    corax_unode_t *snode = level ? node->next : node;
    do {
      recursive_assign_indices(snode->back,
                               tip_clv_index,
                               inner_clv_index,
                               inner_scaler_index,
                               inner_node_index,
                               level + 1);
      snode = snode->next;
    } while (snode != node);

    snode = node;
    do {
      snode->node_index   = (*inner_node_index)++;
      snode->clv_index    = *inner_clv_index;
      snode->scaler_index = *inner_scaler_index;
      if (snode == node && level > 0)
        snode->pmatrix_index = *inner_clv_index;
      else
        snode->pmatrix_index = snode->back->pmatrix_index;
      snode = snode->next;
    } while (snode != node);

    *inner_clv_index += 1;
    *inner_scaler_index += 1;
  }
}

CORAX_EXPORT void corax_utree_reset_template_indices(corax_unode_t *root,
                                                     unsigned int   tip_count)
{
  unsigned int tip_clv_index      = 0;
  unsigned int inner_clv_index    = tip_count;
  unsigned int inner_node_index   = tip_count;
  int          inner_scaler_index = 0;

  if (!root->next) root = root->back;

  recursive_assign_indices(root,
                           &tip_clv_index,
                           &inner_clv_index,
                           &inner_scaler_index,
                           &inner_node_index,
                           0);
}

static void fill_nodes_recursive(corax_unode_t * node,
                                 corax_unode_t **array,
                                 unsigned int    array_size,
                                 unsigned int *  tip_index,
                                 unsigned int *  inner_index,
                                 unsigned int    level)
{
  unsigned int index;
  if (!node->next)
  {
    /* tip node */
    index = *tip_index;
    *tip_index += 1;
  }
  else
  {
    /* inner node */
    corax_unode_t *snode = level ? node->next : node;
    do {
      fill_nodes_recursive(
          snode->back, array, array_size, tip_index, inner_index, level + 1);
      snode = snode->next;
    } while (snode != node);

    index = *inner_index;
    *inner_index += 1;
  }

  assert(index < array_size);
  array[index] = node;
}

static unsigned int utree_count_nodes_recursive(corax_unode_t *node,
                                                unsigned int * tip_count,
                                                unsigned int * inner_count,
                                                unsigned int   level)
{
  if (!node->next)
  {
    *tip_count += 1;
    return 1;
  }
  else
  {
    unsigned int count = 0;

    corax_unode_t *snode = level ? node->next : node;
    do {
      count += utree_count_nodes_recursive(
          snode->back, tip_count, inner_count, level + 1);
      snode = snode->next;
    } while (snode != node);

    *inner_count += 1;

    return count + 1;
  }
}

static unsigned int utree_count_nodes(corax_unode_t *root,
                                      unsigned int * tip_count,
                                      unsigned int * inner_count)
{
  unsigned int count = 0;

  if (tip_count) *tip_count = 0;

  if (inner_count) *inner_count = 0;

  if (!root->next && !root->back->next) return 0;

  if (!root->next) root = root->back;

  count = utree_count_nodes_recursive(root, tip_count, inner_count, 0);

  if (tip_count && inner_count) assert(count == *tip_count + *inner_count);

  return count;
}

CORAX_EXPORT int corax_unode_is_rooted(const corax_unode_t *root)
{
  return (root->next && root->next->next == root) ? 1 : 0;
}

CORAX_EXPORT int corax_utree_is_rooted(const corax_utree_t *tree)
{
  return corax_unode_is_rooted(tree->vroot);
}

static corax_utree_t *utree_wraptree(corax_unode_t *root,
                                     unsigned int   tip_count,
                                     unsigned int   inner_count,
                                     int            binary)
{
  unsigned int node_count;

  corax_utree_t *tree = (corax_utree_t *)malloc(sizeof(corax_utree_t));
  if (!tree)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC, "Unable to allocate enough memory.");
    return CORAX_FAILURE;
  }

  if (tip_count < 3 && tip_count != 0)
  {
    corax_set_error(
        CORAX_ERROR_INVALID_PARAM, "Invalid tip_count value (%u).", tip_count);
    return CORAX_FAILURE;
  }

  if (!root->next) root = root->back;

  if (binary)
  {
    if (tip_count == 0)
    {
      node_count = utree_count_nodes(root, &tip_count, &inner_count);
      if (inner_count != tip_count - 2)
      {
        corax_set_error(CORAX_ERROR_INVALID_PARAM,
                        "Input tree is not strictly bifurcating.");
        return CORAX_FAILURE;
      }
    }
    else
    {
      inner_count = tip_count - 2;
      node_count  = tip_count + inner_count;
    }
  }
  else
  {
    if (tip_count == 0 || inner_count == 0)
      node_count = utree_count_nodes(root, &tip_count, &inner_count);
    else
      node_count = tip_count + inner_count;
  }

  if (!tip_count)
  {
    corax_set_error(CORAX_ERROR_INVALID_PARAM,
                    "Input tree contains no inner nodes.");
    return CORAX_FAILURE;
  }

  tree->nodes = (corax_unode_t **)malloc(node_count * sizeof(corax_unode_t *));
  if (!tree->nodes)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC, "Unable to allocate enough memory.");
    return CORAX_FAILURE;
  }

  unsigned int tip_index   = 0;
  unsigned int inner_index = tip_count;

  fill_nodes_recursive(
      root, tree->nodes, node_count, &tip_index, &inner_index, 0);

  assert(tip_index == tip_count);
  assert(inner_index == tip_count + inner_count);

  tree->tip_count   = tip_count;
  tree->inner_count = inner_count;
  tree->edge_count  = node_count - 1;
  tree->binary =
      (inner_count == tip_count - (corax_unode_is_rooted(root) ? 1 : 2));
  tree->vroot = root;

  return tree;
}

/* wraps/encalupsates the unrooted tree graph into a tree structure
   that contains a list of nodes, number of tips and number of inner
   nodes. If 0 is passed as tip_count, then an additional recrursion
   of the tree structure is done to detect the number of tips */
CORAX_EXPORT corax_utree_t *corax_utree_wraptree(corax_unode_t *root,
                                                 unsigned int   tip_count)
{
  return utree_wraptree(root, tip_count, 0, 1);
}

CORAX_EXPORT corax_utree_t *corax_utree_wraptree_multi(corax_unode_t *root,
                                                       unsigned int   tip_count,
                                                       unsigned int inner_count)
{
  return utree_wraptree(root, tip_count, inner_count, 0);
}

/**
 * @brief Creates a new circular node
 *
 *           n2
 *          / |
 *        n1  |
 *          \ |
 *           n3
 *
 * All parameters are shared among the nodes in the triplet
 *
 * @param clv_index    the clv_index
 * @param scaler_index the scaler index
 * @param label        the node label
 * @param data         the data pointer
 *
 * @return the new node
 */
CORAX_EXPORT corax_unode_t *corax_utree_create_node(unsigned int clv_index,
                                                    int          scaler_index,
                                                    char *       label,
                                                    void *       data)
{
  corax_unode_t *new_node = (corax_unode_t *)calloc(1, sizeof(corax_unode_t));
  new_node->next          = (corax_unode_t *)calloc(1, sizeof(corax_unode_t));
  new_node->next->next    = (corax_unode_t *)calloc(1, sizeof(corax_unode_t));
  if (!(new_node && new_node->next && new_node->next->next))
  {
    if (new_node)
    {
      if (new_node->next)
      {
        free(new_node->next->next);
        free(new_node->next);
      }
      free(new_node);
    }
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Cannot allocate memory for new node\n");
    return NULL;
  }

  new_node->next->next->next = new_node;
  new_node->label            = label;
  new_node->next->label = new_node->next->next->label = new_node->label;
  new_node->next->data = new_node->next->next->data = new_node->data = data;
  new_node->next->length = new_node->next->next->length = new_node->length = 0;
  new_node->next->clv_index    = new_node->next->next->clv_index =
      new_node->clv_index      = clv_index;
  new_node->next->scaler_index = new_node->next->next->scaler_index =
      new_node->scaler_index   = scaler_index;
  new_node->back = new_node->next->back = new_node->next->next->back = NULL;
  return new_node;
}

struct clv_set_data
{
  int *        set_indices;
  unsigned int max_index;
  unsigned int tip_count;
};

static int cb_set_clv_minimal(corax_unode_t *node, void *data)
{
  unsigned int         i, next_index;
  int                  index_found;
  struct clv_set_data *clv_data = (struct clv_set_data *)data;
  int *                v        = 0;

  if (!CORAX_UTREE_IS_TIP(node))
  {
    /* find next free position */
    v           = clv_data->set_indices;
    next_index  = 0;
    index_found = 0;
    for (i = 0; i < clv_data->max_index; ++i)
    {
      if (!v[i])
      {
        index_found = 1;
        next_index  = i;
        v[i]        = 1;
        break;
      }
    }
    assert(index_found);

    /* set clv index */
    node->clv_index = node->next->clv_index = node->next->next->clv_index =
        next_index + clv_data->tip_count;
    /* set scaler index */
    node->scaler_index = node->next->scaler_index =
        node->next->next->scaler_index =
            (int)(next_index + clv_data->tip_count);

    /* free indices from children */
    if (!CORAX_UTREE_IS_TIP(node->next->back))
    {
      v[node->next->back->clv_index - clv_data->tip_count] = 0;
    }
    if (!CORAX_UTREE_IS_TIP(node->next->next->back))
    {
      v[node->next->next->back->clv_index - clv_data->tip_count] = 0;
    }
  }

  /* continue */
  return 1;
}

CORAX_EXPORT int corax_utree_set_clv_minimal(corax_unode_t *root,
                                             unsigned int   tip_count)
{
  unsigned int clv_count   = (unsigned int)ceil(log2(tip_count)) + 2;
  int *        set_indices = (int *)calloc((size_t)clv_count, sizeof(int));
  struct clv_set_data data;
  data.set_indices = set_indices;
  data.max_index   = clv_count;
  data.tip_count   = tip_count;
  corax_utree_traverse_apply(root, 0, 0, cb_set_clv_minimal, (void *)&data);
  free(set_indices);

  return CORAX_SUCCESS;
}
