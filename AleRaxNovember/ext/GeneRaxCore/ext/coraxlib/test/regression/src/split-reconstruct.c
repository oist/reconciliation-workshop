/*
 Copyright (C) 2016 Diego Darriba

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

 Contact: Diego Darriba <Diego.Darriba@h-its.org>,
 Exelixis Lab, Heidelberg Instutute for Theoretical Studies
 Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
 */
#include "rng.h"
#include "common.h"

#include <assert.h>

#define TREEFILE  "testdata/medium.tree"

/*
 * This example reads an input tree, generates the set of splits, shuffles them
 * randomly and reconstructs the tree out of the splits
 */
/* static functions */
static void shuffle(corax_split_t *array, size_t n);
static void print_newick_recurse(corax_unode_t * node);
static void print_newick(corax_unode_t * tree);

static const unsigned int n_iters = 10;
int main (int argc, char * argv[])
{
  unsigned int rf_dist, i, iter;
  char **labels;

  /* tree properties */
  corax_unode_t * tree = NULL;
  unsigned int tip_count;
  unsigned int attributes = get_attributes(argc, argv);

  if (attributes != CORAX_ATTRIB_ARCH_CPU)
  {
    skip_test();
  }

  /* parse the input trees */
  corax_utree_t * parsed_tree = corax_utree_parse_newick (TREEFILE);
  tip_count = parsed_tree->tip_count;
  tree = parsed_tree->nodes[2*tip_count - 3];
  if (!tree)
  {
    fatal("Error %d: %s", corax_errno, corax_errmsg);
  }

  corax_unode_t ** tipnodes = parsed_tree->nodes;

  labels = (char **) malloc(tip_count * sizeof(char *));
  for (i=0; i<tip_count; ++i)
    labels[tipnodes[i]->node_index] = tipnodes[i]->label;

  unsigned int n_splits = tip_count - 3;
  corax_split_t * splits = corax_utree_split_create(tree,
                                                   tip_count,
                                                   NULL);

  for (iter=0; iter<n_iters; ++iter)
  {
    shuffle(splits, n_splits);

    corax_split_system_t split_system;
    split_system.splits = splits;
    split_system.support = 0;
    split_system.split_count = n_splits;
    split_system.max_support = 1.0;

    corax_consensus_utree_t * constree = corax_utree_from_splits(&split_system,
                                                                tip_count,
                                                                labels);

    corax_utree_t * consensus = corax_utree_wraptree(constree->tree, tip_count);
    if (!corax_utree_consistency_set(consensus, parsed_tree))
       fatal("Cannot set trees consistent!");

    corax_utree_show_ascii(constree->tree, CORAX_UTREE_SHOW_CLV_INDEX | CORAX_UTREE_SHOW_LABEL);
    print_newick(constree->tree);

    corax_split_t * splits2 = corax_utree_split_create(constree->tree,
                                                      tip_count,
                                                      NULL);

    corax_utree_split_normalize_and_sort(splits2,
                                          tip_count,
                                          n_splits,
                                          1);

    /* sort splits back */
    corax_utree_split_normalize_and_sort(splits,
                                          tip_count,
                                          n_splits,
                                          0);

    rf_dist = corax_utree_split_rf_distance(splits, splits2, tip_count);
    printf(" RF DIST = %d\n", rf_dist);

    if (rf_dist > 0)
      fatal("Error: Initial and reconstructed trees differ!");

    /* in-loop cleanup */
    corax_utree_consensus_destroy(constree);
    free (consensus->nodes);
    free (consensus);
    corax_utree_split_destroy(splits2);
  }

  /* clean */
  free(labels);
  corax_utree_destroy (parsed_tree, NULL);
  corax_utree_split_destroy(splits);

  return (0);
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

static void shuffle(corax_split_t *array, size_t n)
{
  unsigned int i, j;
  corax_split_t t;
  for (i = 0; i < n - 1; i++)
  {
    j = i + (unsigned int) RAND / (RAND_MAX / (n - i) + 1);
    t = array[j];
    array[j] = array[i];
    array[i] = t;
  }
}

static void print_newick_recurse(corax_unode_t * node)
{
  corax_unode_t * child;
  if (CORAX_UTREE_IS_TIP(node))
  {
    printf("%s", node->label);
    return;
  }

  printf("(");
  child = node->next;
  while(child != node)
  {
    print_newick_recurse(child->back);

    if (child->next != node)
      printf(",");

    child = child->next;
  }
  printf(")");
}

static void print_newick(corax_unode_t * tree)
{
  printf("(");
  print_newick_recurse(tree->back);
  corax_unode_t * child = tree->next;
  while(child != tree)
  {
    printf(",");
    print_newick_recurse(child->back);
    child = child->next;
  }
  printf(");\n");
}
