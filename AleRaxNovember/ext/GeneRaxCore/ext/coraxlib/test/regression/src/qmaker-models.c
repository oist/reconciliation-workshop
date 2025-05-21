/*
    Copyright (C) 2015 Diego Darriba, Tomas Flouri

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
#include "common.h"

#define FASTAFILE "testdata/ribosomal_l5_pf00673.fas"
#define TREEFILE "testdata/ribosomal_l5_pf00673.tree"

#define N_RATE_CATS 4
#define ALPHA 1

#define N_PROT_MODELS 8
#define N_STATES 20

#define FLOAT_PRECISION 5

static const char *prot_model_names[N_PROT_MODELS] = {
    "Q.PFAM", "Q.PFAM_GB", "Q.LG",  "Q.BIRD",
    "Q.INSECT", "Q.MAMMAL", "Q.PLANT", "Q.YEAST"};

int main(int argc, char *argv[])
{
  unsigned int i, cur_model;
  corax_utree_t * tree;
  corax_unode_t * root;
  unsigned int    taxa_count, nodes_count, inner_nodes_count, branch_count;
  int             retval;

  unsigned int     traversal_size, matrix_count, ops_count;
  corax_unode_t **   travbuffer;
  unsigned int *   matrix_indices;
  double *         branch_lengths;
  double           logl;

  corax_partition_t *partition;
  corax_operation_t *operations;
  unsigned int     params_indices[N_RATE_CATS] = {0, 0, 0, 0};

  unsigned int attributes = get_attributes(argc, argv);

  tree = corax_utree_parse_newick (TREEFILE);

  taxa_count        = tree->tip_count;
  root              = tree->vroot;
  inner_nodes_count = tree->inner_count;
  nodes_count       = taxa_count + inner_nodes_count;
  branch_count      = tree->edge_count;

  printf("Creating PLL partition\n");

  partition = parse_msa(FASTAFILE, N_STATES, N_RATE_CATS, 1, tree, attributes);


  double *rate_cats = (double *)malloc(N_RATE_CATS * sizeof(double));

  if (N_RATE_CATS > 1)
  {
    if (corax_compute_gamma_cats(
            ALPHA, N_RATE_CATS, rate_cats, CORAX_GAMMA_RATES_MEAN)
        == CORAX_FAILURE)
    {
      printf("Fail computing the gamma rates\n");
      exit(1);
    }
    corax_set_category_rates(partition, rate_cats);
    free(rate_cats);
  }

  /* build fixed structures */
  travbuffer     = (corax_unode_t **)malloc(nodes_count * sizeof(corax_unode_t *));
  branch_lengths = (double *)malloc(branch_count * sizeof(double));
  matrix_indices = (unsigned int *)malloc(branch_count * sizeof(unsigned int));
  operations =
      (corax_operation_t *)malloc(inner_nodes_count * sizeof(corax_operation_t));

  retval = corax_utree_traverse(root,
                              CORAX_TREE_TRAVERSE_POSTORDER,
                              cb_full_traversal,
                              travbuffer,
                              &traversal_size);

  if (!retval)
  {
    printf("ERROR: corax_utree traversal failed: %s\n", corax_errmsg);
    exit(-1);
  }

  corax_utree_create_operations(travbuffer,
                              traversal_size,
                              branch_lengths,
                              matrix_indices,
                              operations,
                              &matrix_count,
                              &ops_count);


  for (cur_model = 0; cur_model < N_PROT_MODELS; cur_model++)
  {

    printf("\nSetting model %s...\n", prot_model_names[cur_model]);

    corax_subst_model_t * modinfo = corax_util_model_info_protein(prot_model_names[cur_model]);

    if (!modinfo)
    {
      printf("ERROR: corax_util_model_info_protein failed: %s\n", corax_errmsg);
      exit(-1);
    }

    corax_set_subst_params(partition, 0, modinfo->rates);
    corax_set_frequencies(partition, 0, modinfo->freqs);

    printf("Updating prob matrices...\n");

    corax_update_prob_matrices(
        partition, params_indices, matrix_indices, branch_lengths, matrix_count);

    // corax_update_invariant_sites_proportion(partition, 0.17);
    for (i = 0; i < 1; ++i)
    {
      printf("P-matrix for branch length %f\n", branch_lengths[i]);
      corax_show_pmatrix(partition, i, FLOAT_PRECISION);
      printf("\n");
    }

    corax_update_clvs(partition, operations, ops_count);

    logl = corax_compute_edge_loglikelihood(partition,
                                          root->clv_index,
                                          root->scaler_index,
                                          root->back->clv_index,
                                          root->back->scaler_index,
                                          root->pmatrix_index,
                                          params_indices,
                                          NULL);

    printf("Log-L (%s): %.6f\n", prot_model_names[cur_model], logl);

    corax_util_model_destroy(modinfo);
  }

  /* clean */
  free(travbuffer);
  free(branch_lengths);
  free(operations);
  free(matrix_indices);
  corax_partition_destroy(partition);

  corax_utree_destroy(tree, NULL);

  return (0);
}
