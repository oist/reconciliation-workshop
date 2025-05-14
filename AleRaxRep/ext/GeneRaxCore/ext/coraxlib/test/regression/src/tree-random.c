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

#include <string.h>
#include <stdarg.h>
#include <assert.h>
#include <stdlib.h>
#include <time.h>

#define ALPHA 0.5
#define N_STATES 4
#define N_CAT_GAMMA 4

#define N_TAXA_BIG 43
#define N_TAXA_SMALL 5
#define N_SITES 491
#define MAX_TAXA 200

int main (int argc, char * argv[])
{
   unsigned int i, n_taxa;
   char * seq, *header[MAX_TAXA];
   long seq_len, header_len, seqno;
   long n_sites = 0;
   corax_fasta_t * fp;
   unsigned int params_indices[4] = {0, 0, 0, 0};
   unsigned int attributes = get_attributes(argc, argv);

   fp = corax_fasta_open ("testdata/small.fas", corax_map_fasta);
   if (!fp)
   {
     printf (" ERROR opening file (%d): %s\n", corax_errno, corax_errmsg);
     exit (CORAX_FAILURE);
   }

   /* first read for getting number of taxa and headers */
   i = 0;
   while (corax_fasta_getnext (fp, &header[i], &header_len, &seq, &seq_len, &seqno))
   {
     if (!n_sites)
       n_sites = seq_len;
     else if (seq_len != n_sites)
     {
       printf (
           " ERROR: Mismatching sequence length for sequence %d (%ld, and it should be %d)\n",
           i, seq_len - 1, N_SITES);
       exit (CORAX_FAILURE);
     }

     printf ("Header of sequence %d(%ld) %s (%ld sites)\n", i, seqno, header[i],
             seq_len);
     free (seq);
     ++i;
   }
 //  corax_fasta_close (fp);
   n_taxa = i;

   if (corax_errno != CORAX_ERROR_FILE_EOF)
   {
     printf (" ERROR at the end (%d): %s\n", corax_errno, corax_errmsg);
     exit (CORAX_FAILURE);
   }

   corax_utree_t * random_tree = corax_utree_random_create(n_taxa,
                                                          (const char **)header,
                                                          42);
   corax_unode_t * tree = random_tree->nodes[2*n_taxa - 3];

   for (i=0; i<n_taxa; ++i)
     free(header[i]);

   corax_utree_show_ascii(tree, CORAX_UTREE_SHOW_CLV_INDEX | CORAX_UTREE_SHOW_LABEL | CORAX_UTREE_SHOW_PMATRIX_INDEX);

    unsigned int tip_nodes_count, inner_nodes_count, nodes_count, branch_count;
    unsigned int matrix_count, ops_count;
    unsigned int * matrix_indices;
    double * branch_lengths;
    corax_partition_t * partition;
    corax_operation_t * operations;
    corax_unode_t ** travbuffer;

    tip_nodes_count = n_taxa;
    inner_nodes_count = n_taxa - 2;
    nodes_count = tip_nodes_count + inner_nodes_count;
    branch_count = 2 * tip_nodes_count - 3;


    printf("SCALE BUFFERS = %d\n", inner_nodes_count);
  /* create the PLL partition instance */
    partition = corax_partition_create (tip_nodes_count,
                                      inner_nodes_count,
                                      N_STATES,
                                      (unsigned int) n_sites,
                                      1,
                                      branch_count,
                                      N_CAT_GAMMA,
                                      inner_nodes_count,
                                      attributes
                                      );

  /* second read for assigning sequences */
    corax_fasta_rewind(fp);
    i = 0;
    char *tmp_header;
    while (corax_fasta_getnext (fp, &tmp_header, &header_len, &seq, &seq_len, &seqno))
    {
      if (seq_len != n_sites)
      {
        printf (
            " ERROR: Mismatching sequence length for sequence %d (%ld, and it should be %d)\n",
            i, seq_len - 1, N_SITES);
        exit (CORAX_FAILURE);
      }
      corax_set_tip_states (partition, i, corax_map_nt, seq);
      free (tmp_header);
      free (seq);
      ++i;
    }
    corax_fasta_close (fp);

    /* initialize base frequencies */
    double frequencies[4] = {0.25,0.25,0.25,0.25};
    corax_set_frequencies (partition, 0, frequencies);

    /* initialize substitution rates */
    double subst_params[6] = {1,1,1,1,1,1};
    corax_set_subst_params (partition, 0,  subst_params);

    /* compute the discretized category rates from a gamma distribution
       with alpha shape 1 and store them in rate_cats  */
    double rate_cats[N_CAT_GAMMA] = { 0 };
    corax_compute_gamma_cats (1, N_CAT_GAMMA, rate_cats, CORAX_GAMMA_RATES_MEAN);
    corax_set_category_rates (partition, rate_cats);

    travbuffer = (corax_unode_t **) malloc (nodes_count * sizeof(corax_unode_t *));

    branch_lengths = (double *) malloc (branch_count * sizeof(double));
    matrix_indices = (unsigned int *) malloc (
        branch_count * sizeof(unsigned int));
    operations = (corax_operation_t *) malloc (
        inner_nodes_count * sizeof(corax_operation_t));

    /* perform a postorder traversal of the unrooted tree */
    unsigned int traversal_size;
    if (!corax_utree_traverse (tree,
                             CORAX_TREE_TRAVERSE_POSTORDER,
                             cb_full_traversal,
                             travbuffer,
                             &traversal_size))
      fatal ("Function corax_utree_traverse() requires inner nodes as parameters");

    printf("\nTRAVBUFFER: ");
       for (i=0; i<nodes_count; ++i)
         printf("%d/%d/%d  XX ", travbuffer[i]->clv_index, travbuffer[i]->scaler_index, travbuffer[i]->pmatrix_index);
       printf("\n");

    corax_utree_create_operations (travbuffer, traversal_size, branch_lengths,
                                   matrix_indices, operations, &matrix_count,
                                   &ops_count);

    printf("\nMATRICES: ");
    for (i=0; i<matrix_count; ++i)
      printf("%d XX ", matrix_indices[i]);
    printf("\n");
    printf("\nSCALERS: ");
    for (i=0; i<ops_count; ++i)
      printf("%d/%d XX %d/%d XX %d/%d\n", operations[i].child1_clv_index, operations[i].child1_scaler_index, operations[i].child2_clv_index,
             operations[i].child2_scaler_index, operations[i].parent_clv_index, operations[i].parent_scaler_index);
    printf("\n");

    printf ("Traversal size: %d\n", traversal_size);
    printf ("Operations: %d\n", ops_count);
    printf ("Probability Matrices: %d\n", matrix_count);

    corax_update_prob_matrices (partition,
                              params_indices,
                              matrix_indices,
                              branch_lengths,
                              matrix_count);

    corax_update_clvs (partition, operations, ops_count);

    double logl = corax_compute_edge_loglikelihood (partition,
                                                  tree->clv_index,
                                                  tree->scaler_index,
                                                  tree->back->clv_index,
                                                  tree->back->scaler_index,
                                                  tree->pmatrix_index,
                                                  params_indices,
                                                  NULL);

   printf ("Log-L at %s-%s: %f\n", tree->label, tree->back->label, logl);

   free (travbuffer);
   free (branch_lengths);
   free (matrix_indices);
   free (operations);

   corax_partition_destroy (partition);
   corax_utree_destroy(random_tree, NULL);
   return CORAX_SUCCESS;
}
