/*
 Copyright (C) 2017 Diego Darriba, Pierre Barbera

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
#include <search.h>

#define ALPHA        0.5
#define N_STATES      4
#define N_SUBST_RATES 6
#define N_RATE_CATS   4

#define MSA_FILENAME  "testdata/246x4465.fas"
#define TREE_FILENAME "testdata/246x4465.tree"

#define BLOCK_ID_PARTITION -1
#define BLOCK_ID_TREE      -2

unsigned int get_sites_number_scaler(const corax_partition_t * partition,
                                                unsigned int scaler_index)
{
    unsigned int sites = partition->attributes & CORAX_ATTRIB_SITE_REPEATS ?
            partition->repeats->perscale_ids[scaler_index] : 0;
      sites = sites ? sites : partition->sites;
      return sites;
}

static corax_partition_t * parse_msa(unsigned int attributes, corax_utree_t * tree)
{
  unsigned int i;
  unsigned int tip_clv_index;
  char * seq, *header;
  char ** headers, ** seqdata;
  long seq_len    = 0,
       read_len   = 0,
       header_len = 0,
       seqno      = 0;
  corax_fasta_t * fp;
  corax_partition_t * partition;

  corax_unode_t ** tipnodes = tree->nodes;
  unsigned int tip_nodes_count = tree->tip_count;

  headers = (char **)calloc(tip_nodes_count, sizeof(char *));
  seqdata = (char **)calloc(tip_nodes_count, sizeof(char *));

  /* create a libc hash table of size tip_nodes_count */
  hcreate(tip_nodes_count);

  /* populate a libc hash table with tree tip labels */
  unsigned int * data = (unsigned int *) malloc ( tip_nodes_count *
                                                  sizeof(unsigned int) );
  for (i = 0; i < tip_nodes_count; ++i)
  {
    data[i] = i;
    ENTRY entry;
    entry.key = tipnodes[i]->label;
    entry.data = (void *)(data+i);
    hsearch(entry, ENTER);
  }

  fp = corax_fasta_open (MSA_FILENAME, corax_map_fasta);
  if (!fp)
  {
    printf (" ERROR opening file (%d): %s\n", corax_errno, corax_errmsg);
    return NULL;
  }

  i = 0;
  while (corax_fasta_getnext (fp, &header, &header_len, &seq, &read_len, &seqno))
  {
    if (!seq_len)
    {
      seq_len = read_len;
    }
    else if (seq_len != read_len)
    {
      printf (
          " ERROR: Mismatching sequence length for sequence %d (%ld, and it should be %ld)\n",
          i, read_len - 1, seq_len);
      return NULL;
    }
    headers[i] = header;
    seqdata[i] = seq;
    ++i;
  }

  if (corax_errno != CORAX_ERROR_FILE_EOF)
  {
    printf (" ERROR at the end (%d): %s\n", corax_errno, corax_errmsg);
    return NULL;
  }

  corax_fasta_close (fp);

  partition = corax_partition_create(tip_nodes_count,      /* tips */
                                    tip_nodes_count - 2, /* clv buffers */
                                    N_STATES,            /* states */
                                    seq_len,             /* sites */
                                    1,           /* different rate parameters */
                                    2*tip_nodes_count - 3, /* prob matrices */
                                    N_RATE_CATS,           /* rate categories */
                                    tip_nodes_count - 2,   /* scale buffers */
                                    attributes             /* attributes */
                                    );

  for (i = 0; i < tip_nodes_count; ++i)
  {
    ENTRY query;
    query.key = headers[i];
    ENTRY * found = NULL;

    found = hsearch(query,FIND);

    if (!found)
      fatal("Sequence with header %s does not appear in the tree", headers[i]);

    tip_clv_index = *((unsigned int *)(found->data));
    corax_set_tip_states(partition, tip_clv_index, corax_map_nt, seqdata[i]);
  }

  hdestroy();
  free(data);

  for(i = 0; i < tip_nodes_count; ++i)
  {
    free(seqdata[i]);
    free(headers[i]);
  }
  free(seqdata);
  free(headers);

  return partition;
}

typedef struct
{
  int * a;
  int * b;
} test_t;

int write(void * data, size_t size, size_t count, FILE * file)
{
  size_t ret = fwrite(data, size, count, file);
  if (ret != count)
  {
    return CORAX_FAILURE;
  }
  return CORAX_SUCCESS;
}

static long int get_offset(const size_t n_blocks, corax_block_map_t * block_map, int id)
{
  for (size_t i = 0; i<n_blocks; ++i)
  {
    if (block_map[i].block_id == id)
    {
      return block_map[i].block_offset;
    }
  }
  printf("no such block ID in map!\n");
  return -1;
}

int main (int argc, char * argv[])
{
  unsigned int attributes = get_attributes(argc, argv);
  unsigned int matrix_count, ops_count;
  unsigned int tip_nodes_count, inner_nodes_count, nodes_count, branch_count;
  corax_partition_t * old_partition;
  corax_partition_t * partition;
  corax_unode_t * tree;
  double logl, save_logl;

  corax_unode_t ** travbuffer;
  double * branch_lengths;
  unsigned int * matrix_indices;
  corax_operation_t * operations;

  double frequencies[N_STATES]             = { 0.17, 0.19, 0.25, 0.39 };
  double subst_params[N_SUBST_RATES]       = {1,1,1,1,1,1};
  double rate_cats[N_RATE_CATS]            = {0};
  unsigned int params_indices[N_RATE_CATS] = {0, 0, 0, 0};

  unsigned int lk_parent_clv_index, lk_parent_scaler_index,
               lk_child_clv_index, lk_child_scaler_index,
               lk_pmatrix_index;

  corax_utree_t * parsed_tree = corax_utree_parse_newick(TREE_FILENAME);
  tip_nodes_count = parsed_tree->tip_count;
  tree = parsed_tree->nodes[2*tip_nodes_count - 3];

  if (!tree)
  {
    printf ("Error %d parsing tree: %s\n", corax_errno, corax_errmsg);
    return 1;
  }

  inner_nodes_count = tip_nodes_count - 2;
  nodes_count = inner_nodes_count + tip_nodes_count;
  branch_count = nodes_count - 1;

  printf("Number of tip/leaf nodes in tree: %d\n", tip_nodes_count);
  printf("Number of inner nodes in tree: %d\n", inner_nodes_count);
  printf("Total number of nodes in tree: %d\n", nodes_count);
  printf("Number of branches in tree: %d\n", branch_count);

  old_partition = parse_msa (attributes, parsed_tree);
  if (!old_partition)
  {
    printf ("Error creating old_partition\n");
    return 1;
  }

  corax_compute_gamma_cats(ALPHA, N_RATE_CATS, rate_cats, CORAX_GAMMA_RATES_MEAN);
  corax_set_frequencies(old_partition, 0, frequencies);
  corax_set_subst_params(old_partition, 0, subst_params);
  corax_set_category_rates(old_partition, rate_cats);

  travbuffer = (corax_unode_t **)malloc(nodes_count * sizeof(corax_unode_t *));
  branch_lengths = (double *)malloc(branch_count * sizeof(double));
  matrix_indices = (unsigned int *)malloc(branch_count * sizeof(int));
  operations = (corax_operation_t *)malloc(inner_nodes_count *
                                                sizeof(corax_operation_t));

  unsigned int traversal_size;

  corax_unode_t * node = tree;

  if (!corax_utree_traverse(node,
                          CORAX_TREE_TRAVERSE_POSTORDER,
                          cb_full_traversal,
                          travbuffer,
                          &traversal_size))
    fatal("Function corax_utree_traverse() requires inner nodes as parameters");

  corax_utree_create_operations(travbuffer,
                              traversal_size,
                              branch_lengths,
                              matrix_indices,
                              operations,
                              &matrix_count,
                              &ops_count);

  printf("\nComputing logL between CLV %d and %d - "
         "(pmatrix %d with branch length %f)\n",
            node->clv_index,
            node->back->clv_index,
            node->pmatrix_index,
            node->length);

  printf ("Traversal size: %d\n", traversal_size);
  printf ("Operations: %d\n", ops_count);
  printf ("Matrices: %d\n", matrix_count);

  corax_update_prob_matrices(old_partition,
                           params_indices,
                           matrix_indices,
                           branch_lengths,
                           matrix_count);

  corax_update_clvs(old_partition, operations, ops_count);

  lk_parent_clv_index = node->clv_index;
  lk_parent_scaler_index = node->scaler_index;
  lk_child_clv_index = node->back->clv_index;
  lk_child_scaler_index = node->back->scaler_index;
  lk_pmatrix_index = node->pmatrix_index;

  logl = corax_compute_edge_loglikelihood(old_partition,
                                        lk_parent_clv_index,
                                        lk_parent_scaler_index,
                                        lk_child_clv_index,
                                        lk_child_scaler_index,
                                        lk_pmatrix_index,
                                        params_indices,
                                        NULL);

  save_logl = logl;
  printf("Log-L: %f\n", logl);

  FILE * bin_file;
  corax_binary_header_t bin_header;
  const char * bin_fname = "test.bin";

  printf("** create binary file\n");
  bin_file = corax_binary_create(bin_fname,
                               &bin_header,
                               CORAX_BIN_ACCESS_RANDOM,
                               2 + old_partition->tips
                               + old_partition->clv_buffers
                               + old_partition->scale_buffers);

  if (!bin_file)
  {
    fatal("Cannot create binary file: %s\n", bin_fname);
  }

  printf("** dump old_partition\n");

  /* We save the structures in an arbitrary order */

  /* IMPORTANT! Attribute CORAX_BIN_ATTRIB_UPDATE_MAP must be set! */
  const unsigned int dump_attr = CORAX_BIN_ATTRIB_UPDATE_MAP | CORAX_BIN_ATTRIB_PARTITION_DUMP_WGT;

  if (!corax_binary_partition_dump(bin_file,
                            BLOCK_ID_PARTITION,
                            old_partition,
                            dump_attr))
  {
    printf("Error dumping old_partition\n");
  }

  /* dump tree */
  if (!corax_binary_utree_dump(bin_file,
                       BLOCK_ID_TREE,
                       tree,
                       tip_nodes_count,
                       CORAX_BIN_ATTRIB_UPDATE_MAP))
  {
    printf("Error dumping tree\n");
  }

  /* dump tipchars */
  int block_id = 0;
  size_t tip_index = 0;
  if (old_partition->attributes & CORAX_ATTRIB_PATTERN_TIP)
  {
    for (tip_index = 0; tip_index < old_partition->tips; tip_index++)
    {
      if(!corax_binary_custom_dump(bin_file,
                                    block_id++,
                                    old_partition->tipchars[tip_index],
                                    old_partition->sites * sizeof(unsigned char),
                                    dump_attr))
      {
        printf("Couldn't dump tipchar number: %lu\n", tip_index);
        exit(1);
      }
    }
  }

  // dump the clvs
  for ( size_t clv_index = tip_index;
        clv_index < old_partition->tips + old_partition->clv_buffers;
        clv_index++)
  {
    if(!corax_binary_clv_dump( bin_file,
                                block_id++,
                                old_partition,
                                clv_index,
                                dump_attr))
    {
      printf("Couldn't dump clv number: %lu\n", clv_index);
      exit(1);
    }
  }

  // dump the scalers
  for (size_t scaler_index = 0; scaler_index < old_partition->scale_buffers; scaler_index++)
  {
    if(!corax_binary_custom_dump(bin_file,
                                  block_id++,
                                  old_partition->scale_buffer[scaler_index],
                                  get_sites_number_scaler(old_partition, scaler_index)  * sizeof(unsigned int), 
                                  dump_attr))
    {
      printf("Couldn't dump scaler number: %lu\n", scaler_index);
      exit(1);
    }

  }

  printf("number of blocks dumped: %d\n", block_id + 2);

  printf("** close binary file\n");

  corax_binary_close(bin_file);

  // corax_utree_show_ascii(tree, (1<<5)-1);

  /* remember important values of the dumped part */
  const unsigned int old_clv_buffers = old_partition->clv_buffers;
  const unsigned int old_tips = old_partition->tips;
  const unsigned int old_scale_buffers = old_partition->scale_buffers;

  /* clean */
  corax_utree_destroy(parsed_tree, NULL);

  printf("\n\n");


  /* reload */
  printf("** reload data\n");
  corax_binary_header_t input_header;
  unsigned int bin_attributes = CORAX_BIN_ATTRIB_PARTITION_LOAD_SKELETON;
  corax_block_map_t * block_map;
  unsigned int n_blocks;

  bin_file = corax_binary_open(bin_fname, &input_header);

  block_map = corax_binary_get_map(bin_file, &n_blocks);

  printf("There are %d blocks in the map\n", n_blocks);
  int partition_offset = get_offset(n_blocks, block_map, BLOCK_ID_PARTITION);

  /* For the offset we can use the actual offset (from the block_map),
     or CORAX_BIN_ACCESS_SEEK */
  partition = corax_binary_partition_load(bin_file,
                                        BLOCK_ID_PARTITION,
                                        NULL, /* in order to create a new partition */
                                        &bin_attributes,
                                        partition_offset);

  if (!partition)
    printf("Error loading partition\n");

  /* check that num clv/scaler/tipchars are correct */
  if (old_clv_buffers != partition->clv_buffers)
  {
    printf("partition->clv_buffers set incorrectly\n");
    exit(1);
  }

  if (old_tips != partition->tips)
  {
    printf("partition->tips set incorrectly\n");
    exit(1);
  }

  if (old_scale_buffers != partition->scale_buffers)
  {
    printf("partition->scale_buffers set incorrectly\n");
    exit(1);
  }

  /* check that pointers are actually null */
  tip_index = 0;
  for (size_t tip_index = 0;
    tip_index < partition->tips; ++tip_index)
  {
    if (partition->clv[tip_index] != NULL)
    {
      printf("skeleton tipchar pointer number %lu was not NULL\n", tip_index);
      exit(1);
    }
  }

  for (size_t i = tip_index; i < partition->clv_buffers + partition->tips; ++i)
  {
    if (partition->clv[i] != NULL)
    {
      printf("skeleton clv pointer number %lu was not NULL\n", i);
      exit(1);
    }
  }

  for (size_t i = 0; i < partition->scale_buffers; ++i)
  {
    if (partition->scale_buffer[i] != NULL)
    {
      printf("skeleton scale_buffer pointer number %lu was not NULL\n", i);
      exit(1);
    }
  }

  /* compare tipchars */
  block_id = 0;
  tip_index = 0;
  if (old_partition->attributes & CORAX_ATTRIB_PATTERN_TIP)
  {
    for (tip_index = 0; tip_index < old_partition->tips; tip_index++)
    {
      size_t size;
      unsigned int type, attribs;
      unsigned char* ptr = (unsigned char*) corax_binary_custom_load(
                                              bin_file,
                                              0,
                                              &size,
                                              &type,
                                              &attribs,
                                              get_offset(n_blocks, block_map, tip_index));
      if (!ptr)
      {
        printf("Error loading tipchar number: %lu\n", tip_index);
        exit(1);
      }

      partition->tipchars[tip_index] = (unsigned char*)ptr;

      if (memcmp( old_partition->tipchars[tip_index],
                  partition->tipchars[tip_index],
                  partition->sites * sizeof(unsigned char)))
      {
        printf("Error! tipchar #%lu does not agree\n", tip_index);
        printf("saved : %.10s\n", old_partition->tipchars[tip_index]);
        printf("loaded: %.10s\n", partition->tipchars[tip_index]);
        exit(1);
      }
    }
  }

  const size_t clvs_and_tips = partition->tips + partition->clv_buffers;
  const double tol = 1e-14;

  /* compare CLVs */
  for (size_t clv_index = tip_index; clv_index < clvs_and_tips; clv_index++)
  {
    unsigned int attribs;
    size_t clv_size = corax_get_clv_size(partition, clv_index);
    if (partition->attributes & CORAX_ATTRIB_SITE_REPEATS)
    {
      partition->repeats->reallocate_repeats(partition, clv_index, CORAX_SCALE_BUFFER_NONE, 
          corax_get_sites_number(partition, clv_index));
    }
    else
    {
      partition->clv[clv_index] = (double*) corax_aligned_alloc(clv_size*sizeof(double), partition->alignment);
    }
    if (!partition->clv[clv_index])
    {
      printf("Error allocating clv number: %lu\n", clv_index);
      exit(1);
    }

    int err = corax_binary_clv_load(bin_file,
                                      0,
                                      partition,
                                      clv_index,
                                      &attribs,
                                      get_offset(n_blocks, block_map, clv_index));

    if (err != CORAX_SUCCESS)
    {
      printf("Error loading clv number: %lu\n", clv_index);
      exit(1);
    }

    for (size_t i = 0; i < clv_size; ++i)
    {
      if (fabs(old_partition->clv[clv_index][i] - partition->clv[clv_index][i]) > tol)
      {
        printf("CLV index %lu does not agree. %.20f vs %.20f\n",
          clv_index, old_partition->clv[clv_index][i], partition->clv[clv_index][i]);
        exit(1);
      }
    }
  }

  /* compare scalers */
  for (size_t scaler_index = 0; scaler_index < partition->scale_buffers; scaler_index++)
  {
    size_t size;
    unsigned int type, attribs;
    unsigned int* ptr = (unsigned int*) corax_binary_custom_load(
                                            bin_file,
                                            0,
                                            &size,
                                            &type,
                                            &attribs,
                                            get_offset(n_blocks, block_map, clvs_and_tips + scaler_index));
    if (!ptr)
    {
      printf("Error loading scaler number: %lu\n", scaler_index);
      exit(1);
    }

    partition->scale_buffer[scaler_index] = ptr;

    if (memcmp( old_partition->scale_buffer[scaler_index],
                partition->scale_buffer[scaler_index], 
                get_sites_number_scaler(partition, scaler_index) * sizeof(unsigned int)))
    {
      printf("Error! scaler #%lu does not agree\n", scaler_index);
      exit(1);
    }
  }

  /* first check if partition was restored correctly */
  logl = corax_compute_edge_loglikelihood(partition,
                                         lk_parent_clv_index,
                                         lk_parent_scaler_index,
                                         lk_child_clv_index,
                                         lk_child_scaler_index,
                                         lk_pmatrix_index,
                                         params_indices,
                                         NULL);

   printf("Restored Log-L: %f\n", logl);
   if (fabs(logl - save_logl) < 1e-7)
     printf("Likelihoods OK!\n");
   else
     fatal("Error: Saved and loaded logL do not agree!!\n");

  /* new we try with ACCESS_SEEK instead of the value taken from the map */
  tree = corax_binary_utree_load(bin_file,
                               BLOCK_ID_TREE,
                               &bin_attributes,
                               CORAX_BIN_ACCESS_SEEK);
  if (!tree)
    fatal("Error loading tree!\n");

  corax_binary_close(bin_file);

  if (!corax_utree_traverse(tree,
                          CORAX_TREE_TRAVERSE_POSTORDER,
                          cb_full_traversal,
                          travbuffer,
                          &traversal_size))
    fatal("Function corax_utree_traverse() requires inner nodes as parameters");

  corax_utree_create_operations(travbuffer,
                              traversal_size,
                              branch_lengths,
                              matrix_indices,
                              operations,
                              &matrix_count,
                              &ops_count);

  corax_update_prob_matrices(partition,
                           params_indices,
                           matrix_indices,
                           branch_lengths,
                           matrix_count);

  corax_update_clvs(partition, operations, ops_count);

  logl = corax_compute_edge_loglikelihood(partition,
                                        tree->clv_index,
                                        tree->scaler_index,
                                        tree->back->clv_index,
                                        tree->back->scaler_index,
                                        tree->pmatrix_index,
                                        params_indices,
                                        NULL);

  printf("Recomputed Log-L: %f\n", logl);

  if (fabs(logl - save_logl) < 1e-7)
    printf("Likelihoods OK!\n");
  else
    fatal("Error: Saved and loaded logL do not agree!!\n");

  //corax_utree_show_ascii(tree, (1<<5)-1);

  /* clean */
  corax_partition_destroy(partition);
  corax_partition_destroy(old_partition);
  corax_utree_graph_destroy(tree, NULL);

  free(block_map);
  free(travbuffer);
  free(branch_lengths);
  free(matrix_indices);
  free(operations);

  return 0;
}
