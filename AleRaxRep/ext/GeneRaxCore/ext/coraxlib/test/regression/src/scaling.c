/*
    Copyright (C) 2017 Alexey Kozlov, Diego Darriba, Tomas Flouri

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

    Contact: Alexey Kozlov <Alexey.Kozlov@h-its.org>,
    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/
#include "common.h"

#define N_STATES_NT 4
#define N_STATES_AA 20
#define N_STATES_ODD 5
#define N_CAT_GAMMA 4
#define N_SITES 17
#define FLOAT_PRECISION 4

#define TREEFILE "testdata/2000.tree"

#define MIN(a, b) (a < b ? a : b)

#define DATATYPE_NT 0
#define DATATYPE_AA 1
#define DATATYPE_ODD 2

static char nt_alphabet[]  = "ACGT-";
static char aa_alphabet[]  = "GALMFWKQESPVICYHRNDT-";
static char odd_alphabet[] = "ABCDE-";

static double       alphas[]           = {0.05, 0.2, 2.0, 99.0};
static size_t       alpha_count        = sizeof(alphas) / sizeof(double);
static double       pinvars[]          = {0.00, 0.1};
static size_t       pinvar_count       = sizeof(pinvars) / sizeof(double);
static unsigned int n_cat_gamma        = N_CAT_GAMMA;
static unsigned int params_indices[16] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
static double base_freqs_nt[4]   = {0.4, 0.4, 0.1, 0.1};
static double subst_params_nt[6] = {0.1, 10., 10., 0.1, 0.1, 1};

static double base_freqs_odd[5]    = {0.3, 0.25, 0.1, 0.2, 0.15};
static double subst_params_odd[10] = {1.452176,
                                      0.937951,
                                      0.462880,
                                      0.617729,
                                      1.745312,
                                      0.937951,
                                      0.462880,
                                      0.617729,
                                      1.745312,
                                      1.000000};

static corax_utree_t *    tree;
static corax_unode_t *    root;
static corax_partition_t *part_noscale_nt, *part_sitescale_nt, *part_ratescale_nt;
static corax_partition_t *part_noscale_aa, *part_sitescale_aa, *part_ratescale_aa;
static corax_partition_t *part_noscale_odd, *part_sitescale_odd,
    *part_ratescale_odd;
static unsigned int     traversal_size, matrix_count, ops_count;
static corax_unode_t **   travbuffer;
static unsigned int *   matrix_indices;
static double *         branch_lengths;
static corax_operation_t *operations;
static double *         persite_lnl;
static double *         sumtable;

unsigned int scaler_idx(const corax_partition_t *p, unsigned int clv_idx)
{
  return (p->scale_buffers > 0 && clv_idx >= p->tips) ? clv_idx - p->tips
                                                      : CORAX_SCALE_BUFFER_NONE;
}

void show_scaler(const corax_partition_t *p, unsigned int clv_idx)
{
  unsigned int *site_id = 0;
  if (corax_repeats_enabled(p) && p->repeats->pernode_ids[clv_idx])
    site_id = p->repeats->pernode_site_id[clv_idx];
  unsigned int scaler = scaler_idx(p, clv_idx);
  if (scaler != CORAX_SCALE_BUFFER_NONE)
  {
    unsigned int s, i, j;
    printf("scaler %u: [ ", scaler);
    for (s = 0; s < p->sites; ++s)
    {
      i = site_id ? site_id[s] : s;
      if (p->attributes & CORAX_ATTRIB_RATE_SCALERS)
      {
        unsigned int *scalev     = p->scale_buffer[scaler] + i * p->rate_cats;
        unsigned int  min_scaler = 1e6;
        for (j = 0; j < p->rate_cats; ++j)
          if (scalev[j] < min_scaler) min_scaler = scalev[j];

        printf("%d+(", min_scaler);
        for (j = 0; j < p->rate_cats; ++j)
        {
          unsigned int s = scalev[j] - min_scaler;
          printf("%u", MIN(50, s));
          if (j < p->rate_cats - 1) printf(" ");
        }
        printf(")  ");
      }
      else
        printf("%u ", p->scale_buffer[scaler][i]);
    }
    printf("]\n");
  }
}

void show_clv(const corax_partition_t *p, unsigned int clv_idx, unsigned int site)
{
  unsigned int i;
  unsigned int clv_span = p->states * p->rate_cats;
  printf("CLV %u site %u (size=%u): [ ", clv_idx, site, clv_span);
  for (i = (site - 1) * clv_span; i < site * clv_span; ++i)
    printf("%e ", p->clv[clv_idx][i]);
  printf("]\n");
}

corax_partition_t *init_partition(unsigned int attrs, int datatype)
{
  unsigned int i, j;

  unsigned int       states      = 0;
  const corax_state_t *map         = NULL;
  const char *       alphabet    = NULL;
  const double *     base_freqs  = NULL;
  const double *     subst_rates = NULL;

  switch (datatype)
  {
  case DATATYPE_NT:
    states      = N_STATES_NT;
    map         = corax_map_nt;
    alphabet    = nt_alphabet;
    base_freqs  = base_freqs_nt;
    subst_rates = subst_params_nt;
    break;
  case DATATYPE_AA:
    states      = N_STATES_AA;
    map         = corax_map_aa;
    alphabet    = aa_alphabet;
    base_freqs  = corax_aa_freqs_lg;
    subst_rates = corax_aa_rates_lg;
    break;
  case DATATYPE_ODD:
    states      = N_STATES_ODD;
    map         = odd5_map;
    alphabet    = odd_alphabet;
    base_freqs  = base_freqs_odd;
    subst_rates = subst_params_odd;
    break;
  default:
    assert(0);
  }

  corax_partition_t *p = corax_partition_create(tree->tip_count,
                                            tree->inner_count,
                                            states,
                                            N_SITES,
                                            1, /* rate matrices */
                                            2 * tree->tip_count - 3,
                                            N_CAT_GAMMA,
                                            tree->inner_count,
                                            attrs);

  if (!p) fatal("ERROR creating partition: %s\n", corax_errmsg);

  size_t len = strlen(alphabet);
  char * seq = (char *)calloc(N_SITES + 1, sizeof(char));
  for (i = 0; i < tree->tip_count; ++i)
  {
    for (j = 0; j < N_SITES; ++j)
    {
      seq[j] = (i < 1500) ? alphabet[j % len] : alphabet[(i + j) % len];
    }

    corax_set_tip_states(p, tree->nodes[i]->clv_index, map, seq);
  }

  free(seq);

  corax_set_frequencies(p, 0, base_freqs);
  corax_set_subst_params(p, 0, subst_rates);

  return p;
}

void init(unsigned int attrs)
{
  unsigned int i;

  tree = corax_utree_parse_newick(TREEFILE);

  if (!tree) fatal("ERROR reading tree file: %s\n", corax_errmsg);

  part_sitescale_nt = init_partition(attrs, DATATYPE_NT);
  part_ratescale_nt =
      init_partition(attrs | CORAX_ATTRIB_RATE_SCALERS, DATATYPE_NT);

  part_sitescale_aa = init_partition(attrs, DATATYPE_AA);
  part_ratescale_aa =
      init_partition(attrs | CORAX_ATTRIB_RATE_SCALERS, DATATYPE_AA);

  part_sitescale_odd = init_partition(attrs, DATATYPE_ODD);
  part_ratescale_odd =
      init_partition(attrs | CORAX_ATTRIB_RATE_SCALERS, DATATYPE_ODD);

  /* build fixed structures */
  unsigned int nodes_count  = tree->inner_count + tree->tip_count;
  unsigned int branch_count = nodes_count - 1;
  travbuffer     = (corax_unode_t **)malloc(nodes_count * sizeof(corax_unode_t *));
  branch_lengths = (double *)malloc(branch_count * sizeof(double));
  matrix_indices = (unsigned int *)malloc(branch_count * sizeof(unsigned int));
  operations =
      (corax_operation_t *)malloc(tree->inner_count * sizeof(corax_operation_t));
  persite_lnl = (double *)malloc(part_sitescale_aa->sites * sizeof(double));
  sumtable    = (double *)corax_aligned_alloc(
      part_sitescale_aa->sites * part_sitescale_aa->rate_cats
          * part_sitescale_aa->states_padded * sizeof(double),
      part_sitescale_aa->alignment);

  root = tree->nodes[tree->tip_count + tree->inner_count - 1];

  /* get full traversal */
  corax_utree_traverse(root,
                     CORAX_TREE_TRAVERSE_POSTORDER,
                     cb_full_traversal,
                     travbuffer,
                     &traversal_size);

  corax_utree_create_operations(travbuffer,
                              traversal_size,
                              branch_lengths,
                              matrix_indices,
                              operations,
                              &matrix_count,
                              &ops_count);

  for (i = 0; i < matrix_count; ++i)
    branch_lengths[i] = (i % 2 == 0) ? 1.0 : 1e-6;
}

void comp_derivatives(corax_partition_t *partition,
                      corax_unode_t *    r,
                      double           brlen,
                      double *         d1,
                      double *         d2)
{
  if (!corax_update_sumtable(partition,
                           r->clv_index,
                           r->back->clv_index,
                           r->scaler_index,
                           r->back->scaler_index,
                           params_indices,
                           sumtable))
  {
    fatal("ERROR computing sumtable: %s\n", corax_errmsg);
  }

  if (!corax_compute_likelihood_derivatives(partition,
                                          r->scaler_index,
                                          r->back->scaler_index,
                                          brlen,
                                          params_indices,
                                          sumtable,
                                          d1,
                                          d2))
  {
    fatal("ERROR computing derivatives: %s\n", corax_errmsg);
  }
}

int eval(corax_partition_t *partition, double alpha, double pinv)
{
  double       rate_cats[N_CAT_GAMMA];
  unsigned int i;
  double       d_f, dd_f;

  corax_compute_gamma_cats(alpha, n_cat_gamma, rate_cats, CORAX_GAMMA_RATES_MEAN);
  corax_set_category_rates(partition, rate_cats);

  for (i = 0; i < partition->rate_matrices; ++i)
    corax_update_invariant_sites_proportion(partition, i, pinv);

  printf("datatype = ");
  if (partition->states == 4)
    printf("DNA");
  else if (partition->states == 20)
    printf("PROT");
  else
    printf("ODD");

  printf(", scaling = ");
  if (partition->attributes & CORAX_ATTRIB_RATE_SCALERS)
    printf("per-rate");
  else if (partition->scale_buffers > 0)
    printf("per-site");
  else
    printf("OFF");

  printf(", alpha = %lf, pinv = %lf, rates = [ ", alpha, pinv);
  for (i = 0; i < n_cat_gamma; ++i) printf("%lf ", rate_cats[i]);
  printf("]\n");

  printf("recompute P-matrices: %d\n", matrix_count);
  corax_update_prob_matrices(
      partition, params_indices, matrix_indices, branch_lengths, matrix_count);

  printf("recompute CLVs: %d\n", ops_count);
  corax_update_clvs(partition, operations, ops_count);

  show_scaler(partition, root->back->clv_index);
  //  show_clv(partition, root->back->clv_index, 53);
  //  corax_show_pmatrix(partition, root->pmatrix_index, 9);
  //  corax_show_pmatrix(partition, root->pmatrix_index-1, 9);

  // test derivatives
  comp_derivatives(partition, root, 1.0, &d_f, &dd_f);

  double ii_loglh = corax_compute_edge_loglikelihood(partition,
                                                   root->clv_index,
                                                   root->scaler_index,
                                                   root->back->clv_index,
                                                   root->back->scaler_index,
                                                   root->pmatrix_index,
                                                   params_indices,
                                                   persite_lnl);

  printf("per-site logLH: [ ");
  for (i = 0; i < partition->sites; ++i) printf("%.4lf ", persite_lnl[i]);
  printf("]\n");
  printf("logLH INNER-INNER at edge (%u-%u): %.4lf, LH derivatives: %.4lf, "
         "%.4lf\n",
         root->clv_index,
         root->back->clv_index,
         ii_loglh,
         d_f,
         dd_f);

  corax_unode_t *new_root = root->next->next;
  corax_unode_t *tip      = new_root->back;
  assert(tip->clv_index < tree->tip_count);

  // change root CLV orientation
  corax_operation_t op;
  op.parent_clv_index    = new_root->clv_index;
  op.child1_clv_index    = new_root->next->back->clv_index;
  op.child2_clv_index    = new_root->next->next->back->clv_index;
  op.child1_matrix_index = new_root->next->pmatrix_index;
  op.child2_matrix_index = new_root->next->next->pmatrix_index;
  op.parent_scaler_index = scaler_idx(partition, op.parent_clv_index);
  op.child1_scaler_index = scaler_idx(partition, op.child1_clv_index);
  op.child2_scaler_index = scaler_idx(partition, op.child2_clv_index);

  corax_update_clvs(partition, &op, 1);

  // test derivatives
  comp_derivatives(partition, new_root, 1.0, &d_f, &dd_f);

  double ti_loglh = corax_compute_edge_loglikelihood(partition,
                                                   new_root->clv_index,
                                                   new_root->scaler_index,
                                                   tip->clv_index,
                                                   tip->scaler_index,
                                                   tip->pmatrix_index,
                                                   params_indices,
                                                   persite_lnl);

  printf("logLH TIP-INNER at edge (%u-%u): %.4lf, LH derivatives: %.4lf, "
         "%.4lf\n\n",
         new_root->clv_index,
         new_root->back->clv_index,
         ti_loglh,
         d_f,
         dd_f);

  if (fabs(ii_loglh - ti_loglh) > 1e-4)
  {
    printf("ERROR: Likelihood mismatch INNER-INNER vs TIP-INNER (see above)\n");
  }

  return 1;
}

void cleanup()
{
  free(travbuffer);
  free(branch_lengths);
  free(operations);
  free(matrix_indices);
  free(persite_lnl);
  free(sumtable);
  corax_partition_destroy(part_noscale_nt);
  corax_partition_destroy(part_sitescale_nt);
  corax_partition_destroy(part_ratescale_nt);
  corax_partition_destroy(part_noscale_aa);
  corax_partition_destroy(part_sitescale_aa);
  corax_partition_destroy(part_ratescale_aa);
  corax_partition_destroy(part_noscale_odd);
  corax_partition_destroy(part_sitescale_odd);
  corax_partition_destroy(part_ratescale_odd);
  corax_utree_destroy(tree, NULL);
}

int main(int argc, char *argv[])
{
  /* check attributes */
  unsigned int attributes = get_attributes(argc, argv);

  init(attributes);

  unsigned int i, j;
  for (i = 0; i < alpha_count; ++i)
  {
    for (j = 0; j < pinvar_count; ++j)
    {
      double alpha = alphas[i];
      double pinv  = pinvars[j];
      eval(part_sitescale_nt, alpha, pinv);
      eval(part_ratescale_nt, alpha, pinv);
      eval(part_sitescale_aa, alpha, pinv);
      eval(part_ratescale_aa, alpha, pinv);
      eval(part_sitescale_odd, alpha, pinv);
      eval(part_ratescale_odd, alpha, pinv);
    }
  }

  cleanup();

  return (0);
}
