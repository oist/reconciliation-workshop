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

#define N_STATES_NT 4
#define N_CAT_GAMMA 4
#define FLOAT_PRECISION 4

static double       titv               = 2.5;
static double       alpha              = 0.5;
static unsigned int n_cat_gamma        = N_CAT_GAMMA;
unsigned int        params_indices[16] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

void test_bad_alpha()
{
  double rate_cats[N_CAT_GAMMA];

  /* test illegal alpha value */
  double invalid_alpha = 0;
  if (corax_compute_gamma_cats(
          invalid_alpha, N_CAT_GAMMA, rate_cats, CORAX_GAMMA_RATES_MEAN)
      == CORAX_FAILURE)
  {
    if (corax_errno != CORAX_ERROR_INVALID_PARAM)
      printf("Error is %d instead of %d\n", corax_errno, CORAX_ERROR_INVALID_PARAM);
  }
  else
  {
    printf("Computing gamma rates for alpha = %f should have failed\n",
           invalid_alpha);
  }
}

double test_lk(corax_partition_t *partition, int gamma_mode)
{
  double          branch_lengths[4] = {0.1, 0.2, 1, 1};
  unsigned int    matrix_indices[4] = {0, 1, 2, 3};
  double          rate_cats[N_CAT_GAMMA];
  unsigned int    j;
  double          lk_score;
  unsigned int    n_sites     = partition->sites;
  double *        persite_lnl = (double *)malloc(n_sites * sizeof(double));
  double          checksum;
  char            prefix[10];
  corax_operation_t operations[4];

  strcpy(prefix, gamma_mode == CORAX_GAMMA_RATES_MEDIAN ? "MEDIAN" : "MEAN");

  if (corax_compute_gamma_cats(alpha, n_cat_gamma, rate_cats, gamma_mode)
      == CORAX_FAILURE)
  {
    printf("Error %d: %s\n", corax_errno, corax_errmsg);
    fatal("Fail computing gamma cats");
  }

  printf("[%s] Discrete GAMMA rates: ", prefix);
  for (j = 0; j < 4; ++j) { printf("%.6lf ", rate_cats[j]); }
  printf("\n");

  corax_set_category_rates(partition, rate_cats);

  operations[0].parent_clv_index    = 5;
  operations[0].child1_clv_index    = 0;
  operations[0].child2_clv_index    = 1;
  operations[0].child1_matrix_index = 1;
  operations[0].child2_matrix_index = 1;
  operations[0].parent_scaler_index = CORAX_SCALE_BUFFER_NONE;
  operations[0].child1_scaler_index = CORAX_SCALE_BUFFER_NONE;
  operations[0].child2_scaler_index = CORAX_SCALE_BUFFER_NONE;

  operations[1].parent_clv_index    = 6;
  operations[1].child1_clv_index    = 5;
  operations[1].child2_clv_index    = 2;
  operations[1].child1_matrix_index = 0;
  operations[1].child2_matrix_index = 1;
  operations[1].parent_scaler_index = CORAX_SCALE_BUFFER_NONE;
  operations[1].child1_scaler_index = CORAX_SCALE_BUFFER_NONE;
  operations[1].child2_scaler_index = CORAX_SCALE_BUFFER_NONE;

  operations[2].parent_clv_index    = 7;
  operations[2].child1_clv_index    = 3;
  operations[2].child2_clv_index    = 4;
  operations[2].child1_matrix_index = 1;
  operations[2].child2_matrix_index = 1;
  operations[2].parent_scaler_index = CORAX_SCALE_BUFFER_NONE;
  operations[2].child1_scaler_index = CORAX_SCALE_BUFFER_NONE;
  operations[2].child2_scaler_index = CORAX_SCALE_BUFFER_NONE;

  corax_update_prob_matrices(
      partition, params_indices, matrix_indices, branch_lengths, 4);
  corax_update_clvs(partition, operations, 3);

  for (j = 0; j < 4; ++j)
  {
    printf("[%s][%d] P-matrix for branch length %f\n",
           prefix,
           j + 1,
           branch_lengths[j]);
    corax_show_pmatrix(partition, j, FLOAT_PRECISION);
    printf("\n");
  }

  /* show CLVs */
  printf("[%s][5] CLV 5: ", prefix);
  corax_show_clv(partition, 5, CORAX_SCALE_BUFFER_NONE, FLOAT_PRECISION + 1);
  printf("\n[%s][6] CLV 6: ", prefix);
  corax_show_clv(partition, 6, CORAX_SCALE_BUFFER_NONE, FLOAT_PRECISION + 1);
  printf("\n[%s][7] CLV 7: ", prefix);
  corax_show_clv(partition, 7, CORAX_SCALE_BUFFER_NONE, FLOAT_PRECISION + 1);

  lk_score = corax_compute_edge_loglikelihood(partition,
                                            6,
                                            CORAX_SCALE_BUFFER_NONE,
                                            7,
                                            CORAX_SCALE_BUFFER_NONE,
                                            0,
                                            params_indices,
                                            persite_lnl);

  printf("\n");
  printf("[%s] inner-inner logL: %.6f\n", prefix, lk_score);
  printf("[%s] persite logL:     ", prefix);
  checksum = 0.0;
  for (int i = 0; i < n_sites; i++)
  {
    checksum += persite_lnl[i];
    printf("%.7f  ", persite_lnl[i]);
  }
  printf("\n");
  printf("[%s] checksum logL:    %.6f\n", prefix, checksum);

  /* move to tip inner */

  operations[0].parent_clv_index    = 7;
  operations[0].child1_clv_index    = 6;
  operations[0].child2_clv_index    = 3;
  operations[0].child1_matrix_index = 0;
  operations[0].child2_matrix_index = 1;
  operations[0].parent_scaler_index = CORAX_SCALE_BUFFER_NONE;
  operations[0].child1_scaler_index = CORAX_SCALE_BUFFER_NONE;
  operations[0].child2_scaler_index = CORAX_SCALE_BUFFER_NONE;

  corax_update_clvs(partition, operations, 1);

  lk_score = corax_compute_edge_loglikelihood(partition,
                                            7,
                                            CORAX_SCALE_BUFFER_NONE,
                                            4,
                                            CORAX_SCALE_BUFFER_NONE,
                                            1,
                                            params_indices,
                                            persite_lnl);

  printf("[%s] tip-inner logL:   %.6f\n", prefix, lk_score);
  printf("[%s] persite logL:     ", prefix);
  checksum = 0.0;
  for (int i = 0; i < n_sites; i++)
  {
    checksum += persite_lnl[i];
    printf("%.7f  ", persite_lnl[i]);
  }
  printf("\n");
  printf("[%s] checksum logL:    %.6f\n", prefix, checksum);
  printf("\n");

  free(persite_lnl);

  return lk_score;
}

int main(int argc, char *argv[])
{
  unsigned int n_sites = 12;
  unsigned int n_tips  = 5;
  int          return_val;
  double       avg_loglh;
  double       med_loglh;

  /* check attributes */
  unsigned int attributes = get_attributes(argc, argv);

  corax_partition_t *partition;
  partition = corax_partition_create(n_tips,      /* numer of tips */
                                   4,           /* clv buffers */
                                   N_STATES_NT, /* number of states */
                                   n_sites,     /* sequence length */
                                   1,           /* different rate parameters */
                                   2 * n_tips - 3, /* probability matrices */
                                   n_cat_gamma,    /* gamma categories */
                                   0,              /* scale buffers */
                                   attributes);    /* attributes */

  if (!partition)
  {
    printf("Error %d: %s\n", corax_errno, corax_errmsg);
    fatal("Fail creating partition");
  }

  double frequencies[4]  = {0.3, 0.4, 0.1, 0.2};
  double subst_params[6] = {1, titv, 1, 1, titv, 1};

  corax_set_frequencies(partition, 0, frequencies);
  corax_set_subst_params(partition, 0, subst_params);

  return_val = CORAX_SUCCESS;
  return_val &= corax_set_tip_states(partition, 0, corax_map_nt, "WAC-CTA-ATCT");
  return_val &= corax_set_tip_states(partition, 1, corax_map_nt, "CCC-TTA-ATGT");
  return_val &= corax_set_tip_states(partition, 2, corax_map_nt, "A-C-TAG-CTCT");
  return_val &= corax_set_tip_states(partition, 3, corax_map_nt, "CTCTTAA-A-CG");
  return_val &= corax_set_tip_states(partition, 4, corax_map_nt, "CAC-TCA-A-TG");

  if (!return_val) fatal("Error setting tip states");

  test_bad_alpha();

  /* test MEAN discrete GAMMA rates */
  avg_loglh = test_lk(partition, CORAX_GAMMA_RATES_MEAN);

  /* test MEDIAN discrete GAMMA rates */
  med_loglh = test_lk(partition, CORAX_GAMMA_RATES_MEDIAN);

  printf("logL MEAN: %.6lf, logL MEDIAN: %.6lf\n", avg_loglh, med_loglh);

  corax_partition_destroy(partition);

  return (0);
}
