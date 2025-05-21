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

#define N_STATES_AA 20
#define N_CAT_GAMMA 4
#define FLOAT_PRECISION 4

static double       alpha              = 0.5;
static unsigned int n_cat_gamma        = N_CAT_GAMMA;
unsigned int        params_indices[16] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

int main(int argc, char *argv[])
{
  unsigned int     j;
  double           lk_score;
  unsigned int     n_sites = 12;
  unsigned int     n_tips  = 5;
  double           rate_cats[N_CAT_GAMMA];
  corax_operation_t *operations;
  double *         persite_lnl = (double *)malloc(n_sites * sizeof(double));
  double           checksum;
  int              return_val;

  operations = (corax_operation_t *)malloc(4 * sizeof(corax_operation_t));

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

  operations[3].parent_clv_index    = 8;
  operations[3].child1_clv_index    = 7;
  operations[3].child2_clv_index    = 6;
  operations[3].child1_matrix_index = 2;
  operations[3].child2_matrix_index = 3;
  operations[3].parent_scaler_index = CORAX_SCALE_BUFFER_NONE;
  operations[3].child1_scaler_index = CORAX_SCALE_BUFFER_NONE;
  operations[3].child2_scaler_index = CORAX_SCALE_BUFFER_NONE;

  /* check attributes */
  unsigned int attributes = get_attributes(argc, argv);

  corax_partition_t *partition;
  partition = corax_partition_create(n_tips,      /* numer of tips */
                                   4,           /* clv buffers */
                                   N_STATES_AA, /* number of states */
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

  double       branch_lengths[4] = {0.5, 0.5, 0.3, 0.2};
  unsigned int matrix_indices[4] = {0, 1, 2, 3};

  if (corax_compute_gamma_cats(
          alpha, n_cat_gamma, rate_cats, CORAX_GAMMA_RATES_MEAN)
      == CORAX_FAILURE)
  {
    printf("Error %d: %s\n", corax_errno, corax_errmsg);
    fatal("Fail computing gamma cats");
  }

  corax_set_frequencies(partition, 0, corax_aa_freqs_dayhoff);
  corax_set_subst_params(partition, 0, corax_aa_rates_dayhoff);

  return_val = CORAX_SUCCESS;
  return_val &= corax_set_tip_states(partition, 0, corax_map_aa, "PIGLRVTLRRDRMWI");
  return_val &= corax_set_tip_states(partition, 1, corax_map_aa, "IQGMDITIVT-----");
  return_val &= corax_set_tip_states(partition, 2, corax_map_aa, "--AFALLQKIGMPFE");
  return_val &= corax_set_tip_states(partition, 3, corax_map_aa, "MDISIVT------TA");
  return_val &= corax_set_tip_states(partition, 4, corax_map_aa, "GLSEQTVFHEIDQDK");

  if (!return_val) fatal("Error setting tip states");

  corax_set_category_rates(partition, rate_cats);

  corax_update_prob_matrices(
      partition, params_indices, matrix_indices, branch_lengths, 4);
  corax_update_clvs(partition, operations, 4);

  for (j = 0; j < 4; ++j)
  {
    printf("[%d] P-matrix for branch length %f\n", j + 1, branch_lengths[j]);
    corax_show_pmatrix(partition, j, FLOAT_PRECISION);
    printf("\n");
  }

  /* show CLVs */
  printf("[5] CLV 5: ");
  corax_show_clv(partition, 5, CORAX_SCALE_BUFFER_NONE, FLOAT_PRECISION + 1);
  printf("[6] CLV 6: ");
  corax_show_clv(partition, 6, CORAX_SCALE_BUFFER_NONE, FLOAT_PRECISION + 1);
  printf("[7] CLV 7: ");
  corax_show_clv(partition, 7, CORAX_SCALE_BUFFER_NONE, FLOAT_PRECISION + 1);

  lk_score = corax_compute_root_loglikelihood(
      partition, 8, CORAX_SCALE_BUFFER_NONE, params_indices, persite_lnl);

  printf("\n");
  printf("inner-inner logL: %.6f\n", lk_score);
  printf("persite logL:     ");
  checksum = 0.0;
  for (int i = 0; i < n_sites; i++)
  {
    checksum += persite_lnl[i];
    printf("%.7f  ", persite_lnl[i]);
  }
  printf("\n");
  printf("checksum logL:    %.6f\n", checksum);

  /* move to tip inner */

  operations[0].parent_clv_index    = 7;
  operations[0].child1_clv_index    = 6;
  operations[0].child2_clv_index    = 3;
  operations[0].child1_matrix_index = 0;
  operations[0].child2_matrix_index = 1;
  operations[0].parent_scaler_index = CORAX_SCALE_BUFFER_NONE;
  operations[0].child1_scaler_index = CORAX_SCALE_BUFFER_NONE;
  operations[0].child2_scaler_index = CORAX_SCALE_BUFFER_NONE;

  operations[1].parent_clv_index    = 8;
  operations[1].child1_clv_index    = 7;
  operations[1].child2_clv_index    = 4;
  operations[1].child1_matrix_index = 2;
  operations[1].child2_matrix_index = 3;
  operations[1].parent_scaler_index = CORAX_SCALE_BUFFER_NONE;
  operations[1].child1_scaler_index = CORAX_SCALE_BUFFER_NONE;
  operations[1].child2_scaler_index = CORAX_SCALE_BUFFER_NONE;

  corax_update_clvs(partition, operations, 2);

  lk_score = corax_compute_root_loglikelihood(
      partition, 8, CORAX_SCALE_BUFFER_NONE, params_indices, persite_lnl);

  printf("tip-inner logL:   %.6f\n", lk_score);
  printf("persite logL:     ");
  checksum = 0.0;
  for (int i = 0; i < n_sites; i++)
  {
    checksum += persite_lnl[i];
    printf("%.7f  ", persite_lnl[i]);
  }
  printf("\n");
  printf("checksum logL:    %.6f\n", checksum);

  corax_partition_destroy(partition);
  free(persite_lnl);
  free(operations);

  return (0);
}
