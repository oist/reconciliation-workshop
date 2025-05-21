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

#define N_RATE_CATS 4
#define ALPHA 1

#define N_PROT_MODELS 20
#define N_STATES 20

#define FLOAT_PRECISION 5

static const double *prot_matrices[N_PROT_MODELS] = {
    corax_aa_rates_dayhoff,  corax_aa_rates_lg,       corax_aa_rates_dcmut,
    corax_aa_rates_jtt,      corax_aa_rates_mtrev,    corax_aa_rates_wag,
    corax_aa_rates_rtrev,    corax_aa_rates_cprev,    corax_aa_rates_vt,
    corax_aa_rates_blosum62, corax_aa_rates_mtmam,    corax_aa_rates_mtart,
    corax_aa_rates_mtzoa,    corax_aa_rates_pmb,      corax_aa_rates_hivb,
    corax_aa_rates_hivw,     corax_aa_rates_jttdcmut, corax_aa_rates_flu,
    corax_aa_rates_stmtrev,  corax_aa_rates_den};

static const double *prot_freqs[N_PROT_MODELS] = {
    corax_aa_freqs_dayhoff,  corax_aa_freqs_lg,       corax_aa_freqs_dcmut,
    corax_aa_freqs_jtt,      corax_aa_freqs_mtrev,    corax_aa_freqs_wag,
    corax_aa_freqs_rtrev,    corax_aa_freqs_cprev,    corax_aa_freqs_vt,
    corax_aa_freqs_blosum62, corax_aa_freqs_mtmam,    corax_aa_freqs_mtart,
    corax_aa_freqs_mtzoa,    corax_aa_freqs_pmb,      corax_aa_freqs_hivb,
    corax_aa_freqs_hivw,     corax_aa_freqs_jttdcmut, corax_aa_freqs_flu,
    corax_aa_freqs_stmtrev,  corax_aa_freqs_den};

static const char *prot_model_names[N_PROT_MODELS] = {
    "Dayhoff", "LG",   "DCMut",     "JTT",   "MtREV",   "WAG",   "RtREV",
    "CpREV",   "VT",   "Blosum62",  "MtMam", "MtArt",   "MtZoa", "PMB",
    "HIVb",    "HIVw", "JTT-DCMut", "FLU",   "StmtREV", "DEN"};

int main(int argc, char *argv[])
{
  unsigned int i, cur_model;

  corax_partition_t *partition;
  corax_operation_t *operations;
  unsigned int     params_indices[N_RATE_CATS] = {0, 0, 0, 0};

  unsigned int attributes = get_attributes(argc, argv);

  printf("Creating PLL partition\n");

  partition = corax_partition_create(5,           /* tips */
                                   4,           /* clv buffers */
                                   N_STATES,    /* states */
                                   113,         /* sites */
                                   1,           /* different rate parameters */
                                   8,           /* probability matrices */
                                   N_RATE_CATS, /* rate categories */
                                   1,
                                   attributes);

  double       branch_lengths[4] = {0.1, 0.2, 1, 1};
  unsigned int matrix_indices[4] = {0, 1, 2, 3};

  double *rate_cats = (double *)malloc(N_RATE_CATS * sizeof(double));

  if (corax_compute_gamma_cats(
          ALPHA, N_RATE_CATS, rate_cats, CORAX_GAMMA_RATES_MEAN)
      == CORAX_FAILURE)
  {
    printf("Fail computing the gamma rates\n");
    exit(1);
  }
  corax_set_category_rates(partition, rate_cats);
  free(rate_cats);

  corax_set_tip_states(partition,
                     0,
                     corax_map_aa,
                     "PIGLRVTLRRDRMWIFLEKLLNVALPRIRDFRGLN--"
                     "PNSFDGRGNYNLGLREQLIFPEITYDMVDALRGMDIAVVT------TAETDEE----"
                     "------ARALLELLGFPFR");
  corax_set_tip_states(partition,
                     1,
                     corax_map_aa,
                     "PIGLKVTLRGARMYNFLYKLINIVLPKVRDFRGLD--"
                     "PNSFDGRGNYSFGLSEQLVFPELNPDEVRRIQGMDITIVT------TAKTDQE----"
                     "------ARRLLELFGMPFK");
  corax_set_tip_states(partition,
                     2,
                     corax_map_aa,
                     "AIGAKVTLRGKKMYDFLDKLINVALPRVRDFRGVS--"
                     "KTSFDGFGNFYTGIKEQIIFPEVDHDKVIRLRGMDITIVT------SAKTNKE----"
                     "------AFALLQKIGMPFE");
  corax_set_tip_states(partition,
                     3,
                     corax_map_aa,
                     "PIGVMVTLRGDYMYAFLDRLINLSLPRIRDFRGIT--"
                     "AKSFDGRGNYNLGLKEQLIFPEVDYDGIEQIRGMDISIVT------TAKTDQE----"
                     "------GLALLKSLGMPFA");
  corax_set_tip_states(partition,
                     4,
                     corax_map_aa,
                     "PIGTHATLRGDRMWEFLDRLVTLPLPRIRDFRGLS--"
                     "DRQFDGNGNYTFGLSEQTVFHEIDQDKIDRVRGMDITVVT------TAKNDDE----"
                     "------GRALLKALGFPFK");

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

  for (cur_model = 0; cur_model < N_PROT_MODELS; cur_model++)
  {

    printf("\nSetting model %s...\n", prot_model_names[cur_model]);

    corax_set_subst_params(partition, 0, prot_matrices[cur_model]);
    corax_set_frequencies(partition, 0, prot_freqs[cur_model]);

    double sum_freqs = 0.0;
    for (i = 0; i < N_STATES; ++i)
    {
      sum_freqs += partition->frequencies[0][i];
    }
    if (fabs(sum_freqs - 1.0) > 1e-8)
    {
      printf(" WARNING: Freq sum diff: %e\n", sum_freqs - 1.0);
      for (i = 0; i < N_STATES; ++i)
      {
        printf("%f ", partition->frequencies[0][i]);
      }
      printf("\n");
    }

    printf("Updating prob matrices...\n");

    // corax_update_invariant_sites_proportion(partition, 0.17);
    corax_update_prob_matrices(
        partition, params_indices, matrix_indices, branch_lengths, 4);
    for (i = 0; i < 4; ++i)
    {
      printf("P-matrix for branch length %f\n", branch_lengths[i]);
      corax_show_pmatrix(partition, i, FLOAT_PRECISION);
      printf("\n");
    }

    corax_update_clvs(partition, operations, 3);

    printf("CLV 5: ");
    corax_show_clv(partition, 5, CORAX_SCALE_BUFFER_NONE, FLOAT_PRECISION + 1);
    printf("CLV 6: ");
    corax_show_clv(partition, 6, CORAX_SCALE_BUFFER_NONE, FLOAT_PRECISION + 1);
    printf("CLV 7: ");
    corax_show_clv(partition, 7, CORAX_SCALE_BUFFER_NONE, FLOAT_PRECISION + 1);

    double logl = corax_compute_edge_loglikelihood(partition,
                                                 6,
                                                 CORAX_SCALE_BUFFER_NONE,
                                                 7,
                                                 CORAX_SCALE_BUFFER_NONE,
                                                 0,
                                                 params_indices,
                                                 NULL);

    printf("Log-L (%s): %.6f\n", prot_model_names[cur_model], logl);
  }

  free(operations);
  corax_partition_destroy(partition);

  return (0);
}
