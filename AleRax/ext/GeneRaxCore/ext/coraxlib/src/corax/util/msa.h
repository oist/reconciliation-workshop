/*
 Copyright (C) 2016-2021 Diego Darriba, Alexey Kozlov

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
#ifndef CORAX_UTIL_MSA_H_
#define CORAX_UTIL_MSA_H_

#include "corax/core/partition.h"

#define CORAX_MSA_STATS_NONE (0ul)
#define CORAX_MSA_STATS_DUP_TAXA (1ul << 0)
#define CORAX_MSA_STATS_DUP_SEQS (1ul << 1)
#define CORAX_MSA_STATS_GAP_PROP (1ul << 2)
#define CORAX_MSA_STATS_GAP_SEQS (1ul << 3)
#define CORAX_MSA_STATS_GAP_COLS (1ul << 4)
#define CORAX_MSA_STATS_INV_PROP (1ul << 5)
#define CORAX_MSA_STATS_INV_COLS (1ul << 6)
#define CORAX_MSA_STATS_FREQS (1ul << 7)
#define CORAX_MSA_STATS_SUBST_RATES (1ul << 8)
#define CORAX_MSA_STATS_ENTROPY (1ul << 9)
#define CORAX_MSA_STATS_ALL (~0ul)

#define CORAX_MSA_MAX_ERRORS 100

/* multiple sequence alignment */
typedef struct corax_msa_s
{
  int count;
  int length;

  char **sequence;
  char **label;
} corax_msa_t;

typedef struct msa_stats
{
  unsigned int states;

  unsigned long  dup_taxa_pairs_count;
  unsigned long *dup_taxa_pairs;

  unsigned long  dup_seqs_pairs_count;
  unsigned long *dup_seqs_pairs;

  double         gap_prop;
  unsigned long  gap_seqs_count;
  unsigned long *gap_seqs;
  unsigned long  gap_cols_count;
  unsigned long *gap_cols;

  double         inv_prop;
  unsigned long  inv_cols_count;
  unsigned long *inv_cols;

  double *freqs;
  double *subst_rates;

  double *column_entropies;
  double entropy;
} corax_msa_stats_t;

typedef struct msa_errors
{
  unsigned long  invalid_char_count;
  char *         invalid_chars;
  unsigned long *invalid_char_seq;
  unsigned long *invalid_char_pos;
  int            status;
} corax_msa_errors_t;

#ifdef __cplusplus
extern "C"
{
#endif

  CORAX_EXPORT double *
  corax_msa_empirical_frequencies(const corax_partition_t *partition);
  CORAX_EXPORT double *
  corax_msa_empirical_subst_rates(const corax_partition_t *partition);
  CORAX_EXPORT double
  corax_msa_empirical_invariant_sites(corax_partition_t *partition);

  /**
     * Takes a multiple sequence alignment, the number of states and a tipmap and
     * computes the Shannon entropies per alignment site (in bits).
     *
     * @param msa Multiple Sequence Alignment
     * @param states Number of states (e.g., DNA=4, AA=20 etc.)
     * @param tipmap Mapping from chars to states (e.g., corax_map_nt for DNA)
     * @return Array of per-site Shannon entropies (in bits). Each per-site entropy is >= 0.
   */
  CORAX_EXPORT double *
  corax_msa_column_entropies(const corax_msa_t *   msa,
                             unsigned int          states,
                             const corax_state_t * tipmap
  );

  /**
   * Takes a multiple sequence alignment, the number of states and a tipmap and
   * compute the Shannon entropy of the alignment.
   *
   * @param msa Multiple Sequence Alignment
   * @param states Number of states (e.g., DNA=4, AA=20 etc.)
   * @param tipmap Mapping from chars to states (e.g., corax_map_nt for DNA)
   * @return Shannon entropy of the multiple sequence alignment (in bits). The entropy is >= 0.
   */
  CORAX_EXPORT double
  corax_msa_entropy(const corax_msa_t *         msa,
                    unsigned int          states,
                    const corax_state_t * tipmap
  );

  /**
   * Takes a multiple sequence alignment and a tipmap and computes an entropy like
   * metric based on the number and frequency of patterns in the MSA.
   *
   * @param msa Multiple Sequence Alignment. Note that the MSA object is modified during
   *            the computation of the multinomial statistic.
   * @param tipmap Mapping from chars to states (e.g., corax_map_nt for DNA)
   * @return Pattern entropy. A pattern-based entropy like metric. The pattern entropy is >= 0.
   */
  CORAX_EXPORT double
  corax_msa_pattern_entropy(const corax_msa_t   *msa,
                            const unsigned int  *site_weights,
                            const corax_state_t *tipmap
  );

  /**
   * Takes a multiple sequence alignment and a tipmap and computes the bollback multinomial
   * statistic according to Bollback, JP: Bayesian model adequacy and choice in phylogenetics (2002).
   *
   * @param msa Multiple Sequence Alignment. Note that the MSA object is modified during
   *            the computation of the multinomial statistic.
   * @param tipmap Mapping from chars to states (e.g., corax_map_nt for DNA)
   * @return Bollback multinomial statistic. The bollback multinomial statistic is always <= 0.
   */
  CORAX_EXPORT double
  corax_msa_bollback_multinomial(const corax_msa_t   *msa,
                                 const unsigned int  *site_weights,
                                 const corax_state_t *tipmap
  );

  CORAX_EXPORT corax_msa_errors_t *corax_msa_check(const corax_msa_t *  msa,
                                                   const corax_state_t *tipmap);

  CORAX_EXPORT void corax_msa_destroy_errors(corax_msa_errors_t *errs);

  CORAX_EXPORT corax_msa_stats_t *
               corax_msa_compute_stats(const corax_msa_t *  msa,
                                       unsigned int         states,
                                       const corax_state_t *tipmap,
                                       const unsigned int * weights,
                                       unsigned long        stats_mask);

  CORAX_EXPORT void corax_msa_destroy_stats(corax_msa_stats_t *stats);

  CORAX_EXPORT corax_msa_t *corax_msa_filter(corax_msa_t *  msa,
                                             unsigned long *remove_seqs,
                                             unsigned long  remove_seqs_count,
                                             unsigned long *remove_cols,
                                             unsigned long  remove_cols_count,
                                             unsigned int   inplace);

  CORAX_EXPORT corax_msa_t **corax_msa_split(const corax_msa_t * msa,
                                             const unsigned int *site_part,
                                             unsigned int        part_count);

#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CORAX_UTIL_MSA_H_ */
