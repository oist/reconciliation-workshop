/*
 Copyright (C) 2016 Diego Darriba, Alexey Kozlov

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

/**
 * @file corax_msa.c
 *
 * @brief Operations on multiple sequence alignments
 *
 * @author Diego Darriba
 * @author Alexey Kozlov
 */

#ifdef USE_POSIX_SEARCH
#include <search.h>
#endif

#include "msa.h"

#include "corax/corax_core.h"
#include "corax/corax_model.h"

CORAX_EXPORT double *
corax_msa_empirical_frequencies(const corax_partition_t *partition)
{
  unsigned int         i, j, k, n;
  unsigned int         states        = partition->states;
  unsigned int         states_padded = partition->states_padded;
  unsigned int         sites         = partition->sites;
  unsigned int         rate_cats     = partition->rate_cats;
  unsigned int         tips          = partition->tips;
  const corax_state_t *tipmap        = partition->tipmap;
  const unsigned int  *w             = partition->pattern_weights;
  double              *frequencies;

  if ((frequencies = (double *)calloc((size_t)states, sizeof(double))) == NULL)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Cannot allocate memory for empirical frequencies");
    return NULL;
  }

  if (partition->attributes & CORAX_ATTRIB_PATTERN_TIP)
  {
    if (states == 4)
    {
      for (i = 0; i < tips; ++i)
      {
        const unsigned char *tipchars = partition->tipchars[i];
        for (n = 0; n < sites; ++n)
        {
          unsigned int state    = (unsigned int)tipchars[n];
          double       sum_site = 1.0 * CORAX_POPCNT32(state);
          for (k = 0; k < states; ++k)
          {
            if (state & 1) frequencies[k] += w[n] / sum_site;
            state >>= 1;
          }
        }
      }
    }
    else
    {
      for (i = 0; i < tips; ++i)
      {
        const unsigned char *tipchars = partition->tipchars[i];
        for (n = 0; n < sites; ++n)
        {
          corax_state_t state    = tipmap[(int)tipchars[n]];
          double        sum_site = 1.0 * CORAX_STATE_POPCNT(state);
          for (k = 0; k < states; ++k)
          {
            if (state & 1) frequencies[k] += w[n] / sum_site;
            state >>= 1;
          }
        }
      }
    }
  }
  else
  {
    for (i = 0; i < tips; ++i)
    {
      unsigned int *site_to_id = 0;
      if ((partition->attributes & CORAX_ATTRIB_SITE_REPEATS)
          && partition->repeats->pernode_ids[i])
      {
        site_to_id = partition->repeats->pernode_site_id[i];
      }

      for (n = 0; n < sites; ++n)
      {
        j = site_to_id ? site_to_id[n] : n;
        j *= (states_padded * rate_cats);
        double sum_site = 0.0;
        for (k = 0; k < states; ++k) sum_site += partition->clv[i][j + k];

        for (k = 0; k < states; ++k)
          frequencies[k] += w[n] * partition->clv[i][j + k] / sum_site;
      }
    }
  }

  /* IMPORTANT: we must use the original number of sites in alignment (before
   * pattern compression), since we previously multiplied our base counts with
   * the respective column weights! */
  unsigned int uncomp_sites = 0;
  for (i = 0; i < sites; ++i) uncomp_sites += w[i];

  for (k = 0; k < states; ++k) frequencies[k] /= uncomp_sites * tips;

#ifndef NDEBUG
  double sum_test = 0.0;
  for (k = 0; k < states; ++k) { sum_test += frequencies[k]; }
  assert(fabs(sum_test - 1) < 1e-6);
#endif

  return frequencies;
}

void compute_pair_rates(unsigned int         states,
                        unsigned int         tips,
                        unsigned long        sites,
                        unsigned char      **tipchars,
                        const unsigned int  *w,
                        const corax_state_t *tipmap,
                        size_t              *state_freq,
                        size_t              *pair_rates)
{
  unsigned int  i, j, k;
  unsigned long n;
  corax_state_t undef_state = (corax_state_t)(pow(2, states)) - 1;

  for (n = 0; n < sites; ++n)
  {
    memset(state_freq, 0, sizeof(unsigned) * (states));
    for (i = 0; i < tips; ++i)
    {
      const unsigned int c     = (unsigned int)tipchars[i][n];
      corax_state_t      state = tipmap ? tipmap[c] : c;
      if (state == undef_state) continue;
      for (k = 0; k < states; ++k)
      {
        if (state & 1) state_freq[k]++;
        state >>= 1;
      }
    }

    for (i = 0; i < states; i++)
    {
      if (state_freq[i] == 0) continue;
      for (j = i + 1; j < states; j++)
      {
        pair_rates[i * states + j] += state_freq[i] * state_freq[j] * w[n];
      }
    }
  }
}

CORAX_EXPORT double *
corax_msa_empirical_subst_rates(const corax_partition_t *partition)
{
  unsigned int         i, j, k, n;
  unsigned int         states        = partition->states;
  unsigned int         states_padded = partition->states_padded;
  unsigned int         sites         = partition->sites;
  unsigned int         tips          = partition->tips;
  unsigned int         rate_cats     = partition->rate_cats;
  const corax_state_t *tipmap        = partition->tipmap;
  const unsigned int  *w             = partition->pattern_weights;
  unsigned char      **tipchars      = partition->tipchars;

  unsigned int n_subst_rates = (states * (states - 1) / 2);
  double *subst_rates = (double *)calloc((size_t)n_subst_rates, sizeof(double));

  size_t *pair_rates = (size_t *)calloc(states * states, sizeof(size_t));
  size_t *state_freq = (size_t *)malloc(states * sizeof(size_t));

  if (!(subst_rates && pair_rates && state_freq))
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Cannot allocate memory for empirical subst rates");
    if (subst_rates) free(subst_rates);
    if (pair_rates) free(pair_rates);
    if (state_freq) free(state_freq);
    return NULL;
  }

  if (partition->attributes & CORAX_ATTRIB_PATTERN_TIP)
  {
    compute_pair_rates(
        states, tips, sites, tipchars, w, tipmap, state_freq, pair_rates);
  }
  else
  {
    unsigned int cur_site = 0;
    for (n = 0; n < sites * states_padded * rate_cats;
         n += (states_padded * rate_cats))
    {
      memset(state_freq, 0, sizeof(unsigned) * (states));
      for (i = 0; i < tips; ++i)
      {
        int unstate = 1;
        for (k = 0; k < states; ++k)
          if (partition->clv[i][n + k] < 1e-7)
          {
            unstate = 0;
            break;
          }
        if (unstate) continue;
        for (k = 0; k < states; ++k)
        {
          if (partition->clv[i][n + k] > 0) { state_freq[k]++; }
        }
      }

      for (i = 0; i < states; i++)
      {
        if (state_freq[i] == 0) continue;
        for (j = i + 1; j < states; j++)
        {
          pair_rates[i * states + j] +=
              state_freq[i] * state_freq[j] * w[cur_site];
        }
      }
      cur_site++;
    }
  }

  k                = 0;
  double last_rate = pair_rates[(states - 2) * states + states - 1];
  if (last_rate < 1e-7) last_rate = 1;
  for (i = 0; i < states - 1; i++)
  {
    for (j = i + 1; j < states; j++)
    {
      subst_rates[k++] = pair_rates[i * states + j] / last_rate;
      if (subst_rates[k - 1] < 0.01) subst_rates[k - 1] = 0.01;
      if (subst_rates[k - 1] > 50.0) subst_rates[k - 1] = 50.0;
    }
  }
  subst_rates[k - 1] = 1.0;

  free(state_freq);
  free(pair_rates);

  return subst_rates;
}

CORAX_EXPORT double
corax_msa_empirical_invariant_sites(corax_partition_t *partition)
{
  unsigned int        n;
  unsigned int        n_inv = 0;
  unsigned int        sites = partition->sites;
  const unsigned int *w     = partition->pattern_weights;

  /* reset errno */
  corax_errno = 0;

  if (!partition->invariant)
    if (!corax_update_invariant_sites(partition)) return (double)-INFINITY;

  const int *invariant = partition->invariant;

  unsigned int uncomp_sites = 0;
  for (n = 0; n < sites; ++n)
  {
    if (invariant[n] > -1) n_inv += w[n];
    uncomp_sites += w[n];
  }

  double empirical_pinv = (double)1.0 * n_inv / uncomp_sites;
  return empirical_pinv;
}

CORAX_EXPORT double *corax_msa_column_entropies(const corax_msa_t   *msa,
                                                unsigned int         states,
                                                const corax_state_t *tipmap)
{
  if (!msa)
  {
    corax_set_error(CORAX_ERROR_INVALID_PARAM, "MSA structure is NULL");
    return CORAX_FAILURE;
  }

  if (!tipmap)
  {
    corax_set_error(CORAX_ERROR_INVALID_PARAM,
                    "Character-to-state mapping (charmap) is NULL");
    return CORAX_FAILURE;
  }

  const unsigned long msa_count  = (unsigned long)msa->count;
  const unsigned long msa_length = (unsigned long)msa->length;

  unsigned long i, j, k;

  corax_state_t gap_state = 0;

  double *const column_entropies = (double *)calloc(msa_length, sizeof(double));

  /* gap state has always all bits set to one */
  gap_state = (1ul << (states)) - 1;

  for (j = 0; j < msa_length; ++j) // site loop msa_length
  {
    // compute the column wise entropy
    // for this first count the proportion of each state in this column
    int state_counts[states];
    memset(state_counts, 0, sizeof(unsigned) * (states));

    int count_non_gaps = 0;

    for (i = 0; i < msa_count; ++i) // taxa loop msa_count
    {
      const corax_state_t state  = tipmap[(int)msa->sequence[i][j]];
      const int           is_gap = state == gap_state ? 1 : 0;

      if (is_gap) { continue; }

      for (k = 0; k < states; ++k)
      {
        if (state & (1ll << k)) { state_counts[k]++; }
      }
      count_non_gaps++;
    }

    double column_entropy = 0.0;
    for (k = 0; k < states; ++k)
    {
      if (state_counts[k] == 0) { continue; }

      double state_probability = state_counts[k] / (double)count_non_gaps;
      column_entropy += -state_probability * log2(state_probability);
    }

    column_entropies[j] = column_entropy;
  }

  return column_entropies;
}

CORAX_EXPORT double corax_msa_entropy(const corax_msa_t   *msa,
                                      unsigned int         states,
                                      const corax_state_t *tipmap)
{
  if (!msa)
  {
    corax_set_error(CORAX_ERROR_INVALID_PARAM, "MSA structure is NULL");
    return CORAX_FAILURE;
  }

  if (!tipmap)
  {
    corax_set_error(CORAX_ERROR_INVALID_PARAM,
                    "Character-to-state mapping (charmap) is NULL");
    return CORAX_FAILURE;
  }

  const unsigned long msa_length = (unsigned long)msa->length;

  unsigned long i;

  double *column_entropies = corax_msa_column_entropies(msa, states, tipmap);

  double sum = 0.0;

  for (i = 0; i < msa_length; ++i) { sum += column_entropies[i]; }

  free(column_entropies);

  return sum / msa_length;
}

CORAX_EXPORT double corax_msa_pattern_entropy(const corax_msa_t   *msa,
                                              const unsigned int  *site_weights,
                                              const corax_state_t *tipmap)
{
  if (!msa)
  {
    corax_set_error(CORAX_ERROR_INVALID_PARAM, "MSA structure is NULL");
    return CORAX_FAILURE;
  }

  if (!site_weights)
  {
    corax_set_error(CORAX_ERROR_INVALID_PARAM,
                    "Site weights structure is NULL");
    return CORAX_FAILURE;
  }

  if (!tipmap)
  {
    corax_set_error(CORAX_ERROR_INVALID_PARAM,
                    "Character-to-state mapping (charmap) is NULL");
    return CORAX_FAILURE;
  }

  int i;

  double sum = 0.0;

  for (i = 0; i < msa->length; ++i)
  {
    unsigned int site_weight = site_weights[i];

    if (site_weight == 0) { continue; }

    double factor = site_weight * log(site_weight);
    sum += factor;
  }

  return sum;
}

CORAX_EXPORT double corax_msa_bollback_multinomial(const corax_msa_t   *msa,
                                                   const unsigned int  *site_weights,
                                                   const corax_state_t *tipmap)
{
  int number_of_sites = 0;
  int i;

  for (i = 0; i < msa->length; ++i)
  {
    unsigned int site_weight = site_weights[i];

    if (site_weight == 0) { continue; }

    number_of_sites += site_weight;
  }
  double pattern_entropy = corax_msa_pattern_entropy(msa, site_weights, tipmap);

  return pattern_entropy - number_of_sites * log(number_of_sites);
}

#ifdef USE_POSIX_SEARCH
/* Find duplicates using hash table from search.h. This works best for short
 * strings, so we use this method for checking taxa names */
static int find_duplicate_strings_htable(char **const    strings,
                                         unsigned long   string_count,
                                         unsigned long **duplicates,
                                         unsigned long  *duplicate_count)
{
  if (!strings) return CORAX_FAILURE;

  unsigned long *data =
      (unsigned long *)malloc(string_count * sizeof(unsigned long));

  unsigned long *tmpdup =
      (unsigned long *)malloc(string_count * 2 * sizeof(unsigned long));

  if (!data || !tmpdup)
  {
    if (data) free(data);
    if (tmpdup) free(tmpdup);
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Cannot allocate memory for duplicates array");
    return CORAX_FAILURE;
  }

  hcreate(string_count);

  unsigned long i;
  for (i = 0; i < string_count; ++i)
  {
    data[i] = i;
    ENTRY entry;

#ifdef __APPLE__
    entry.key = strdup(strings[i]);
#else
    entry.key = (char *)strings[i];
#endif

    entry.data   = (void *)(data + i);
    ENTRY *found = hsearch(entry, ENTER);
    if (found->data != entry.data)
    {
      const unsigned long idx1 = *((unsigned long *)found->data);
      const unsigned long idx2 = *((unsigned long *)entry.data);

      /* duplicate found, save it */
      tmpdup[(*duplicate_count) * 2]     = idx1;
      tmpdup[(*duplicate_count) * 2 + 1] = idx2;
      (*duplicate_count)++;
    }
  }

  hdestroy();
  free(data);

  if (*duplicate_count > 0)
  {
    *duplicates = (unsigned long *)realloc(
        tmpdup, (*duplicate_count) * 2 * sizeof(unsigned long));
  }
  else
  {
    free(tmpdup);
    *duplicates = NULL;
  }

  return CORAX_SUCCESS;
}
#endif

/* Find duplicates using custom hash function optimized for long low-variance
 * strings - this method is used for detecting identical sequences */
static int find_duplicate_strings(char **const    strings,
                                  unsigned long   string_count,
                                  unsigned long   string_len,
                                  unsigned long **duplicates,
                                  unsigned long  *duplicate_count)
{
  if (!strings) return CORAX_FAILURE;

  unsigned long *tmpdup =
      (unsigned long *)malloc(string_count * 2 * sizeof(unsigned long));

  unsigned long *hash =
      (unsigned long *)calloc(string_count, sizeof(unsigned long));

  unsigned char *dupflag =
      (unsigned char *)calloc(string_count, sizeof(unsigned char));

  if (!hash || !tmpdup || !dupflag)
  {
    if (hash) free(hash);
    if (tmpdup) free(tmpdup);
    if (dupflag) free(dupflag);
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Cannot allocate memory for duplicates array");
    return CORAX_FAILURE;
  }

  *duplicate_count = 0;

  unsigned long i;
  unsigned long j;

  for (i = 0; i < string_count; ++i)
  {
    for (j = 0; j < (string_len ? string_len : strlen(strings[i])); ++j)
      hash[i] = hash[i] + (j + 1) * (unsigned long)strings[i][j];
  }

  int coll = 0;

  for (i = 0; i < string_count; ++i)
  {
    if (dupflag[i]) continue;
    for (j = i + 1; j < string_count; ++j)
    {
      if (hash[i] == hash[j])
      {
        coll++;
        if (strcmp(strings[i], strings[j]) == 0)
        {
          /* duplicate found, save it */
          tmpdup[(*duplicate_count) * 2]     = i;
          tmpdup[(*duplicate_count) * 2 + 1] = j;
          (*duplicate_count)++;
          dupflag[j] = 1;
        }
      }
    }
  }

  //  printf("Collisions: %d, duplicates: %lu\n", coll, *duplicate_count);

  free(hash);
  free(dupflag);

  if (*duplicate_count > 0)
  {
    *duplicates = (unsigned long *)realloc(
        tmpdup, (*duplicate_count) * 2 * sizeof(unsigned long));

    if (!(*duplicates))
    {
      corax_set_error(CORAX_ERROR_MEM_ALLOC,
                      "Cannot allocate memory for duplicates array");
      return CORAX_FAILURE;
    }
  }
  else
  {
    free(tmpdup);
    *duplicates = NULL;
  }

  return CORAX_SUCCESS;
}

CORAX_EXPORT corax_msa_errors_t *corax_msa_check(const corax_msa_t   *msa,
                                                 const corax_state_t *tipmap)
{
  unsigned long i, j;

  if (!msa)
  {
    corax_set_error(CORAX_ERROR_INVALID_PARAM, "MSA structure is NULL");
    return CORAX_FAILURE;
  }

  if (!tipmap)
  {
    corax_set_error(CORAX_ERROR_INVALID_PARAM,
                    "Character-to-state mapping (charmap) is NULL");
    return CORAX_FAILURE;
  }

  corax_msa_errors_t *errs =
      (corax_msa_errors_t *)calloc(1, sizeof(corax_msa_errors_t));

  if (!errs)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Cannot allocate memory for MSA error structure");
    return NULL;
  }

  errs->status = CORAX_SUCCESS;

  const unsigned long msa_count  = (unsigned long)msa->count;
  const unsigned long msa_length = (unsigned long)msa->length;
  for (i = 0; i < msa_count; ++i)
  {
    const unsigned char *seqchars = (unsigned char *)msa->sequence[i];
    for (j = 0; j < msa_length; ++j)
    {
      const int           c     = (int)seqchars[j];
      const corax_state_t state = tipmap[c];

      if (!state)
      {
        if (!errs->invalid_chars)
        {
          errs->status = CORAX_FAILURE;
          errs->invalid_chars =
              (char *)calloc(CORAX_MSA_MAX_ERRORS, sizeof(char));
          errs->invalid_char_seq = (unsigned long *)calloc(
              CORAX_MSA_MAX_ERRORS, sizeof(unsigned long));
          errs->invalid_char_pos = (unsigned long *)calloc(
              CORAX_MSA_MAX_ERRORS, sizeof(unsigned long));
        }
        errs->invalid_chars[errs->invalid_char_count]    = (char)c;
        errs->invalid_char_seq[errs->invalid_char_count] = i;
        errs->invalid_char_pos[errs->invalid_char_count] = j;
        errs->invalid_char_count++;
        if (errs->invalid_char_count >= CORAX_MSA_MAX_ERRORS) return errs;
      }
    }
  }

  return errs;
}

CORAX_EXPORT void corax_msa_destroy_errors(corax_msa_errors_t *errs)
{
  if (!errs) return;

  if (errs->invalid_chars) free(errs->invalid_chars);

  if (errs->invalid_char_seq) free(errs->invalid_char_seq);

  if (errs->invalid_char_pos) free(errs->invalid_char_pos);

  free(errs);
}

/**
 *  Compute diverse alignment statistics (see @param stats_mask for details)
 *
 *  @param msa Multiple Sequence Alignment
 *  @param states Number of states (e.g., DNA=4, AA=20 etc.)
 *  @param charmap Mapping from chars to states (e.g., corax_map_nt for DNA)
 *  @param weights Alignment site weights, NULL=equal weights
 *  @param stats_mask Statistics to be computed, any combination of:
 *      CORAX_MSA_STATS_DUP_TAXA   duplicate taxon names
 *      CORAX_MSA_STATS_DUP_SEQS   duplicate/identical sequences
 *      CORAX_MSA_STATS_GAP_PROP   proportion of gaps
 *      CORAX_MSA_STATS_GAP_SEQS   fully undetermined sequences (=all-gap rows)
 *      CORAX_MSA_STATS_GAP_COLS   fully undetermined sites (=all-gap columns)
 *      CORAX_MSA_STATS_INV_PROP   proportion of invariant sites
 *      CORAX_MSA_STATS_INV_COLS   invariant columns
 *      CORAX_MSA_STATS_FREQS      state frequencies (NB: gaps are ignored!)
 *      CORAX_MSA_STATS_ALL        all of the above
 * */
CORAX_EXPORT corax_msa_stats_t *
corax_msa_compute_stats(const corax_msa_t   *msa,
                        unsigned int         states,
                        const corax_state_t *tipmap,
                        const unsigned int  *weights,
                        unsigned long        stats_mask)
{
  if (!msa)
  {
    corax_set_error(CORAX_ERROR_INVALID_PARAM, "MSA structure is NULL");
    return CORAX_FAILURE;
  }

  if (!tipmap)
  {
    corax_set_error(CORAX_ERROR_INVALID_PARAM,
                    "Character-to-state mapping (charmap) is NULL");
    return CORAX_FAILURE;
  }

  corax_msa_stats_t *stats =
      (corax_msa_stats_t *)calloc(1, sizeof(corax_msa_stats_t));

  if (!stats)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Cannot allocate memory for MSA statistics");
    return NULL;
  }

  const unsigned int  msa_count  = msa->count;
  const unsigned long msa_length = (unsigned long)msa->length;

  unsigned long  i, j, k;
  unsigned long  sum_weights     = 0;
  unsigned long  total_gap_count = 0;
  unsigned long *col_gap_weight  = NULL;
  unsigned long *seq_gap_weight  = NULL;
  size_t        *pair_rates      = NULL;
  size_t        *col_state_freq  = NULL;

  corax_state_t *inv_state  = NULL;
  unsigned long  inv_weight = 0;

  corax_state_t gap_state = 0;

  double *column_entropies = NULL;
  double  entropy          = 0;

  /* gap state has always all bits set to one */
  for (i = 0; i < states; ++i)
  {
    gap_state <<= 1;
    gap_state |= 1;
  }

  stats->states         = states;
  stats->gap_cols_count = 0;
  stats->gap_seqs_count = 0;
  stats->inv_cols_count = 0;

  /* search for duplicate taxa names (=sequence labels) */
  if (stats_mask & CORAX_MSA_STATS_DUP_TAXA)
  {
#ifdef USE_POSIX_SEARCH
    int retval = find_duplicate_strings_htable(msa->label,
                                               msa_count,
                                               &stats->dup_taxa_pairs,
                                               &stats->dup_taxa_pairs_count);
#else
    int retval = find_duplicate_strings(msa->label,
                                        msa_count,
                                        0,
                                        &stats->dup_taxa_pairs,
                                        &stats->dup_taxa_pairs_count);
#endif

    if (!retval)
    {
      corax_set_error(CORAX_ERROR_MEM_ALLOC, "Error finding duplicated taxa");
      goto error_exit;
    }
  }

  /* search for duplicate sequences */
  if (stats_mask & CORAX_MSA_STATS_DUP_SEQS)
  {
    int retval = find_duplicate_strings(msa->sequence,
                                        msa_count,
                                        msa_length,
                                        &stats->dup_seqs_pairs,
                                        &stats->dup_seqs_pairs_count);
    if (!retval)
    {
      corax_set_error(CORAX_ERROR_MEM_ALLOC,
                      "Error finding duplicated sequences");
      goto error_exit;
    }
  }

  /* compute empirical substitution rates */
  if (stats_mask & CORAX_MSA_STATS_SUBST_RATES)
  {
    size_t n_subst_rates = CORAX_SUBST_RATE_COUNT(states);
    stats->subst_rates   = (double *)calloc(n_subst_rates, sizeof(double));
    pair_rates           = (size_t *)calloc(states * states, sizeof(size_t));
    col_state_freq       = (size_t *)calloc(states, sizeof(size_t));

    if (!pair_rates || !stats->subst_rates || !col_state_freq)
    {
      corax_set_error(CORAX_ERROR_MEM_ALLOC,
                      "Cannot allocate memory for MSA statistics");
      goto error_exit;
    }

    compute_pair_rates(states,
                       msa_count,
                       msa_length,
                       (unsigned char **)msa->sequence,
                       weights,
                       tipmap,
                       col_state_freq,
                       pair_rates);

    k                = 0;
    double last_rate = pair_rates[(states - 2) * states + states - 1];
    if (last_rate < 1e-7) last_rate = 1;
    for (i = 0; i < states - 1; i++)
    {
      for (j = i + 1; j < states; j++)
      {
        stats->subst_rates[k++] = pair_rates[i * states + j] / last_rate;
        if (stats->subst_rates[k - 1] < 0.01) stats->subst_rates[k - 1] = 0.01;
        if (stats->subst_rates[k - 1] > 50.0) stats->subst_rates[k - 1] = 50.0;
      }
    }
    stats->subst_rates[k - 1] = 1.0;

    assert(k == n_subst_rates);

    free(col_state_freq);
    free(pair_rates);
    col_state_freq = pair_rates = NULL;
  }

  /* if we were asked to find duplicates only, no need to loop over sites */
  if (!(stats_mask
        & ~(CORAX_MSA_STATS_DUP_SEQS | CORAX_MSA_STATS_DUP_TAXA
            | CORAX_MSA_STATS_SUBST_RATES)))
    return stats;

  if (stats_mask & CORAX_MSA_STATS_FREQS)
  {
    stats->freqs = (double *)calloc(states, sizeof(double));

    if (!stats->freqs)
    {
      corax_set_error(CORAX_ERROR_MEM_ALLOC,
                      "Cannot allocate memory for empirical frequencies");
      goto error_exit;
    }
  }

  if (stats_mask & CORAX_MSA_STATS_GAP_COLS)
  {
    col_gap_weight = (unsigned long *)calloc(msa_length, sizeof(unsigned long));
  }

  if (stats_mask & CORAX_MSA_STATS_GAP_SEQS)
  {
    seq_gap_weight = (unsigned long *)calloc(msa_count, sizeof(unsigned long));
  }

  if (stats_mask & (CORAX_MSA_STATS_INV_COLS | CORAX_MSA_STATS_INV_PROP))
  {
    inv_state = (corax_state_t *)calloc(msa_length, sizeof(corax_state_t));
    stats->inv_cols =
        (unsigned long *)calloc(msa_length, sizeof(unsigned long));
    if (!inv_state || !stats->inv_cols)
    {
      corax_set_error(CORAX_ERROR_MEM_ALLOC,
                      "Cannot allocate memory for computing invariant sites");
      goto error_exit;
    }

    /* initialize all elements to the gap state */
    for (i = 0; i < msa_length; ++i) inv_state[i] = gap_state;
  }

  /* check memory allocation */
  if (((stats_mask & CORAX_MSA_STATS_FREQS) && !stats->freqs)
      || ((stats_mask & CORAX_MSA_STATS_GAP_COLS) && !col_gap_weight)
      || ((stats_mask & CORAX_MSA_STATS_GAP_SEQS) && !seq_gap_weight))
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Cannot allocate memory for MSA statistics");
    goto error_exit;
  }

  for (i = 0; i < msa_count; ++i)
  {
    const char *seqchars = msa->sequence[i];
    for (j = 0; j < msa_length; ++j)
    {
      const corax_state_t state       = tipmap[(int)seqchars[j]];
      const unsigned int  site_states = CORAX_STATE_POPCNT(state);
      const int           is_gap      = state == gap_state ? 1 : 0;
      const unsigned int  w           = weights ? weights[j] : 1;

      if (!state)
      {
        if (seqchars[j] == -1)
        {
          /* most likely sequence was already encoded and the original character
             is unknown at this point */
          corax_set_error(CORAX_ERROR_MSA_MAP_INVALID,
                          "Unknown state in sequence %d",
                          i + 1);
        }
        else
        {
          corax_set_error(CORAX_ERROR_MSA_MAP_INVALID,
                          "Unknown state %c at sequence %d position %d",
                          seqchars[j],
                          i + 1,
                          j + 1);
        }
        goto error_exit;
      }

      /* compute sum of weights (=alignment length before pattern compression)*/
      if (i == 0) sum_weights += w;

      if (stats_mask & (CORAX_MSA_STATS_GAP_PROP | CORAX_MSA_STATS_FREQS))
      {
        if (is_gap) total_gap_count += w;
      }

      if (stats_mask & CORAX_MSA_STATS_GAP_COLS)
      {
        if (is_gap) col_gap_weight[j] += 1;
        if (i == msa_count - 1 && col_gap_weight[j] == msa_count)
          stats->gap_cols_count++;
      }

      if (stats_mask & CORAX_MSA_STATS_GAP_SEQS)
      {
        if (is_gap) seq_gap_weight[i] += w;
        if (j == msa_length - 1 && seq_gap_weight[i] == sum_weights)
          stats->gap_seqs_count++;
      }

      if (stats_mask & (CORAX_MSA_STATS_INV_COLS | CORAX_MSA_STATS_INV_PROP))
      {
        inv_state[j] &= state;
        if (i == msa_count - 1 && CORAX_STATE_POPCNT(inv_state[j]) == 1)
        {
          inv_weight += w;
          stats->inv_cols[stats->inv_cols_count++] = j;
        }
      }

      if (stats_mask & CORAX_MSA_STATS_FREQS)
      {
        /* ignore gap sites when computing the base freqs */
        if (!is_gap)
        {
          double state_prob = ((double)w) / site_states;
          for (k = 0; k < states; ++k)
          {
            if (state & (1ll << k)) stats->freqs[k] += state_prob;
          }
        }
      }

    } // site loop
  }   // sequence loop

  /* compute proportion of invariant sites */
  if (stats_mask & (CORAX_MSA_STATS_INV_COLS | CORAX_MSA_STATS_INV_PROP))
  {
    stats->inv_prop = ((double)inv_weight) / sum_weights;
    stats->inv_cols =
        realloc(stats->inv_cols, stats->inv_cols_count * sizeof(unsigned long));
    free(inv_state);
    inv_state = NULL;
  }

  const size_t total_chars = sum_weights * msa_count;

  /* normalize frequencies */
  if (stats_mask & CORAX_MSA_STATS_FREQS)
  {
    for (k = 0; k < states; ++k)
      stats->freqs[k] /= total_chars - total_gap_count;
  }

  /* compute proportion of gaps */
  if (stats_mask & CORAX_MSA_STATS_GAP_PROP)
    stats->gap_prop = ((double)total_gap_count) / total_chars;

  /* detect gap-only columns */
  if ((stats_mask & CORAX_MSA_STATS_GAP_COLS))
  {
    if (stats->gap_cols_count > 0)
    {
      stats->gap_cols =
          (unsigned long *)calloc(stats->gap_cols_count, sizeof(unsigned long));

      if (!stats->gap_cols)
      {
        corax_set_error(CORAX_ERROR_MEM_ALLOC,
                        "Cannot allocate memory for gap columns");
        goto error_exit;
      }

      unsigned long c = 0;
      for (j = 0; j < msa_length; ++j)
      {
        if (col_gap_weight[j] == msa_count) stats->gap_cols[c++] = j;
      }
      assert(c == stats->gap_cols_count);
    }
    free(col_gap_weight);
    col_gap_weight = NULL;
  }

  /* detect gap-only sequences */
  if ((stats_mask & CORAX_MSA_STATS_GAP_SEQS))
  {
    if (stats->gap_seqs_count > 0)
    {
      stats->gap_seqs =
          (unsigned long *)calloc(stats->gap_seqs_count, sizeof(unsigned long));

      if (!stats->gap_seqs)
      {
        corax_set_error(CORAX_ERROR_MEM_ALLOC,
                        "Cannot allocate memory for gap sequences");
        goto error_exit;
      }

      unsigned long c = 0;
      for (i = 0; i < msa_count; ++i)
      {
        if (seq_gap_weight[i] == sum_weights) stats->gap_seqs[c++] = i;
      }
      assert(c == stats->gap_seqs_count);
    }
    free(seq_gap_weight);
    seq_gap_weight = NULL;
  }

  /* compute MSA entropy */
  if (stats_mask & CORAX_MSA_STATS_ENTROPY)
  {
    column_entropies        = corax_msa_column_entropies(msa, states, tipmap);
    for (i = 0; i < msa_length; ++i)
    {
      entropy += weights ? weights[i] * column_entropies[i] : column_entropies[i];
    }
    entropy /= sum_weights;

    stats->column_entropies = column_entropies;
    stats->entropy          = entropy;
  }

  return stats;

error_exit:
  if (inv_state) free(inv_state);
  if (col_gap_weight) free(col_gap_weight);
  if (seq_gap_weight) free(seq_gap_weight);
  if (col_state_freq) free(col_state_freq);
  if (pair_rates) free(pair_rates);
  corax_msa_destroy_stats(stats);

  return NULL;
}

CORAX_EXPORT void corax_msa_destroy_stats(corax_msa_stats_t *stats)
{
  if (stats->dup_taxa_pairs) free(stats->dup_taxa_pairs);

  if (stats->dup_seqs_pairs) free(stats->dup_seqs_pairs);

  if (stats->gap_seqs) free(stats->gap_seqs);

  if (stats->gap_cols) free(stats->gap_cols);

  if (stats->inv_cols) free(stats->inv_cols);

  if (stats->freqs) free(stats->freqs);

  if (stats->subst_rates) free(stats->subst_rates);

  if (stats->column_entropies) free(stats->column_entropies);

  free(stats);
}

/**
 * Filter MSA by removing the specified sequences and/or columns
 *
 * @param msa Multiple Sequence Alignment
 * @param remove_seqs 0-based indices of sequences to remove
 * @param remove_seqs_count size of @param remove_seqs array
 * @param remove_cols 0-based indices of columns to remove
 * @param remove_cols_count size of @param remove_cols array
 * @param inplace create new MSA structure for the filtered alignment (0)
 *                or re-use the original one (1)
 */
CORAX_EXPORT corax_msa_t *corax_msa_filter(corax_msa_t   *msa,
                                           unsigned long *remove_seqs,
                                           unsigned long  remove_seqs_count,
                                           unsigned long *remove_cols,
                                           unsigned long  remove_cols_count,
                                           unsigned int   inplace)
{
  if (!msa)
  {
    corax_set_error(CORAX_ERROR_INVALID_PARAM, "MSA structure is NULL");
    return NULL;
  }

  if (remove_seqs_count && !remove_seqs)
  {
    corax_set_error(CORAX_ERROR_INVALID_PARAM,
                    "List of sequences to remove (remove_seqs) is NULL");
    return NULL;
  }

  if (remove_cols_count && !remove_cols)
  {
    corax_set_error(CORAX_ERROR_INVALID_PARAM,
                    "List of columns to remove (remove_seqs) is NULL");
    return NULL;
  }

  const unsigned long old_count  = (unsigned long)msa->count;
  const unsigned long old_length = (unsigned long)msa->length;

  unsigned long i, j;

  unsigned char *seqflag = NULL;
  unsigned char *colflag = NULL;
  corax_msa_t   *new_msa = NULL;

  if (remove_seqs_count)
  {
    seqflag = (unsigned char *)calloc(old_count, sizeof(unsigned char));
    for (i = 0; i < remove_seqs_count; ++i)
    {
      if (remove_seqs[i] < old_count)
        seqflag[remove_seqs[i]] = 1;
      else
      {
        corax_set_error(CORAX_ERROR_INVALID_PARAM,
                        "Invalid sequence number in remove list: %lu",
                        remove_seqs[i]);
        goto error_exit;
      }
    }
  }

  if (remove_cols_count)
  {
    colflag = (unsigned char *)calloc(old_length, sizeof(unsigned char));
    for (i = 0; i < remove_cols_count; ++i)
    {
      if (remove_cols[i] < old_length)
        colflag[remove_cols[i]] = 1;
      else
      {
        corax_set_error(CORAX_ERROR_INVALID_PARAM,
                        "Invalid column number in remove list: %lu",
                        remove_cols[i]);
        goto error_exit;
      }
    }
  }

  const unsigned long new_count  = old_count - remove_seqs_count;
  const unsigned long new_length = old_length - remove_cols_count;

  if (inplace)
    new_msa = msa;
  else
  {
    new_msa           = (corax_msa_t *)calloc(1, sizeof(corax_msa_t));
    new_msa->count    = (int)new_count;
    new_msa->length   = (int)new_length;
    new_msa->label    = (char **)calloc(new_count, sizeof(char *));
    new_msa->sequence = (char **)calloc(new_count, sizeof(char *));
    if (!msa->label || !new_msa->sequence) goto error_exit;
  }

  unsigned long seq_idx = 0;
  for (i = 0; i < old_count; ++i)
  {
    /* check if we should skip this sequence */
    if (seqflag && seqflag[i]) continue;

    if (inplace)
      new_msa->label[seq_idx] = msa->label[i];
    else
    {
      new_msa->label[seq_idx] = strdup(msa->label[i]);
      new_msa->sequence[seq_idx] =
          (char *)malloc((new_length + 1) * sizeof(char));
      if (!new_msa->label[seq_idx] || !new_msa->sequence[seq_idx])
        goto error_exit;
    }

    if (!colflag)
    {
      /* no columns to remove, just copy/assign the old sequence*/
      if (inplace)
        new_msa->sequence[seq_idx] = msa->sequence[i];
      else
        strcpy(new_msa->sequence[seq_idx], msa->sequence[i]);
    }
    else
    {
      unsigned long col_idx = 0;
      for (j = 0; j < old_length; ++j)
      {
        if (colflag[j]) continue;

        new_msa->sequence[seq_idx][col_idx++] = msa->sequence[i][j];
      }
      assert(col_idx == new_length);
      new_msa->sequence[seq_idx][new_length] = '\0';

      /* trim sequence to the new size */
      if (inplace)
      {
        new_msa->sequence[seq_idx] = (char *)realloc(
            new_msa->sequence[seq_idx], (new_length + 1) * sizeof(char));
      }
    }

    seq_idx++;
  }
  assert(seq_idx == new_count);

  if (inplace)
  {
    /* set new dimensions */
    new_msa->count  = (int)new_count;
    new_msa->length = (int)new_length;
    /* trim arrays to the new size */
    new_msa->sequence =
        (char **)realloc(new_msa->sequence, new_count * sizeof(char *));
    new_msa->label =
        (char **)realloc(new_msa->label, new_count * sizeof(char *));
  }

  if (seqflag) free(seqflag);
  if (colflag) free(colflag);

  return new_msa;

error_exit:
  if (seqflag) free(seqflag);
  if (colflag) free(colflag);
  if (new_msa)
  {
    for (i = 0; i < (unsigned long)new_msa->count; ++i)
    {
      if (msa->sequence && msa->sequence[i]) free(msa->sequence[i]);
      if (msa->label && msa->label[i]) free(msa->label[i]);
    }
    if (msa->sequence) free(msa->sequence);
    if (msa->label) free(msa->label);

    free(new_msa);
  }
  corax_set_error(CORAX_ERROR_MEM_ALLOC,
                  "Cannot allocate memory needed for MSA filtering");
  return NULL;
}

/**
 * Split MSA into several partitions (sub-MSAs)
 *
 * @param msa original MSA to be splitted
 * @param site_part array with 1-based partition indices for each column in
 * original MSA
 * @param part_count number of partitions
 *
 * @return list of per-partition MSA objects
 *
 *  Example:
 *
 *  corax_msa_t msa;
 *  msa.length = 5;
 *  // further init of msa ...
 *  unsigned int part_site = {1, 1, 2, 2, 1};
 *  corax_msa_t ** part_msa = corax_msa_split(msa, part_site, 2);
 *  // now part_msa contains 2 MSA objects:
 *  // part_msa[0] = 1st, 2nd and 5th columns of msa
 *  // part_msa[1] = 3rd and 4th columns of msa
 */
CORAX_EXPORT corax_msa_t **corax_msa_split(const corax_msa_t  *msa,
                                           const unsigned int *site_part,
                                           unsigned int        part_count)
{
  unsigned int  p;
  unsigned long i, j;

  corax_msa_t **part_msa_list =
      (corax_msa_t **)calloc(part_count, sizeof(corax_msa_t *));
  unsigned int *part_len =
      (unsigned int *)calloc(part_count, sizeof(unsigned int));

  if (!part_msa_list || !part_len) goto malloc_error;

  /* compute partition sizes */
  for (j = 0; j < (unsigned long)msa->length; j++)
  {
    if (site_part[j])
    {
      p = site_part[j] - 1;
      if (p < part_count)
        part_len[p]++;
      else
      {
        corax_set_error(
            CORAX_ERROR_INVALID_INDEX, "Partition index out of bounds: %u", p);
        goto error_exit;
      }
    }
  }

  for (p = 0; p < part_count; ++p)
  {
    part_msa_list[p] = (corax_msa_t *)calloc(1, sizeof(corax_msa_t));
    if (!part_msa_list[p]) goto malloc_error;

    part_msa_list[p]->count  = msa->count;
    part_msa_list[p]->length = 0;
    part_msa_list[p]->label  = NULL;
    part_msa_list[p]->sequence =
        (char **)calloc((size_t)msa->count, sizeof(char *));
    if (!part_msa_list[p]->sequence) goto malloc_error;

    for (i = 0; i < (unsigned long)msa->count; i++)
    {
      part_msa_list[p]->sequence[i] =
          (char *)malloc(part_len[p] * sizeof(char));
      if (!part_msa_list[p]->sequence[i]) goto malloc_error;
    }
  }

  for (j = 0; j < (unsigned long)msa->length; j++)
  {
    /* partition index of 0 means that site should be skipped */
    if (site_part[j])
    {
      p = site_part[j] - 1;
      assert(p < part_count);
      const unsigned int len = (unsigned int)part_msa_list[p]->length;
      part_msa_list[p]->length++;
      for (i = 0; i < (unsigned long)msa->count; i++)
        part_msa_list[p]->sequence[i][len] = msa->sequence[i][j];
    }
  }

  free(part_len);

  return part_msa_list;

malloc_error:
  corax_set_error(CORAX_ERROR_MEM_ALLOC,
                  "Cannot allocate memory needed for MSA splitting");
error_exit:
  if (part_len) free(part_len);
  if (part_msa_list)
  {
    for (p = 0; p < part_count; ++p)
    {
      if (part_msa_list[p])
      {
        if (part_msa_list[p]->sequence)
        {
          for (i = 0; i < (unsigned long)msa->count; i++)
          {
            if (part_msa_list[p]->sequence[i])
              free(part_msa_list[p]->sequence[i]);
          }
          free(part_msa_list[p]->sequence);
        }
        free(part_msa_list[p]);
      }
    }
    free(part_msa_list);
  }
  return NULL;
}
