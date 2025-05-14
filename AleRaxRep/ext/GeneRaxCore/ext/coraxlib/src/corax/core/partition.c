/*
    Copyright (C) 2015-2021 Tomas Flouri, Diego Darriba, Alexey Kozlov

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

#include "corax/corax.h"

__thread corax_hardware_t corax_hardware = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

static void dealloc_partition_data(corax_partition_t *partition);

static void dealloc_partition_data(corax_partition_t *partition)
{
  unsigned int i;

  if (!partition) return;

  free(partition->rates);
  free(partition->rate_weights);
  free(partition->eigen_decomp_valid);
  if (partition->prop_invar) free(partition->prop_invar);
  if (partition->invariant) free(partition->invariant);
  if (!partition->pattern_weights) free(partition->pattern_weights);

  if (partition->scale_buffer)
    for (i = 0; i < partition->scale_buffers; ++i)
      free(partition->scale_buffer[i]);
  free(partition->scale_buffer);

  if (partition->tipchars)
    for (i = 0; i < partition->tips; ++i)
      corax_aligned_free(partition->tipchars[i]);
  free(partition->tipchars);

  if (partition->ttlookup) corax_aligned_free(partition->ttlookup);

  if (partition->charmap) free(partition->charmap);

  if (partition->tipmap) free(partition->tipmap);

  if (partition->clv)
  {
    unsigned int start = (partition->attributes & CORAX_ATTRIB_PATTERN_TIP)
                             ? partition->tips
                             : 0;
    for (i = start; i < partition->clv_buffers + partition->tips; ++i)
      corax_aligned_free(partition->clv[i]);
  }
  free(partition->clv);

  if (partition->pmatrix)
  {
    // for (i = 0; i < partition->prob_matrices; ++i)
    corax_aligned_free(partition->pmatrix[0]);
  }
  free(partition->pmatrix);

  if (partition->subst_params)
    for (i = 0; i < partition->rate_matrices; ++i)
      corax_aligned_free(partition->subst_params[i]);
  free(partition->subst_params);

  if (partition->eigenvecs)
    for (i = 0; i < partition->rate_matrices; ++i)
      corax_aligned_free(partition->eigenvecs[i]);
  free(partition->eigenvecs);

  if (partition->inv_eigenvecs)
    for (i = 0; i < partition->rate_matrices; ++i)
      corax_aligned_free(partition->inv_eigenvecs[i]);
  free(partition->inv_eigenvecs);

  if (partition->eigenvals)
    for (i = 0; i < partition->rate_matrices; ++i)
      corax_aligned_free(partition->eigenvals[i]);
  free(partition->eigenvals);

  if (partition->frequencies)
    for (i = 0; i < partition->rate_matrices; ++i)
      corax_aligned_free(partition->frequencies[i]);
  free(partition->frequencies);

  if (partition->pattern_weights) free(partition->pattern_weights);

  if (partition->repeats)
  {
    corax_repeats_t *repeats = partition->repeats;
    for (i = 0; i < partition->nodes; ++i)
    {
      free(repeats->pernode_site_id[i]);
      free(repeats->pernode_id_site[i]);
    }
    free(repeats->pernode_site_id);
    free(repeats->pernode_id_site);
    free(repeats->pernode_ids);
    free(repeats->perscale_ids);
    free(repeats->pernode_allocated_clvs);
    free(repeats->lookup_buffer);
    free(repeats->toclean_buffer);
    free(repeats->id_site_buffer);
    free(repeats->bclv_buffer);
    free(repeats->charmap);
    free(repeats);
  }

  free(partition);
}

static int update_charmap(corax_partition_t *  partition,
                          const corax_state_t *map)
{
  unsigned int  i, j, k;
  unsigned int  new_states_count = 0;
  corax_state_t mapcopy[CORAX_ASCII_SIZE];

  memcpy(mapcopy, map, CORAX_ASCII_SIZE * sizeof(corax_state_t));

  /* find maximum value in charmap table */
  k = 0;
  while (partition->tipmap[k]) ++k;

  /* compute the number of new states in the map */
  for (i = 0; i < CORAX_ASCII_SIZE; ++i)
  {
    if (mapcopy[i])
    {
      /* check whether state map[i] already exists in the tipmap */
      for (j = 0; j < k; ++j)
        if (mapcopy[i] == partition->tipmap[j]) break;

      /* if it does not exist */
      if (j == k)
      {
        /* check whether it is the first time we find it */
        for (j = 0; j < i; ++j)
          if (mapcopy[j] == mapcopy[i]) break;

        /* if first time, increase number of new states */
        if (j == i) new_states_count++;
      }
    }
  }

  /* erase old charmap */
  memset(partition->charmap, 0, CORAX_ASCII_SIZE * sizeof(unsigned char));

  /* using this map we will have more than 256 states, so return an error */
  if (new_states_count + k >= CORAX_ASCII_SIZE)
  {
    corax_set_error(
        CORAX_ERROR_MSA_MAP_INVALID,
        "Cannot specify 256 or more states with CORAX_ATTRIB_PATTERN_TIP.");
    return CORAX_FAILURE;
  }

  /* traverse the new map */
  for (i = 0; i < CORAX_ASCII_SIZE; ++i)
  {
    if (mapcopy[i])
    {
      unsigned int code = 0;

      /* check whether state map[i] already exists in the tipmap */
      for (j = 0; j < k; ++j)
        if (mapcopy[i] == partition->tipmap[j]) break;

      if (j == k)
      {
        /* if it does not exist */
        code                    = k;
        partition->tipmap[code] = mapcopy[i];
        ++k;
      }
      else
      {
        /* if it exists already then j is the index to the tipmap */
        code = j;
      }

      partition->charmap[i] = (unsigned char)code;

      /* find all characters with the same state in the map */
      for (j = i + 1; j < CORAX_ASCII_SIZE; ++j)
      {
        if (mapcopy[i] == mapcopy[j])
        {
          partition->charmap[j] = (unsigned char)code;
          mapcopy[j]            = 0;
        }
      }
    }
  }

  /* set maximum number of states (including ambiguities), its logarithm,
     and the logarithm of states */
  if (new_states_count)
  {
    /* special cases which do not use remapping */
    if (partition->states == 4)
    {
      for (k = 0, i = 0; partition->tipmap[i]; ++i)
      {
        if (partition->tipmap[i] > k) k = (unsigned int)partition->tipmap[i];
      }

      partition->maxstates = k + 1;
    }
    else
      partition->maxstates += new_states_count;

    unsigned int l2_maxstates = (unsigned int)ceil(log2(partition->maxstates));

    /* allocate space for the precomputed tip-tip likelihood vector */
    size_t alloc_size = (1 << (2 * l2_maxstates))
                        * (partition->states_padded * partition->rate_cats);

    /* for AVX we do not need to reallocate ttlookup as it has fixed size */
    if ((partition->states == 4)
        && (partition->attributes & CORAX_ATTRIB_ARCH_AVX)
        && CORAX_HAS_CPU_FEATURE(avx_present))
      return CORAX_SUCCESS;

    free(partition->ttlookup);
    partition->ttlookup =
        corax_aligned_alloc(alloc_size * sizeof(double), partition->alignment);
    if (!partition->ttlookup)
    {
      corax_set_error(
          CORAX_ERROR_MEM_ALLOC,
          "Cannot allocate space for storing precomputed tip-tip CLVs.");
      return CORAX_FAILURE;
    }
  }

  return CORAX_SUCCESS;
}

/* create a bijective mapping from states to the range <1,maxstates> where
   maxstates is the maximum number of states (including ambiguities). It is
   needed to index the precomputed conditional likelihoods for each pair of
   states. The sequences are then encoded using this charmap, and we store
   the precomputated CLV for a charmapped pair i and j, at index:

   (i << ceil(log(maxstate)) + j) << log(states) << log(rates) */
static int create_charmap(corax_partition_t *  partition,
                          const corax_state_t *usermap)
{
  unsigned int  i, j, k = 0;
  corax_state_t m = 0;
  corax_state_t map[CORAX_ASCII_SIZE];

  /* If ascertainment bias correction attribute is set, CLVs will be allocated
     with additional sites for each state */
  unsigned int sites_alloc = partition->asc_bias_alloc
                                 ? partition->sites + partition->states
                                 : partition->sites;

  // memcpy(map, partition->map, CORAX_ASCII_SIZE * sizeof(unsigned int));
  memcpy(map, usermap, CORAX_ASCII_SIZE * sizeof(corax_state_t));

  if (!(partition->charmap =
            (unsigned char *)calloc(CORAX_ASCII_SIZE, sizeof(unsigned char))))
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Cannot allocate charmap for tip-tip precomputation.");
    return CORAX_FAILURE;
  }

  if (!(partition->tipmap =
            (corax_state_t *)calloc(CORAX_ASCII_SIZE, sizeof(corax_state_t))))
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Cannot allocate tipmap for tip-tip precomputation.");
    return CORAX_FAILURE;
  }

  /* create charmap (remapped table of ASCII characters to range 0,|states|)
     and tipmap which is a (1,|states|) -> state */
  for (i = 0; i < CORAX_ASCII_SIZE; ++i)
  {
    if (map[i])
    {
      if (map[i] > m) m = map[i];

      partition->charmap[i] = (unsigned char)k;
      partition->tipmap[k]  = map[i];
      for (j = i + 1; j < CORAX_ASCII_SIZE; ++j)
      {
        if (map[i] == map[j])
        {
          partition->charmap[j] = (unsigned char)k;
          map[j]                = 0;
        }
      }
      ++k;
    }
  }

  /* For all state settings for which remapping will not be done, we need to
     increment maxstates by one to account for a fictive state 0 which will
     never be used */
  if (partition->states == 4) k = (unsigned int)m + 1;

  /* set maximum number of states (including ambiguities), its logarithm,
     and the logarithm of states */
  partition->maxstates = k;

  unsigned int l2_maxstates = (unsigned int)ceil(log2(partition->maxstates));

  /* allocate space for the precomputed tip-tip likelihood vector */
  size_t alloc_size = (1 << (2 * l2_maxstates))
                      * (partition->states_padded * partition->rate_cats);

  /* dedicated 4x4 function  - if AVX is not used we can allocate less space
     in case not all 16 possible ambiguities are present */
  if ((partition->states == 4)
      && (partition->attributes & CORAX_ATTRIB_ARCH_AVX)
      && CORAX_HAS_CPU_FEATURE(avx_present))
  {
    partition->ttlookup = corax_aligned_alloc(
        1024 * partition->rate_cats * sizeof(double), partition->alignment);
    if (!partition->ttlookup)
    {
      corax_set_error(
          CORAX_ERROR_MEM_ALLOC,
          "Cannot allocate space for storing precomputed tip-tip CLVs.");
      return CORAX_FAILURE;
    }
  }
  else
  {
    partition->ttlookup =
        corax_aligned_alloc(alloc_size * sizeof(double), partition->alignment);
    if (!partition->ttlookup)
    {
      corax_set_error(
          CORAX_ERROR_MEM_ALLOC,
          "Cannot allocate space for storing precomputed tip-tip CLVs.");
      return CORAX_FAILURE;
    }
  }

  /* allocate tip character arrays */
  partition->tipchars =
      (unsigned char **)calloc(partition->tips, sizeof(unsigned char *));
  if (!partition->tipchars)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Cannot allocate space for storing tip characters.");
    return CORAX_FAILURE;
  }

  for (i = 0; i < partition->tips; ++i)
  {
    partition->tipchars[i] =
        (unsigned char *)malloc(sites_alloc * sizeof(unsigned char));
    if (!partition->tipchars[i])
    {
      corax_set_error(CORAX_ERROR_MEM_ALLOC,
                      "Cannot allocate space for storing tip characters.");
      return CORAX_FAILURE;
    }
  }

  return CORAX_SUCCESS;
}

CORAX_EXPORT corax_partition_t *
             corax_partition_create(unsigned int tips,
                                    unsigned int clv_buffers,
                                    unsigned int states,
                                    unsigned int sites,
                                    unsigned int rate_matrices,
                                    unsigned int prob_matrices,
                                    unsigned int rate_cats,
                                    unsigned int scale_buffers,
                                    unsigned int attributes)
{
  unsigned int i;
  unsigned int sites_alloc;

  // Multiple specified ARCHs is no problem, we will autoselect the best available later.
  // "best" is defined as AVX2 > AVX > SSE3
  //if (CORAX_POPCNT32(attributes & CORAX_ATTRIB_ARCH_MASK) > 1)
  //{
  //  corax_set_error(CORAX_ERROR_INVALID_PARAM,
  //                  "Multiple architecture flags specified.");
  //  return CORAX_FAILURE;
  //}

  /* disable repeats if there are to few sites */
  if (sites < 16 && (attributes & CORAX_ATTRIB_SITE_REPEATS))
  {
    attributes &= ~CORAX_ATTRIB_SITE_REPEATS;
  }

  /* allocate partition */
  corax_partition_t *partition =
      (corax_partition_t *)malloc(sizeof(corax_partition_t));
  if (!partition)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Cannot allocate memory for partition.");
    return CORAX_FAILURE;
  }

  /* extract architecture and set vectorization parameters */
  partition->alignment     = CORAX_ALIGNMENT_CPU;
  partition->attributes    = attributes;
  partition->states_padded = states;
#ifdef HAVE_SSE3
  if (attributes & CORAX_ATTRIB_ARCH_SSE && CORAX_HAS_CPU_FEATURE(sse3_present))
  {
    partition->alignment     = CORAX_ALIGNMENT_SSE;
    partition->states_padded = (states + 1) & 0xFFFFFFFE;
  }
#endif
#ifdef HAVE_AVX
  if (attributes & CORAX_ATTRIB_ARCH_AVX && CORAX_HAS_CPU_FEATURE(avx_present))
  {
    partition->alignment     = CORAX_ALIGNMENT_AVX;
    partition->states_padded = (states + 3) & 0xFFFFFFFC;
  }
#endif
#ifdef HAVE_AVX2
  if (attributes & CORAX_ATTRIB_ARCH_AVX2 && CORAX_HAS_CPU_FEATURE(avx2_present))
  {
    partition->alignment     = CORAX_ALIGNMENT_AVX;
    partition->states_padded = (states + 3) & 0xFFFFFFFC;
  }
#endif

  unsigned int states_padded = partition->states_padded;

  /* initialize properties */

  partition->tips               = tips;
  partition->clv_buffers        = clv_buffers;
  partition->nodes              = tips + clv_buffers;
  partition->states             = states;
  partition->sites              = sites;
  partition->pattern_weight_sum = sites;

  partition->rate_matrices = rate_matrices;
  partition->prob_matrices = prob_matrices;
  partition->rate_cats     = rate_cats;
  partition->scale_buffers = scale_buffers;

  partition->prop_invar      = NULL;
  partition->invariant       = NULL;
  partition->pattern_weights = NULL;

  partition->eigenvecs     = NULL;
  partition->inv_eigenvecs = NULL;
  partition->eigenvals     = NULL;

  partition->rates              = NULL;
  partition->rate_weights       = NULL;
  partition->subst_params       = NULL;
  partition->scale_buffer       = NULL;
  partition->frequencies        = NULL;
  partition->eigen_decomp_valid = 0;

  partition->ttlookup = NULL;
  partition->tipchars = NULL;
  partition->charmap  = NULL;
  partition->tipmap   = NULL;

  partition->repeats = NULL;

  /* If ascertainment bias correction attribute is set, CLVs will be allocated
     with additional sites for each state */
  partition->asc_bias_alloc =
      (partition->attributes & (CORAX_ATTRIB_AB_MASK | CORAX_ATTRIB_AB_FLAG))
      > 0;
  partition->asc_additional_sites = (partition->asc_bias_alloc ? states : 0);
  sites_alloc = (unsigned int)partition->asc_additional_sites + sites;

  /* allocate structures */

  /* eigen_decomp_valid */
  partition->eigen_decomp_valid =
      (int *)calloc(partition->rate_matrices, sizeof(int));
  if (!partition->eigen_decomp_valid)
  {
    dealloc_partition_data(partition);
    corax_errno = CORAX_ERROR_MEM_ALLOC;
    corax_set_error(CORAX_ERROR_MEM_ALLOC, "Unable to allocate enough memory.");
    return CORAX_FAILURE;
  }
  /* clv */
  partition->clv = (double **)calloc(partition->nodes, sizeof(double *));
  if (!partition->clv)
  {
    dealloc_partition_data(partition);
    corax_errno = CORAX_ERROR_MEM_ALLOC;
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Unable to allocate enough memory for CLVs.");
    return CORAX_FAILURE;
  }

  /* if site repeats are enabled, we allocate CLVs dynamically */
  if (!corax_repeats_enabled(partition))
  {
    /* if tip pattern precomputation is enabled, then do not allocate CLV space
       for the tip nodes */
    unsigned int start = (partition->attributes & CORAX_ATTRIB_PATTERN_TIP)
                             ? partition->tips
                             : 0;

    for (i = start; i < partition->tips + partition->clv_buffers; ++i)
    {
      partition->clv[i] = corax_aligned_alloc(sites_alloc * states_padded
                                                  * rate_cats * sizeof(double),
                                              partition->alignment);
      if (!partition->clv[i])
      {
        dealloc_partition_data(partition);
        corax_errno = CORAX_ERROR_MEM_ALLOC;
        corax_set_error(CORAX_ERROR_MEM_ALLOC,
                        "Unable to allocate enough memory for CLVs.");
        return CORAX_FAILURE;
      }
      /* zero-out CLV vectors to avoid valgrind warnings when using odd number
         of states with vectorized code */
      memset(partition->clv[i],
             0,
             (size_t)sites_alloc * states_padded * rate_cats * sizeof(double));
    }
  }
  /* pmatrix */
  partition->pmatrix =
      (double **)calloc(partition->prob_matrices, sizeof(double *));
  if (!partition->pmatrix)
  {
    dealloc_partition_data(partition);
    corax_errno = CORAX_ERROR_MEM_ALLOC;
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Unable to allocate enough memory for p-matrix.");
    return CORAX_FAILURE;
  }

  /* allocate transition probability matrices in contiguous space, in order
     to save the 'displacement' amount of memory per matrix, which is
     required for updating CLVs when the number of states is not a multiple
     of states_padded. */
  size_t displacement =
      (states_padded - states) * (states_padded) * sizeof(double);
  partition->pmatrix[0] =
      corax_aligned_alloc(partition->prob_matrices * states * states_padded
                                  * rate_cats * sizeof(double)
                              + displacement,
                          partition->alignment);
  if (!partition->pmatrix[0])
  {
    dealloc_partition_data(partition);
    corax_errno = CORAX_ERROR_MEM_ALLOC;
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Unable to allocate enough memory for p-matrix.");
    return CORAX_FAILURE;
  }
  for (i = 1; i < partition->prob_matrices; ++i)
    partition->pmatrix[i] =
        partition->pmatrix[i - 1] + states * states_padded * rate_cats;

  /* zero-out p-matrices to avoid valgrind warnings when using odd number of
     states with vectorized code */
  memset(partition->pmatrix[0],
         0,
         partition->prob_matrices * states * states_padded * rate_cats
                 * sizeof(double)
             + displacement);

  /* eigenvecs */
  partition->eigenvecs =
      (double **)calloc(partition->rate_matrices, sizeof(double *));
  if (!partition->eigenvecs)
  {
    dealloc_partition_data(partition);
    corax_errno = CORAX_ERROR_MEM_ALLOC;
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Unable to allocate enough memory for eigenvectors.");
    return CORAX_FAILURE;
  }
  for (i = 0; i < partition->rate_matrices; ++i)
  {
    partition->eigenvecs[i] = corax_aligned_alloc(
        states * states_padded * sizeof(double), partition->alignment);
    if (!partition->eigenvecs[i])
    {
      dealloc_partition_data(partition);
      corax_errno = CORAX_ERROR_MEM_ALLOC;
      corax_set_error(CORAX_ERROR_MEM_ALLOC,
                      "Unable to allocate enough memory for eigenvectors.");
      return CORAX_FAILURE;
    }
    memset(partition->eigenvecs[i], 0, states * states_padded * sizeof(double));
    /* TODO: don't forget to add code for SSE/AVX */
  }

  /* inv_eigenvecs */
  partition->inv_eigenvecs =
      (double **)calloc(partition->rate_matrices, sizeof(double *));
  if (!partition->inv_eigenvecs)
  {
    dealloc_partition_data(partition);
    corax_errno = CORAX_ERROR_MEM_ALLOC;
    corax_set_error(
        CORAX_ERROR_MEM_ALLOC,
        "Unable to allocate enough memory for inverse eigenvectors.");
    return CORAX_FAILURE;
  }
  for (i = 0; i < partition->rate_matrices; ++i)
  {
    partition->inv_eigenvecs[i] = corax_aligned_alloc(
        states * states_padded * sizeof(double), partition->alignment);
    if (!partition->inv_eigenvecs[i])
    {
      dealloc_partition_data(partition);
      corax_set_error(
          CORAX_ERROR_MEM_ALLOC,
          "Unable to allocate enough memory for inverse eigenvectors.");
      return CORAX_FAILURE;
    }
    memset(partition->inv_eigenvecs[i],
           0,
           states * states_padded * sizeof(double));
    /* TODO: don't forget to add code for SSE/AVX */
  }

  /* eigenvals */
  partition->eigenvals =
      (double **)calloc(partition->rate_matrices, sizeof(double *));
  if (!partition->eigenvals)
  {
    dealloc_partition_data(partition);
    corax_errno = CORAX_ERROR_MEM_ALLOC;
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Unable to allocate enough memory for eigenvalues.");
    return CORAX_FAILURE;
  }
  for (i = 0; i < partition->rate_matrices; ++i)
  {
    partition->eigenvals[i] = corax_aligned_alloc(
        states_padded * sizeof(double), partition->alignment);
    if (!partition->eigenvals[i])
    {
      dealloc_partition_data(partition);
      corax_set_error(CORAX_ERROR_MEM_ALLOC,
                      "Unable to allocate enough memory for eigenvalues.");
      return CORAX_FAILURE;
    }
    memset(partition->eigenvals[i], 0, states_padded * sizeof(double));
    /* TODO: don't forget to add code for SSE/AVX */
  }

  /* subst_params */
  partition->subst_params =
      (double **)calloc(partition->rate_matrices, sizeof(double *));
  if (!partition->subst_params)
  {
    dealloc_partition_data(partition);
    corax_set_error(
        CORAX_ERROR_MEM_ALLOC,
        "Unable to allocate enough memory for substitution parameters.");
    return CORAX_FAILURE;
  }
  for (i = 0; i < partition->rate_matrices; ++i)
  {
    size_t subst_count = 0;
    if (partition->attributes & CORAX_ATTRIB_NONREV)
    {
      subst_count = states * states - states;
    }
    else
    {
      subst_count = (states * states - states) / 2;
    }
    partition->subst_params[i] = corax_aligned_alloc(subst_count *
        sizeof(double), partition->alignment);
    if (!partition->subst_params[i])
    {
      dealloc_partition_data(partition);
      corax_set_error(
          CORAX_ERROR_MEM_ALLOC,
          "Unable to allocate enough memory for substitution parameters.");
      return CORAX_FAILURE;
    }
    /* TODO: don't forget to add code for SSE/AVX */
  }

  /* frequencies */
  partition->frequencies =
      (double **)calloc(partition->rate_matrices, sizeof(double *));
  if (!partition->frequencies)
  {
    dealloc_partition_data(partition);
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Unable to allocate enough memory for frequencies.");
    return CORAX_FAILURE;
  }
  for (i = 0; i < partition->rate_matrices; ++i)
  {
    partition->frequencies[i] = corax_aligned_alloc(
        states_padded * sizeof(double), partition->alignment);
    if (!partition->frequencies[i])
    {
      dealloc_partition_data(partition);
      corax_set_error(CORAX_ERROR_MEM_ALLOC,
                      "Unable to allocate enough memory for frequencies.");
      return CORAX_FAILURE;
    }
    /* TODO: don't forget to add code for SSE/AVX */
    memset(partition->frequencies[i], 0, states_padded * sizeof(double));
  }

  /* rates */
  partition->rates = (double *)calloc(partition->rate_cats, sizeof(double));
  if (!partition->rates)
  {
    dealloc_partition_data(partition);
    corax_set_error(
        CORAX_ERROR_MEM_ALLOC,
        "Unable to allocate enough memory for heterogeneity rates.");
    return CORAX_FAILURE;
  }

  /* rate weights */
  partition->rate_weights =
      (double *)calloc(partition->rate_cats, sizeof(double));
  if (partition->rate_weights)
  {
    /* initialize to 1/n_rates */
    for (i = 0; i < partition->rate_cats; ++i)
      partition->rate_weights[i] = 1.0 / partition->rate_cats;
  }
  else
  {
    dealloc_partition_data(partition);
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Unable to allocate enough memory for rate weights.");
    return CORAX_FAILURE;
  }

  /* proportion of invariant sites */
  partition->prop_invar =
      (double *)calloc(partition->rate_matrices, sizeof(double));
  if (!partition->prop_invar)
  {
    dealloc_partition_data(partition);
    corax_set_error(
        CORAX_ERROR_MEM_ALLOC,
        "Unable to allocate enough memory for invar sites proportion.");
    return CORAX_FAILURE;
  }

  /* site weights */
  partition->pattern_weights =
      (unsigned int *)malloc(sites_alloc * sizeof(unsigned int));
  if (!partition->pattern_weights)
  {
    dealloc_partition_data(partition);
    corax_set_error(
        CORAX_ERROR_MEM_ALLOC,
        "Unable to allocate enough memory for site pattern weights.");
    return CORAX_FAILURE;
  }
  for (i = 0; i < partition->sites; ++i) partition->pattern_weights[i] = 1;
  /* additional positions if asc_bias is set are initialized to zero */
  for (i = sites; i < sites_alloc; ++i) partition->pattern_weights[i] = 0;

  /* scale_buffer */
  partition->scale_buffer =
      (unsigned int **)calloc(partition->scale_buffers, sizeof(unsigned int *));
  if (!partition->scale_buffer)
  {
    dealloc_partition_data(partition);
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Unable to allocate enough memory for scale buffers.");
    return CORAX_FAILURE;
  }
  /* if we use site repeats, we allocate scales dynamically (later) */
  if (!corax_repeats_enabled(partition))
  {
    for (i = 0; i < partition->scale_buffers; ++i)
    {
      size_t scaler_size = (attributes & CORAX_ATTRIB_RATE_SCALERS)
                               ? sites_alloc * rate_cats
                               : sites_alloc;
      partition->scale_buffer[i] =
          (unsigned int *)calloc(scaler_size, sizeof(unsigned int));
      if (!partition->scale_buffer[i])
      {
        dealloc_partition_data(partition);
        corax_set_error(CORAX_ERROR_MEM_ALLOC,
                        "Unable to allocate enough memory for scale buffers.");
        return CORAX_FAILURE;
      }
    }
  }

  if (corax_repeats_enabled(partition))
  {
    if (CORAX_FAILURE == corax_repeats_initialize(partition))
    {
      dealloc_partition_data(partition);
      return CORAX_FAILURE;
    }
  }
  return partition;
}

CORAX_EXPORT void corax_partition_destroy(corax_partition_t *partition)
{
  dealloc_partition_data(partition);
}

static int set_tipchars_4x4(corax_partition_t *  partition,
                            unsigned int         tip_index,
                            const corax_state_t *map,
                            const char *         sequence)
{
  corax_state_t c;
  unsigned int  i;

  /* iterate through sites */
  for (i = 0; i < partition->sites; ++i)
  {
    if ((c = map[(int)sequence[i]]) == 0)
    {
      corax_set_error(CORAX_ERROR_TIPDATA_ILLEGALSTATE,
                      "Illegal state code in tip \"%c\"",
                      sequence[i]);
      return CORAX_FAILURE;
    }

    /* store states as the remapped characters from charmap */
    partition->tipchars[tip_index][i] = (unsigned char)c;
  }

  /* if asc_bias is set, we initialize the additional positions */
  if (partition->asc_bias_alloc)
  {
    for (i = 0; i < partition->states; ++i)
    {
      partition->tipchars[tip_index][partition->sites + i] = (unsigned char)1
                                                             << i;
    }
  }

  /* tipmap is never used in the 4x4 case except create and update_charmap */

  return CORAX_SUCCESS;
}

static int set_tipchars(corax_partition_t *  partition,
                        unsigned int         tip_index,
                        const corax_state_t *map,
                        const char *         sequence)
{
  corax_state_t  c;
  unsigned int   i;
  unsigned char *tipchars = partition->tipchars[tip_index];

  /* iterate through sites */
  for (i = 0; i < partition->sites; ++i)
  {
    if ((c = map[(int)sequence[i]]) == 0)
    {
      corax_set_error(CORAX_ERROR_TIPDATA_ILLEGALSTATE,
                      "Illegal state code in tip \"%c\"",
                      sequence[i]);
      return CORAX_FAILURE;
    }

    /* store states as the remapped characters from charmap */
    tipchars[i] = partition->charmap[(int)(sequence[i])];
  }

  /* if asc_bias is set, we initialize the additional positions */
  if (partition->asc_bias_alloc)
  {
    tipchars += partition->sites;
    memset(tipchars, 0, partition->states);

    /* tip chars should go in the same order as expected, or the pattern
       weights for the invariant sites would not match the correct character.
       For example, the expected order of amino acids is A,R,N,..., and the
       tipchars order is 1,16,13,... (i.e., not sequential)  */
    for (i = 0; i < partition->maxstates; ++i)
    {
      corax_state_t state = partition->tipmap[i];
      if (CORAX_STATE_POPCNT(state) == 1)
      {
        unsigned int pos = CORAX_STATE_CTZ(state);
        assert(pos < partition->states);
        tipchars[pos] = i;
      }
    }
  }
  return CORAX_SUCCESS;
}

static int set_tipclv(corax_partition_t *  partition,
                      unsigned int         tip_index,
                      const corax_state_t *map,
                      const char *         sequence)
{
  corax_state_t c;
  unsigned int  i, j;
  double *      tipclv = partition->clv[tip_index];

  corax_repeats_t *repeats     = partition->repeats;
  int              use_repeats = corax_repeats_enabled(partition);
  unsigned int     ids =
      use_repeats ? repeats->pernode_ids[tip_index] : partition->sites;
  /* iterate through sites */
  for (i = 0; i < ids; ++i)
  {
    unsigned int index =
        use_repeats ? repeats->pernode_id_site[tip_index][i] : i;
    if ((c = map[(int)sequence[index]]) == 0)
    {
      corax_set_error(CORAX_ERROR_TIPDATA_ILLEGALSTATE,
                      "Illegal state code in tip \"%c\"",
                      sequence[index]);
      return CORAX_FAILURE;
    }

    /* decompose basecall into the encoded residues and set the appropriate
       positions in the tip vector */
    for (j = 0; j < partition->states; ++j)
    {
      tipclv[j] = c & 1;
      c >>= 1;
    }

    /* fill in the entries for the other gamma values */
    tipclv += partition->states_padded;
    for (j = 0; j < partition->rate_cats - 1; ++j)
    {
      memcpy(tipclv,
             tipclv - partition->states_padded,
             partition->states * sizeof(double));
      tipclv += partition->states_padded;
    }
  }

  /* if asc_bias is set, we initialize the additional positions */
  if (partition->asc_bias_alloc)
  {
    for (i = 0; i < partition->states; ++i)
    {
      for (j = 0; j < partition->states; ++j) { tipclv[j] = j == i; }

      /* fill in the entries for the other gamma values */
      tipclv += partition->states_padded;
      for (j = 0; j < partition->rate_cats - 1; ++j)
      {
        memcpy(tipclv,
               tipclv - partition->states_padded,
               partition->states * sizeof(double));
        tipclv += partition->states_padded;
      }
    }
  }

  return CORAX_SUCCESS;
}

CORAX_EXPORT int corax_set_tip_states(corax_partition_t *  partition,
                                      unsigned int         tip_index,
                                      const corax_state_t *map,
                                      const char *         sequence)
{
  int rc;

  if (corax_repeats_enabled(partition))
  {
    if (CORAX_FAILURE
        == corax_update_repeats_tips(partition, tip_index, map, sequence))
      return CORAX_FAILURE;
  }
  if (partition->attributes & CORAX_ATTRIB_PATTERN_TIP)
  {
    /* create (or update) character map for tip-tip precomputations */
    if (partition->tipchars) { update_charmap(partition, map); }
    else
    {
      if (!create_charmap(partition, map))
      {
        dealloc_partition_data(partition);
        return CORAX_FAILURE;
      }
    }

    if (partition->states == 4)
      rc = set_tipchars_4x4(partition, tip_index, map, sequence);
    else
      rc = set_tipchars(partition, tip_index, map, sequence);
  }
  else
    rc = set_tipclv(partition, tip_index, map, sequence);

  return rc;
}

// TODO: <DOC> We should account for padding before calling this function
CORAX_EXPORT int corax_set_tip_clv(corax_partition_t *partition,
                                   unsigned int       tip_index,
                                   const double *     clv,
                                   int                padding)
{
  unsigned int i, j, k;

  if (partition->attributes & CORAX_ATTRIB_PATTERN_TIP)
  {
    corax_set_error(
        CORAX_ERROR_TIPDATA_ILLEGALFUNCTION,
        "Cannot use corax_set_tip_clv with CORAX_ATTRIB_PATTERN_TIP.");
    return CORAX_FAILURE;
  }

  double *tipclv = partition->clv[tip_index];

  for (i = 0; i < partition->sites; ++i)
  {
    for (j = 0; j < partition->rate_cats; ++j)
    {
      memcpy(tipclv, clv, partition->states * sizeof(double));
      tipclv += partition->states_padded;
    }
    clv += padding ? partition->states_padded : partition->states;
  }

  /* if asc_bias is set, we initialize the additional positions */
  if (partition->asc_bias_alloc)
  {
    for (i = 0; i < partition->states; ++i)
    {
      for (j = 0; j < partition->rate_cats; ++j)
      {
        for (k = 0; k < partition->states; ++k)
        {
          tipclv[k] = k == i ? 1.0 : 0.0;
        }
        tipclv += partition->states_padded;
      }
    }
  }

  return CORAX_SUCCESS;
}

CORAX_EXPORT void corax_set_pattern_weights(corax_partition_t * partition,
                                            const unsigned int *pattern_weights)
{
  unsigned int i;
  memcpy(partition->pattern_weights,
         pattern_weights,
         sizeof(unsigned int) * partition->sites);

  /* recompute the sum of weights */
  partition->pattern_weight_sum = 0;
  for (i = 0; i < partition->sites; ++i)
    partition->pattern_weight_sum += pattern_weights[i];
}

CORAX_EXPORT void corax_set_frequencies(corax_partition_t *partition,
                                        unsigned int       freqs_index,
                                        const double *     frequencies)
{
  unsigned int i;
  double       sum = 0.;

  memcpy(partition->frequencies[freqs_index],
         frequencies,
         partition->states * sizeof(double));

  /* make sure frequencies sum up to 1.0 */
  for (i = 0; i < partition->states; ++i)
    sum += partition->frequencies[freqs_index][i];

  if (fabs(sum - 1.0) > CORAX_MISC_EPSILON)
  {
    for (i = 0; i < partition->states; ++i)
      partition->frequencies[freqs_index][i] /= sum;
  }

  partition->eigen_decomp_valid[freqs_index] = 0;
}

CORAX_EXPORT void corax_set_category_rates(corax_partition_t *partition,
                                           const double *     rates)
{
  memcpy(partition->rates, rates, partition->rate_cats * sizeof(double));
}

CORAX_EXPORT void corax_set_category_weights(corax_partition_t *partition,
                                             const double *     rate_weights)
{
  memcpy(partition->rate_weights,
         rate_weights,
         partition->rate_cats * sizeof(double));
}

CORAX_EXPORT void corax_set_subst_params(corax_partition_t *partition,
                                         unsigned int       params_index,
                                         const double *     params)
{
  unsigned int count = partition->attributes & CORAX_ATTRIB_NONREV
                           ? CORAX_SUBST_RATE_COUNT_NONREV(partition->states)
                           : CORAX_SUBST_RATE_COUNT(partition->states);

  memcpy(partition->subst_params[params_index], params, count * sizeof(double));
  partition->eigen_decomp_valid[params_index] = 0;

  /* NOTE: For protein models PLL/RAxML do a rate scaling by 10.0/max_rate */
}

CORAX_EXPORT int corax_set_asc_bias_type(corax_partition_t *partition,
                                         int                asc_bias_type)
{
  unsigned int i;
  int          prop_invar    = 0;
  int          asc_bias_attr = asc_bias_type & CORAX_ATTRIB_AB_MASK;

  /* If the partition was created **without** ascertainment bias correction,
     CLVs do not have space allocated for the invariant states, and setting
     ascertaiment bias will likely produce a segfault later. */
  if (!partition->asc_bias_alloc)
  {
    corax_set_error(
        CORAX_ERROR_AB_NOSUPPORT,
        "Partition was not created with ascertainment bias support");
    return CORAX_FAILURE;
  }

  /* check that there is no proportion of invariant sites */
  for (i = 0; i < partition->rate_matrices; ++i)
    prop_invar |= (partition->prop_invar[i] > 0);
  if (asc_bias_type != 0 && prop_invar)
  {
    corax_set_error(
        CORAX_ERROR_INVAR_INCOMPAT,
        "Invariant sites are not compatible with asc bias correction");
    return CORAX_FAILURE;
  }

  /* check that asc_bias_type is either 0 or a valid type */
  if (asc_bias_attr != asc_bias_type)
  {
    corax_set_error(CORAX_ERROR_AB_INVALIDMETHOD,
                    "Illegal ascertainment bias algorithm \"%d\"",
                    asc_bias_type);
    return CORAX_FAILURE;
  }

  /* reset current ascertainment bias type (if any) */
  partition->attributes &= (unsigned int)~CORAX_ATTRIB_AB_MASK;

  /* set new ascertainment bias type */
  partition->attributes |= (unsigned int)asc_bias_attr;

  return CORAX_SUCCESS;
}

CORAX_EXPORT void corax_set_asc_state_weights(corax_partition_t * partition,
                                              const unsigned int *state_weights)
{
  assert(partition->asc_bias_alloc);
  memcpy(partition->pattern_weights + partition->sites,
         state_weights,
         sizeof(unsigned int) * partition->states);
}

CORAX_EXPORT void corax_fill_parent_scaler(unsigned int        scaler_size,
                                           unsigned int *      parent_scaler,
                                           const unsigned int *left_scaler,
                                           const unsigned int *right_scaler)
{
  unsigned int i;

  if (!left_scaler && !right_scaler)
    memset(parent_scaler, 0, sizeof(unsigned int) * scaler_size);
  else if (left_scaler && right_scaler)
  {
    memcpy(parent_scaler, left_scaler, sizeof(unsigned int) * scaler_size);
    for (i = 0; i < scaler_size; ++i) parent_scaler[i] += right_scaler[i];
  }
  else
  {
    if (left_scaler)
      memcpy(parent_scaler, left_scaler, sizeof(unsigned int) * scaler_size);
    else
      memcpy(parent_scaler, right_scaler, sizeof(unsigned int) * scaler_size);
  }
}
