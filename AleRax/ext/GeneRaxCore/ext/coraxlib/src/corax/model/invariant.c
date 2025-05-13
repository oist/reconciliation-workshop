#include "corax/corax.h"

CORAX_EXPORT int corax_update_invariant_sites_proportion(
    corax_partition_t *partition, unsigned int params_index, double prop_invar)
{

  /* check that there is no ascertainment bias correction */
  if (prop_invar != 0.0 && (partition->attributes & CORAX_ATTRIB_AB_MASK))
  {
    corax_set_error(
        CORAX_ERROR_INVAR_INCOMPAT,
        "Invariant sites are not compatible with asc bias correction");
    return CORAX_FAILURE;
  }

  /* validate new invariant sites proportion */
  if (prop_invar < 0 || prop_invar >= 1)
  {
    corax_set_error(CORAX_ERROR_INVAR_PROPORTION,
                    "Invalid proportion of invariant sites (%f)",
                    prop_invar);
    return CORAX_FAILURE;
  }

  if (params_index > partition->rate_matrices)
  {
    corax_set_error(CORAX_ERROR_INVAR_PARAMINDEX,
                    "Invalid params index (%u)",
                    params_index);
    return CORAX_FAILURE;
  }

  if (prop_invar > 0.0 && !partition->invariant)
  {
    if (!corax_update_invariant_sites(partition))
    {
      corax_set_error(CORAX_ERROR_INVAR_NONEFOUND, "No invariant sites found");
      return CORAX_FAILURE;
    }
  }

  partition->prop_invar[params_index] = prop_invar;

  return CORAX_SUCCESS;
}

CORAX_EXPORT unsigned int
corax_count_invariant_sites(const corax_partition_t *partition,
                            unsigned int *           state_inv_count)
{
  unsigned int  i, j, k;
  unsigned int  invariant_count = 0;
  unsigned int  tips            = partition->tips;
  unsigned int  sites           = partition->sites;
  unsigned int  states          = partition->states;
  corax_state_t gap_state       = 0;
  corax_state_t cur_state;
  int *         invariant = partition->invariant;
  double *      tipclv;

  /* gap state has always all bits set to one */
  for (i = 0; i < states; ++i)
  {
    gap_state <<= 1;
    gap_state |= 1;
  }

  if (state_inv_count)
    memset(state_inv_count, 0, states * sizeof(unsigned int));

  if (invariant)
  {
    /* count the invariant sites for each state */
    for (i = 0; i < sites; ++i)
    {
      if (invariant[i] > -1)
      {
        cur_state = (corax_state_t)invariant[i];
        /* since the invariant sites array is generated in the library,
           it should not contain invalid values */
        assert(cur_state < states);

        /* increase the counter and per-state count */
        invariant_count += partition->pattern_weights[i];
        if (state_inv_count) state_inv_count[cur_state]++;
      }
    }
  }
  else
  {
    if (partition->attributes & CORAX_ATTRIB_PATTERN_TIP)
    {
      for (j = 0; j < sites; ++j)
      {
        cur_state = gap_state;
        for (i = 0; i < tips; ++i)
        {
          cur_state &= ((unsigned int)(partition->tipchars[i][j]));
          if (!cur_state) { break; }
        }
        if (CORAX_STATE_POPCNT(cur_state) == 1)
        {
          invariant_count += partition->pattern_weights[j];
          if (state_inv_count) state_inv_count[CORAX_STATE_CTZ(cur_state)]++;
        }
      }
    }
    else
    {
      /* warning: note that this operation traverses the clvs by columns, and
         hence it may be slow. If CORAX_ATTRIB_PATTERN_TIP is not set, I suggest
         to call corax_update_invariant_sites() before calling this function in
         order to populate partition->invariant beforehand. It can be freed
         afterwards. */
      unsigned int span_padded =
          partition->rate_cats * partition->states_padded;

      for (j = 0; j < sites; ++j)
      {
        unsigned int clv_shift = j * span_padded;
        tipclv                 = partition->clv[0] + clv_shift;
        corax_state_t state    = gap_state;
        for (i = 0; i < tips; ++i)
        {
          tipclv    = partition->clv[i] + clv_shift;
          cur_state = 0;
          for (k = 0; k < states; ++k)
          {
            cur_state |= ((corax_state_t)tipclv[k] << k);
          }
          state &= cur_state;
          if (!state) { break; }
        }
        if (CORAX_STATE_POPCNT(state) == 1)
        {
          invariant_count += partition->pattern_weights[j];
          if (state_inv_count) state_inv_count[CORAX_STATE_CTZ(state)]++;
        }
      }
    }
  }
  return invariant_count;
}

CORAX_EXPORT int corax_update_invariant_sites(corax_partition_t *partition)
{
  unsigned int   i, j, k;
  corax_state_t  state;
  unsigned int   states        = partition->states;
  unsigned int   states_padded = partition->states_padded;
  unsigned int   sites         = partition->sites;
  unsigned int   tips          = partition->tips;
  unsigned int   rate_cats     = partition->rate_cats;
  corax_state_t  gap_state     = 0;
  corax_state_t *invariant;
  double *       tipclv;

  /* gap state has always all bits set to one */
  for (i = 0; i < states; ++i)
  {
    gap_state <<= 1;
    gap_state |= 1;
  }

  /* allocate array (on first call) denoting the frequency index for invariant
     sites, or -1 for variant sites */
  if (!partition->invariant)
  {
    partition->invariant = (int *)malloc(sites * sizeof(int));
  }

  invariant = (corax_state_t *)malloc(sites * sizeof(corax_state_t));

  if (!invariant || !partition->invariant)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Cannot allocate charmap for invariant sites array.");
    return CORAX_FAILURE;
  }

  /* initialize all elements to the gap state */
  for (i = 0; i < partition->sites; ++i) invariant[i] = gap_state;

  /* depending on the attribute flag, fill each element of the invariant array
     with the bitwise AND of gap and all states in the corresponding site */
  if (partition->attributes & CORAX_ATTRIB_PATTERN_TIP)
  {
    if (states == 4)
    {
      for (i = 0; i < tips; ++i)
        for (j = 0; j < sites; ++j)
        {
          state = (unsigned int)(partition->tipchars[i][j]);
          invariant[j] &= state;
        }
    }
    else
    {
      for (i = 0; i < tips; ++i)
        for (j = 0; j < sites; ++j)
        {
          state = partition->tipmap[(int)(partition->tipchars[i][j])];
          invariant[j] &= state;
        }
    }
  }
  else
  {
    unsigned int span_padded = rate_cats * states_padded;
    for (i = 0; i < tips; ++i)
    {
      const unsigned int *site_id = NULL;
      if (partition->repeats && partition->repeats->pernode_ids[i])
      {
        site_id = partition->repeats->pernode_site_id[i];
      }
      for (j = 0; j < sites; ++j)
      {
        unsigned int site = site_id ? site_id[j] : j;
        tipclv            = partition->clv[i] + span_padded * site;
        state             = 0;
        for (k = 0; k < states; ++k)
        {
          state |= ((corax_state_t)tipclv[k] << k);
        }
        invariant[j] &= state;
      }
    }
  }

  /* if all basecalls at current site are the same and not degenerate set the
     index in invariant to the frequency index of the basecall, otherwise -1 */
  for (i = 0; i < partition->sites; ++i)
  {
    if (invariant[i] == 0 || CORAX_STATE_POPCNT(invariant[i]) > 1)
      partition->invariant[i] = -1;
    else
      partition->invariant[i] = CORAX_STATE_CTZ(invariant[i]);
  }

  free(invariant);

  return CORAX_SUCCESS;
}
