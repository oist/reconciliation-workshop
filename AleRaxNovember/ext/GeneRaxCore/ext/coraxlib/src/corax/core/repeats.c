/*
    Copyright (C) 2015 Tomas Flouri, Diego Darriba

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

    Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include "corax/corax.h"

const unsigned int EMPTY_ELEMENT = (unsigned int)-1;

// map in charmap each char to a unique char identifier, according to map
static void repeats_fill_charmap(const corax_state_t *map, char *charmap)
{
  unsigned int i, j;
  char         maxChar = 0;
  for (i = 0; i < CORAX_ASCII_SIZE; ++i)
  {
    for (j = 0; j < i; ++j)
    {
      if (map[i] == map[j])
      {
        charmap[i] = charmap[j];
        break;
      }
    }
    if (!charmap[i]) charmap[i] = ++maxChar;
  }
}

CORAX_EXPORT int corax_repeats_enabled(const corax_partition_t *partition)
{
  return CORAX_ATTRIB_SITE_REPEATS & partition->attributes;
}

CORAX_EXPORT void corax_resize_repeats_lookup(corax_partition_t *partition,
                                              unsigned int       size)
{
  if (!size) return;
  partition->repeats->lookup_buffer_size = size;
  free(partition->repeats->lookup_buffer);
  partition->repeats->lookup_buffer = malloc(size * sizeof(unsigned int));
  memset(partition->repeats->lookup_buffer,
         EMPTY_ELEMENT,
         partition->repeats->lookup_buffer_size * sizeof(unsigned int));
}

CORAX_EXPORT unsigned int
corax_get_sites_number(const corax_partition_t *partition,
                       unsigned int             clv_index)
{
  unsigned int sites = partition->attributes & CORAX_ATTRIB_SITE_REPEATS
                           ? partition->repeats->pernode_ids[clv_index]
                           : 0;
  sites              = sites ? sites : partition->sites;
  sites += partition->asc_bias_alloc ? partition->states : 0;
  return sites;
}

CORAX_EXPORT unsigned int corax_get_clv_size(const corax_partition_t *partition,
                                             unsigned int             clv_index)
{
  return corax_get_sites_number(partition, clv_index) * partition->states_padded
         * partition->rate_cats;
}

CORAX_EXPORT unsigned int *corax_get_site_id(const corax_partition_t *partition,
                                             unsigned int             clv_index)
{
  unsigned int *site_id = 0;
  if (corax_repeats_enabled(partition)
      && partition->repeats->pernode_ids[clv_index])
    site_id = partition->repeats->pernode_site_id[clv_index];
  return site_id;
}

CORAX_EXPORT unsigned int *corax_get_id_site(const corax_partition_t *partition,
                                             unsigned int             clv_index)
{
  unsigned int *id_site = 0;
  if (corax_repeats_enabled(partition)
      && partition->repeats->pernode_ids[clv_index])
    id_site = partition->repeats->pernode_id_site[clv_index];
  return id_site;
}

CORAX_EXPORT unsigned int corax_default_enable_repeats(
    corax_partition_t *partition, unsigned int left_clv, unsigned int right_clv)
{
  corax_repeats_t *  repeats = partition->repeats;
  unsigned long long min_size =
      (unsigned long long)repeats->pernode_ids[left_clv]
      * (unsigned long long)repeats->pernode_ids[right_clv];
  return !(!min_size
           || ((unsigned long long)repeats->lookup_buffer_size <= min_size)
           || (repeats->pernode_ids[left_clv] > (partition->sites / 2))
           || (repeats->pernode_ids[right_clv] > (partition->sites / 2)));
}

CORAX_EXPORT int corax_repeats_initialize(corax_partition_t *partition)
{
  unsigned int sites_alloc =
      (unsigned int)partition->asc_additional_sites + partition->sites;
  unsigned int i;
  partition->repeats = malloc(sizeof(corax_repeats_t));
  if (!partition->repeats)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Unable to allocate enough memory for repeats structure.");
    return CORAX_FAILURE;
  }
  memset(partition->repeats, 0, sizeof(corax_repeats_t));
  corax_repeats_t *repeats    = partition->repeats;
  repeats->enable_repeats     = corax_default_enable_repeats;
  repeats->reallocate_repeats = corax_default_reallocate_repeats;
  repeats->pernode_site_id = calloc(partition->nodes, sizeof(unsigned int *));
  repeats->pernode_id_site = calloc(partition->nodes, sizeof(unsigned int *));
  if (!repeats->pernode_site_id || !repeats->pernode_id_site)
  {
    corax_set_error(
        CORAX_ERROR_MEM_ALLOC,
        "Unable to allocate enough memory for repeats identifiers.");
    return CORAX_FAILURE;
  }
  for (i = 0; i < partition->nodes; ++i)
  {
    repeats->pernode_site_id[i] = calloc(sites_alloc, sizeof(unsigned int));
    repeats->pernode_id_site[i] = calloc(sites_alloc, sizeof(unsigned int));
    if (!repeats->pernode_site_id[i])
    {
      corax_set_error(
          CORAX_ERROR_MEM_ALLOC,
          "Unable to allocate enough memory for repeats identifiers.");
      return CORAX_FAILURE;
    }
  }
  repeats->pernode_ids = calloc(partition->nodes, sizeof(unsigned int));
  repeats->perscale_ids =
      calloc(partition->scale_buffers, sizeof(unsigned int));
  repeats->pernode_allocated_clvs =
      calloc(partition->nodes, sizeof(unsigned int));
  repeats->lookup_buffer      = 0;
  repeats->lookup_buffer_size = 0;
  repeats->toclean_buffer     = malloc(sites_alloc * sizeof(unsigned int));
  repeats->id_site_buffer     = malloc(sites_alloc * sizeof(unsigned int));
  repeats->bclv_buffer =
      corax_aligned_alloc(sites_alloc * partition->rate_cats
                              * partition->states_padded * sizeof(double),
                          partition->alignment);
  repeats->charmap = calloc(CORAX_ASCII_SIZE, sizeof(char));
  if (!(repeats->pernode_ids && repeats->pernode_allocated_clvs
        && repeats->bclv_buffer && repeats->toclean_buffer
        && repeats->id_site_buffer && repeats->charmap))
  {
    corax_set_error(
        CORAX_ERROR_MEM_ALLOC,
        "Unable to allocate enough memory for one of the repeats buffer.");
    return CORAX_FAILURE;
  }
  return CORAX_SUCCESS;
}

CORAX_EXPORT int corax_update_repeats_tips(corax_partition_t *  partition,
                                           unsigned int         tip_index,
                                           const corax_state_t *map,
                                           const char *         sequence)
{
  if (!partition->repeats->lookup_buffer)
    corax_resize_repeats_lookup(partition, CORAX_REPEATS_LOOKUP_SIZE);

  unsigned int     s;
  corax_repeats_t *repeats = partition->repeats;
  unsigned int **  id_site = repeats->pernode_id_site;
  unsigned int     additional_sites =
      partition->asc_bias_alloc ? partition->states : 0;

  repeats_fill_charmap(map, repeats->charmap);
  repeats->pernode_ids[tip_index] = 0;
  unsigned int curr_id            = 0;
  /* fill pernode_site_id */
  for (s = 0; s < partition->sites; ++s)
  {
    unsigned int index_lookup =
        (unsigned int)repeats->charmap[(int)sequence[s]];
    if (EMPTY_ELEMENT == repeats->lookup_buffer[index_lookup])
    {
      repeats->toclean_buffer[curr_id]     = index_lookup;
      repeats->id_site_buffer[curr_id]     = s;
      repeats->lookup_buffer[index_lookup] = curr_id++;
    }
    repeats->pernode_site_id[tip_index][s] =
        repeats->lookup_buffer[index_lookup];
  }
  unsigned int ids                = curr_id;
  repeats->pernode_ids[tip_index] = ids;
  free(id_site[tip_index]);
  id_site[tip_index] = malloc(sizeof(unsigned int) * (ids + additional_sites));
  for (s = 0; s < ids; ++s)
  {
    id_site[tip_index][s] = repeats->id_site_buffer[s];
    repeats->lookup_buffer[repeats->toclean_buffer[s]] = EMPTY_ELEMENT;
  }
  for (s = 0; s < additional_sites; ++s)
  {
    id_site[tip_index][ids + s] = partition->sites + s;
    repeats->pernode_site_id[tip_index][partition->sites + s] = ids + s;
  }
  unsigned int sizealloc = (ids + additional_sites) * partition->states_padded
                           * partition->rate_cats * sizeof(double);
  free(partition->clv[tip_index]);

  partition->clv[tip_index] =
      corax_aligned_alloc(sizealloc, partition->alignment);
  if (!partition->clv[tip_index])
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Unable to allocate enough memory for repeats structure.");
    return CORAX_FAILURE;
  }
  /* zero-out CLV vectors to avoid valgrind warnings when using odd number of
       states with vectorized code */
  repeats->pernode_allocated_clvs[tip_index] = ids;
  memset(partition->clv[tip_index], 0, sizealloc);
  return CORAX_SUCCESS;
}

CORAX_EXPORT void corax_default_reallocate_repeats(corax_partition_t *partition,
                                                   unsigned int       parent,
                                                   int          scaler_index,
                                                   unsigned int sites_to_alloc)
{
  corax_repeats_t *repeats = partition->repeats;
  if (sites_to_alloc == repeats->pernode_allocated_clvs[parent]) return;
  repeats->pernode_allocated_clvs[parent] = sites_to_alloc;
  unsigned int **id_site                  = repeats->pernode_id_site;
  // reallocate clvs
  corax_aligned_free(partition->clv[parent]);
  partition->clv[parent] =
      corax_aligned_alloc(sites_to_alloc * partition->states_padded
                              * partition->rate_cats * sizeof(double),
                          partition->alignment);

  if (!partition->clv[parent])
  {
    corax_errno = CORAX_ERROR_MEM_ALLOC;
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Unable to allocate enough memory for repeats structure.");
    return;
  }
  // reallocate scales
  if (CORAX_SCALE_BUFFER_NONE != scaler_index)
  {
    unsigned int scaler_size = sites_to_alloc;
    if (partition->attributes & CORAX_ATTRIB_RATE_SCALERS)
      scaler_size *= partition->rate_cats;
    free(partition->scale_buffer[scaler_index]);
    partition->scale_buffer[scaler_index] =
        calloc(scaler_size, sizeof(unsigned int));
  }
  // reallocate id to site lookup
  free(id_site[parent]);
  id_site[parent] = malloc(sites_to_alloc * sizeof(unsigned int));
  // avoid valgrind errors
  memset(partition->clv[parent], 0, sites_to_alloc);
}

/* Fill the repeat structure in partition for the parent node of op */
CORAX_EXPORT void corax_update_repeats(corax_partition_t *      partition,
                                       const corax_operation_t *op)
{
  if (!partition->repeats->lookup_buffer)
    corax_resize_repeats_lookup(partition, CORAX_REPEATS_LOOKUP_SIZE);

  corax_repeats_t *   repeats        = partition->repeats;
  unsigned int        left           = op->child1_clv_index;
  unsigned int        right          = op->child2_clv_index;
  unsigned int        parent         = op->parent_clv_index;
  unsigned int **     site_ids       = repeats->pernode_site_id;
  unsigned int *      site_id_parent = site_ids[parent];
  const unsigned int *site_id_left   = site_ids[left];
  const unsigned int *site_id_right  = site_ids[right];
  const unsigned int  ids_left       = repeats->pernode_ids[left];
  unsigned int **     id_site        = repeats->pernode_id_site;
  unsigned int *      toclean_buffer = repeats->toclean_buffer;
  unsigned int *      id_site_buffer = repeats->id_site_buffer;
  unsigned int        curr_id        = 0;
  unsigned int        additional_sites =
      partition->asc_bias_alloc ? partition->states : 0;
  unsigned int sites_to_alloc;
  unsigned int s;
  unsigned int ids = 0;
  // in case site repeats is activated but not used for this node
  if (!partition->repeats->enable_repeats(partition, left, right))
  {
    sites_to_alloc               = partition->sites + additional_sites;
    repeats->pernode_ids[parent] = 0;
    if (op->parent_scaler_index != CORAX_SCALE_BUFFER_NONE)
      repeats->perscale_ids[op->parent_scaler_index] = 0;
  }
  else
  {
    // fill the parent repeats identifiers
    for (s = 0; s < partition->sites; ++s)
    {
      unsigned int index_lookup = site_id_left[s] + site_id_right[s] * ids_left;
      unsigned int id           = repeats->lookup_buffer[index_lookup];
      if (EMPTY_ELEMENT == id)
      {
        toclean_buffer[curr_id]              = index_lookup;
        id_site_buffer[curr_id]              = s;
        id                                   = curr_id;
        repeats->lookup_buffer[index_lookup] = curr_id++;
      }
      site_id_parent[s] = id;
    }
    ids = curr_id;
    for (s = 0; s < additional_sites; ++s)
    {
      site_id_parent[s + partition->sites] = ids + s;
    }
    repeats->pernode_ids[parent] = ids;
    if (op->parent_scaler_index != CORAX_SCALE_BUFFER_NONE)
      repeats->perscale_ids[op->parent_scaler_index] = ids;
    sites_to_alloc = ids + additional_sites;
  }

  repeats->reallocate_repeats(
      partition, op->parent_clv_index, op->parent_scaler_index, sites_to_alloc);

  // there is no repeats. Set pernode_ids to 0
  // to force the core functions not to use repeats
  if (sites_to_alloc >= partition->sites + additional_sites)
  {
    repeats->pernode_ids[parent] = 0;
    if (op->parent_scaler_index != CORAX_SCALE_BUFFER_NONE)
      repeats->perscale_ids[op->parent_scaler_index] = 0;
  }

  // set id to site lookups
  for (s = 0; s < ids; ++s)
  {
    id_site[parent][s]                        = id_site_buffer[s];
    repeats->lookup_buffer[toclean_buffer[s]] = EMPTY_ELEMENT;
  }
  for (s = 0; s < additional_sites; ++s)
  {
    id_site[parent][s + ids] = partition->sites + s;
  }
}

CORAX_EXPORT void corax_disable_bclv(corax_partition_t *partition)
{
  if (!corax_repeats_enabled(partition)) return;
  corax_aligned_free(partition->repeats->bclv_buffer);
  partition->repeats->bclv_buffer = 0;
}

CORAX_EXPORT void
corax_fill_parent_scaler_repeats(unsigned int        sites,
                                 unsigned int *      parent_scaler,
                                 const unsigned int *psites,
                                 const unsigned int *left_scaler,
                                 const unsigned int *lids,
                                 const unsigned int *right_scaler,
                                 const unsigned int *rids)
{
  // no repeats
  if (!lids && !rids)
  {
    corax_fill_parent_scaler(sites, parent_scaler, left_scaler, right_scaler);
    return;
  }

  // no scalers
  if (!left_scaler && !right_scaler)
  {
    memset(parent_scaler, 0, sizeof(unsigned int) * sites);
    return;
  }

  unsigned int i;
  if (!psites)
  {
    memset(parent_scaler, 0, sizeof(unsigned int) * sites);
    if (left_scaler)
    {
      if (lids)
      {
        for (i = 0; i < sites; ++i) parent_scaler[i] += left_scaler[lids[i]];
      }
      else
      {
        for (i = 0; i < sites; ++i) parent_scaler[i] += left_scaler[i];
      }
    }
    if (right_scaler)
    {
      if (rids)
      {
        for (i = 0; i < sites; ++i) parent_scaler[i] += right_scaler[rids[i]];
      }
      else
      {
        for (i = 0; i < sites; ++i) parent_scaler[i] += right_scaler[i];
      }
    }
  }
  else
  {
    if (left_scaler && right_scaler)
    {
      for (i = 0; i < sites; ++i)
        parent_scaler[i] =
            left_scaler[lids[psites[i]]] + right_scaler[rids[psites[i]]];
    }
    else if (left_scaler)
    {
      for (i = 0; i < sites; ++i)
        parent_scaler[i] = left_scaler[lids[psites[i]]];
    }
    else
    {
      for (i = 0; i < sites; ++i)
        parent_scaler[i] = right_scaler[rids[psites[i]]];
    }
  }
}

CORAX_EXPORT void
corax_fill_parent_scaler_repeats_per_rate(unsigned int        sites,
                                          unsigned int        rates,
                                          unsigned int *      parent_scaler,
                                          const unsigned int *psites,
                                          const unsigned int *left_scaler,
                                          const unsigned int *lids,
                                          const unsigned int *right_scaler,
                                          const unsigned int *rids)
{
  unsigned int total_size     = sites * rates;
  unsigned int cpy_size       = rates * sizeof(unsigned int);
  unsigned int total_cpy_size = total_size * sizeof(unsigned int);
  // no repeats
  if (!lids && !rids)
  {
    corax_fill_parent_scaler(
        total_size, parent_scaler, left_scaler, right_scaler);
    return;
  }

  // no scalers
  if (!left_scaler && !right_scaler)
  {
    memset(parent_scaler, 0, total_cpy_size);
    return;
  }

  unsigned int i, j;
  if (!psites)
  {
    memset(parent_scaler, 0, total_cpy_size);
    if (left_scaler)
    {
      if (lids)
      {
        for (i = 0; i < sites; ++i)
          memcpy(&parent_scaler[i * rates],
                 &left_scaler[lids[i] * rates],
                 cpy_size);
      }
      else
      {
        memcpy(parent_scaler, left_scaler, total_cpy_size);
      }
    }
    if (right_scaler)
    {
      if (rids)
      {
        for (i = 0; i < sites; ++i)
          for (j = 0; j < rates; ++j)
            parent_scaler[i * rates + j] += right_scaler[rids[i] * rates + j];
      }
      else
      {
        for (i = 0; i < total_size; ++i) parent_scaler[i] += right_scaler[i];
      }
    }
  }
  else
  {
    if (left_scaler && right_scaler)
    {
      for (i = 0; i < sites; ++i)
        for (j = 0; j < rates; ++j)
          parent_scaler[i * rates + j] =
              left_scaler[lids[psites[i]] * rates + j]
              + right_scaler[rids[psites[i]] * rates + j];
    }
    else if (left_scaler)
    {
      for (i = 0; i < sites; ++i)
        memcpy(&parent_scaler[i * rates],
               &left_scaler[lids[psites[i]] * rates],
               cpy_size);
    }
    else
    {
      for (i = 0; i < sites; ++i)
        memcpy(&parent_scaler[i * rates],
               &right_scaler[rids[psites[i]] * rates],
               cpy_size);
    }
  }
}
