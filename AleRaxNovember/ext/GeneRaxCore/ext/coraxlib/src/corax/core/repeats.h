#ifndef CORAX_CORE_REPEATS_H_
#define CORAX_CORE_REPEATS_H_

#include "corax/core/partition.h"

#define CORAX_REPEATS_LOOKUP_SIZE 2000000

/** @defgroup corax_repeats_t corax_repeats_t
 * Module relating to the repeats structure
 */

/**
 *  Site repeats is a technique that, for each node of a tree,
 *  compresses all the site CLVs that are expected to be equal,
 *  in order to save memory and computations.
 *
 *  Let u be a node, and let i and j be two sites. If the sites
 *  i and j are equal in all the sequences under the node u, then
 *  their CLVs are also equal, and thus do not need to be computed/stored
 *  twice. We say that they belong to the same repeat class.
 *
 *  For a given node, the site repeats technique identifies all the
 *  different repeat classes, and associates to each of them a unique
 *  class identifier (starting from 1 for each node), which is required
 *  to find the location of the CLV of this repeat class
 *
 *  Note that, if scaling is enabled, scalers are also (similarly to
 *  clvs) compressed by site repeats
 *
 *  @ingroup corax_repeats_t
 */
typedef struct corax_repeats
{
  /**
   * `pernode_site_id[u->clv_index][i]` is equal to the class identifier of the
   * site `i` at node `u` (all the sites that have the same class identifier
   * under the same node belong to the same repeat class)
   */
  unsigned int **pernode_site_id;

  /**
   * `pernode_id_site[u->clv_index][id]` is the first site which belongs to the
   * repeat class with identifier `id` at node `u`.
   */
  unsigned int **pernode_id_site;

  /**
   * `pernode_ids[u->clv_index]` is the number of different repeat classes (and
   * thus of different class identifiers) at node `u`
   */
  unsigned int *pernode_ids;

  /**
   * `perscale_ids[scaler_index]` is the number of different repeat classes for
   * the scaler associated with the `id` scaler_index.  For a given
   * `corax_operation_t *op`, if scalers are enabled:
   * `pernode_ids[op->parent_clv_index] ==
   * perscale_ids[op->parent_scaler_index]`
   */
  unsigned int *perscale_ids;

  /** `pernode_allocated_clvs[u->clv_index] * size_of_a_site_clv` is the total
   * size of the vector that was allocated for the (compressed) CLV of the node
   * `u`, where `size_of_a_site_clv` is the size of CLV chunk corresponding to
   * one site and one node. Note that this size might be larger than required
   * (its similar to the capacity of an STL container that can be larger than
   * its size)
   */
  unsigned int *pernode_allocated_clvs;

  /**
   * Returns `true` if site repeats compression should be applied on the parent
   * node of `left_clv` and `right_clv`. In particular, applying site repeats on
   * nodes that are close to the (virtual) root of the tree is often
   * counterproductive.  This function can be redefined, and its default
   * definition is is `corax_default_enable_repeats`
   */
  unsigned int (*enable_repeats)(struct corax_partition *partition,
                                 unsigned int            left_clv,
                                 unsigned int            right_clv);

  /**
   * callback called when repeats are updated. Reallocate the CLV
   * and scaler vector for the node whose clv_index is parent.
   * sites_to_alloc indicates the number of "unique sites",
   * (or rather class identifiers) for this node.
   *
   * This function can be redefined and its default definition
   * is corax_default_reallocate_repeats
   *
   * By default, we always reallocate the exact required size, in
   * order to save memory. An alternative strategy could consist
   * in preallocating the maximum size (assuming that there is no
   * repeat) once at the first call, and then not doing anything
   * at the next calls (to avoid deallocation/reallocation).
   */
  void (*reallocate_repeats)(struct corax_partition *partition,
                             unsigned int            parent,
                             int                     scaler_index,
                             unsigned int            sites_to_alloc);

  /*
   * The `lookup_buffer` corresponds to the "matrix M" in the original site
   * repeats publication. It is used to compute the repeat classes at a given
   * node `u`. Let `v` and `w` be the children of `u`, and let `nv` and `nw` be
   * their respective numbers of class repeats. If `nv*nw > lookup_buffer_size`,
   * site repeats cannot be computed, and we disable site repeats for node `u`.
   * Thus, `lookup_buffer_size` should be large enough to allow as much nodes as
   * possible to benefit from site repeats, without costing too much memory
   * (remember, there might be many partitions...). Default size is
   * `CORAX_REPEATS_LOOKUP_SIZE` and can be changed with
   * `corax_resize_repeats_lookup`
   */
  unsigned int *lookup_buffer;
  unsigned int  lookup_buffer_size;

  /**
   * Map each character (representing a state) to a unique identifier
   */
  char *charmap;

  /**
   * Those vectors are pre-allocated buffers for the site repeats algorithm (or
   * for using site repeats in some kernels functions). They are only relevant
   * in the scope of the function in which they are being used.
   */
  unsigned int *toclean_buffer;
  unsigned int *id_site_buffer;
  double *      bclv_buffer;

} corax_repeats_t;

#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in repeats.c */

  #define CORAX_GET_ID(site_id, site) ((site_id) ? ((site_id)[(site)]) : (site))
  #define CORAX_GET_SITE(id_site, site) ((id_site) ? ((id_site)[(site)]) : (site))

  CORAX_EXPORT int corax_repeats_enabled(const corax_partition_t *partition);

  CORAX_EXPORT void corax_resize_repeats_lookup(corax_partition_t *partition,
                                                unsigned int       size);

  CORAX_EXPORT unsigned int
  corax_get_sites_number(const corax_partition_t *partition,
                         unsigned int             clv_index);

  CORAX_EXPORT unsigned int *
  corax_get_site_id(const corax_partition_t *partition, unsigned int clv_index);

  CORAX_EXPORT unsigned int *
  corax_get_id_site(const corax_partition_t *partition, unsigned int clv_index);

  CORAX_EXPORT unsigned int
  corax_get_clv_size(const corax_partition_t *partition,
                     unsigned int             clv_index);

  CORAX_EXPORT unsigned int
  corax_default_enable_repeats(corax_partition_t *partition,
                               unsigned int       left_clv,
                               unsigned int       right_clv);

  CORAX_EXPORT void
  corax_default_reallocate_repeats(corax_partition_t *partition,
                                   unsigned int       parent,
                                   int                scaler_index,
                                   unsigned int       sites_to_alloc);

  CORAX_EXPORT int corax_repeats_initialize(corax_partition_t *partition);

  CORAX_EXPORT int corax_update_repeats_tips(corax_partition_t *  partition,
                                             unsigned int         tip_index,
                                             const corax_state_t *map,
                                             const char *         sequence);

  CORAX_EXPORT void corax_update_repeats(corax_partition_t *      partition,
                                         const corax_operation_t *op);

  CORAX_EXPORT void corax_disable_bclv(corax_partition_t *partition);

  CORAX_EXPORT void
  corax_fill_parent_scaler_repeats(unsigned int        sites,
                                   unsigned int *      parent_scaler,
                                   const unsigned int *psites,
                                   const unsigned int *left_scaler,
                                   const unsigned int *lids,
                                   const unsigned int *right_scaler,
                                   const unsigned int *rids);

  CORAX_EXPORT void
  corax_fill_parent_scaler_repeats_per_rate(unsigned int        sites,
                                            unsigned int        rates,
                                            unsigned int *      parent_scaler,
                                            const unsigned int *psites,
                                            const unsigned int *left_scaler,
                                            const unsigned int *lids,
                                            const unsigned int *right_scaler,
                                            const unsigned int *rids);

#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CORAX_CORE_REPEATS_H_ */
