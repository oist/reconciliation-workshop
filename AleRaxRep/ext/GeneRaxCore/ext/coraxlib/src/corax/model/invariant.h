#ifndef CORAX_MODEL_INVARIANT_H_
#define CORAX_MODEL_INVARIANT_H_

#include "corax/corax_core.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in invariant.c */

  CORAX_EXPORT unsigned int
  corax_count_invariant_sites(const corax_partition_t *partition,
                              unsigned int *           state_inv_count);

  CORAX_EXPORT int corax_update_invariant_sites(corax_partition_t *partition);

  CORAX_EXPORT int
  corax_update_invariant_sites_proportion(corax_partition_t *partition,
                                          unsigned int       params_index,
                                          double             prop_invar);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CORAX_MODEL_INVARIANT_H_ */
