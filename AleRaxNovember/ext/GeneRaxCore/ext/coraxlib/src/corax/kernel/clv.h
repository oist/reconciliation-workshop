#ifndef CORAX_KERNEL_CLV_H_
#define CORAX_KERNEL_CLV_H_

#include "corax/core/partition.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in clvs.c */

  /**
   * Computes the CLVS for all the trees in the nodes specified in the
   * `operations` array.
   *
   * @param[in,out] partition The partition for which the operations will be
   * computed.
   *
   * @param operations The list of operations. Typically this will be generated
   * using corax_utree_create_operations
   *
   * @param count Number of elements in the operations buffer.
   *
   * @ingroup corax_partition_t
   * @ingroup corax_operation_t
   */
  CORAX_EXPORT void corax_update_clvs(corax_partition_t *      partition,
                                      const corax_operation_t *operations,
                                      unsigned int             count);

  CORAX_EXPORT void corax_update_clvs_rep(corax_partition_t *      partition,
                                          const corax_operation_t *operations,
                                          unsigned int             count,
                                          unsigned int update_repeats);

#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CORAX_KERNEL_CLV_H_ */
