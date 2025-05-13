#ifndef CORAX_KERNEL_PMATRIX_H_
#define CORAX_KERNEL_PMATRIX_H_

#include "corax/core/partition.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /**
   * Update the probability matrices of partition.
   *
   * @param params_index Index of the parameters to use, as in rate categories.
   *
   * @param matrix_indices An array of indices into the `prob_matrices` in
   * `corax_partition_t`. These are the locations in which the matrices will be
   * stored. The best way to get these is via `corax_utree_create_operations`.
   *
   * @param branch_lengths A list of branch lengths which will be used for
   * computing the probability matrices. The best way to get these is via
   * `corax_utree_create_operations`.
   *
   * @param count The number of matrices to update.
   *
   * @ingroup corax_partition_t
   */
  CORAX_EXPORT int
  corax_update_prob_matrices(corax_partition_t * partition,
                             const unsigned int *params_index,
                             const unsigned int *matrix_indices,
                             const double *      branch_lengths,
                             unsigned int        count);

#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CORAX_KERNEL_PMATRIX_H_ */
