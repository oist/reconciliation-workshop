#ifndef CORAX_KERNEL_DERIVATIVES_H_
#define CORAX_KERNEL_DERIVATIVES_H_

#include "corax/core/partition.h"

#ifdef __cplusplus
extern "C"
{
#endif

  /* functions in derivatives.c */

  /**
   * Function which computes a "sumtable". This sumtable can be used to compute
   * the derivative of the likelihood with respect to a branch length.
   *
   * @param partition The partition for which the sumtable will be computed.
   *
   * @param parent_clv_index Parent CLV index of the edge in question.
   *
   * @param child_clv_index Child CLV index of the edge in question.
   *
   * @param params_indices A list of the indices for each rate category present
   * in the partition.
   *
   * @param[out] sumtable Buffer for the resulting sumtable. Should be allocated
   * with `rates * states_padded` elements.
   *
   * @ingroup corax_partition_t
   */
  CORAX_EXPORT int corax_update_sumtable(corax_partition_t *partition,
                                         unsigned int       parent_clv_index,
                                         unsigned int       child_clv_index,
                                         int                parent_scaler_index,
                                         int                child_scaler_index,
                                         const unsigned int *params_indices,
                                         double *            sumtable);

  /**
   * Computes the first and second derivative with respect to a specific branch
   * length.
   *
   * @param partition Partition that the derivative is computed for.
   *
   * @param parent_scaler_index Scaler index for the parent of the edge in
   * question.
   *
   * @param child_scaler_index Scaler index for the child of the edge in
   * question.
   *
   * @param branch_length Value at which to evaluate the derivative at.
   *
   * @param sumtable Sumbtable from `corax_udate_sumtable`.
   *
   * @param[out] d_f Buffer to store the first derivative. Only a single double.
   *
   * @param[out] dd_f Buffer to store the second derivative. Only a single
   * double.
   *
   * @ingroup corax_partition_t
   */
  CORAX_EXPORT int
  corax_compute_likelihood_derivatives(corax_partition_t * partition,
                                       int                 parent_scaler_index,
                                       int                 child_scaler_index,
                                       double              branch_length,
                                       const unsigned int *params_indices,
                                       const double *      sumtable,
                                       double *            d_f,
                                       double *            dd_f);

#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CORAX_KERNEL_DERIVATIVES_H_ */
