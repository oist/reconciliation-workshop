#ifndef CORAX_OPTIMIZE_BRANCHES_H_
#define CORAX_OPTIMIZE_BRANCHES_H_

#include "opt_generic.h"

/* Custom parameters structure provided by PLL for the
 * high level optimization functions (Newton-Raphson). */
typedef struct
{
  corax_partition_t * partition;
  corax_unode_t *     tree;
  const unsigned int *params_indices;
  double *            sumtable;
  double              branch_length_min;
  double              branch_length_max;
  double              tolerance;
  int                 max_newton_iters;
  int                 opt_method; /* see CORAX_OPT_BLO_* constants above */
} corax_newton_tree_params_t;

typedef struct
{
  corax_unode_t *     tree;
  unsigned int        partition_count;
  corax_partition_t **partitions;
  unsigned int **     params_indices;
  double **           precomp_buffers;
  double **           brlen_buffers;
  double *            brlen_orig;
  double *            brlen_guess;
  int *               converged;
  double *            brlen_scalers;
  double              branch_length_min;
  double              branch_length_max;
  double              tolerance;
  int                 max_newton_iters;
  int                 opt_method; /* see CORAX_OPT_BLO_* constants above */
  int                 brlen_linkage;
  void *              parallel_context;
  void (*parallel_reduce_cb)(void *, double *, size_t, int);
} corax_newton_tree_params_multi_t;

#ifdef __cplusplus
extern "C"
{
#endif

  CORAX_EXPORT void corax_opt_derivative_func(void *  parameters,
                                              double  proposal,
                                              double *df,
                                              double *ddf);

  /* high level optimization functions */

  CORAX_EXPORT double
  corax_opt_optimize_branch_lengths_iterative(corax_partition_t * partition,
                                              corax_unode_t *     tree,
                                              const unsigned int *params_indices,
                                              double branch_length_min,
                                              double branch_length_max,
                                              double tolerance,
                                              int    smoothings,
                                              int    keep_update);

  CORAX_EXPORT double
  corax_opt_optimize_branch_lengths_local(corax_partition_t * partition,
                                          corax_unode_t *     tree,
                                          const unsigned int *params_indices,
                                          double              branch_length_min,
                                          double              branch_length_max,
                                          double              tolerance,
                                          int                 smoothings,
                                          int                 radius,
                                          int                 keep_update);

  CORAX_EXPORT double corax_opt_optimize_branch_lengths_local_multi(
      corax_partition_t **partitions,
      size_t              partition_count,
      corax_unode_t *     tree,
      unsigned int **     params_indices,
      double **           sumtable_buffers,
      double **           brlen_buffers,
      double *            brlen_scalers,
      double              branch_length_min,
      double              branch_length_max,
      double              lh_epsilon,
      int                 max_iters,
      int                 radius,
      int                 keep_update,
      int                 opt_method,
      int                 brlen_linkage,
      void *              parallel_context,
      void (*parallel_reduce_cb)(void *, double *, size_t, int));

  // Optimizes branches around the local quartet defined by the root of the tree.
  // it also works for multi-partition data.
  // This functions assumes that p-matrix indices are up-to-date.
  CORAX_EXPORT double corax_opt_optimize_branch_lengths_local_multi_quartet(
      corax_partition_t **partitions,
      size_t              partition_count,
      corax_unode_t *     tree,
      unsigned int **     params_indices,
      double **           sumtable_buffers,
      double **           brlen_buffers,
      double *            brlen_scalers,
      double              branch_length_min,
      double              branch_length_max,
      double              lh_epsilon,
      int                 max_iters,
      int                 keep_update,
      int                 opt_method,
      int                 brlen_linkage,
      void *              parallel_context,
    void (*parallel_reduce_cb)(void *, double *, size_t, int));


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CORAX_OPTIMIZE_BRANCHES_H_ */
