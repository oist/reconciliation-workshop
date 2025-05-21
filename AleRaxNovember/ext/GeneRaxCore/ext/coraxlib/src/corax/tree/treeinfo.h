#ifndef CORAX_TREE_TREEINFO_H_
#define CORAX_TREE_TREEINFO_H_

#include "corax/tree/utree.h"
#include "corax/tree/utree_split.h"

#define CORAX_TREEINFO_PARTITION_ALL -1

typedef struct treeinfo_edge
{
  unsigned int left_index;
  unsigned int right_index;
  unsigned int pmatrix_index;
  double       brlen;
} corax_treeinfo_edge_t;

typedef struct treeinfo_topology
{
  unsigned int           edge_count;
  unsigned int           brlen_set_count;
  unsigned int           root_index;
  corax_treeinfo_edge_t *edges;
  double **              branch_lengths;
} corax_treeinfo_topology_t;

typedef struct treeinfo
{
  // dimensions
  unsigned int tip_count;
  unsigned int partition_count;

  /* 0 = linked/shared, 1 = linked with scaler, 2 = unlinked */
  int     brlen_linkage;
  double *linked_branch_lengths;

  corax_unode_t *root;
  corax_utree_t *tree;

  unsigned int    subnode_count;
  corax_unode_t **subnodes;

  // partitions & partition-specific stuff
  corax_partition_t **partitions;
  double *            alphas;
  int *gamma_mode; /* discrete GAMMA rates computation mode (mean, median) */
  unsigned int **param_indices;
  int **         subst_matrix_symmetries;
  double **      branch_lengths;
  double *       brlen_scalers;
  double *       partition_loglh;
  int *          params_to_optimize;

  // partition that have been initialized (useful for parallelization)
  unsigned int        init_partition_count;
  unsigned int *      init_partition_idx;
  corax_partition_t **init_partitions;

  /* tree topology constraint */
  unsigned int * constraint;         /* legacy: constraint check vector */
  corax_split_set_t * cons_splits;   /* constraint tree splits */
  corax_split_set_t * tree_splits;   /* current tree splits for constraint check */

  /* precomputation buffers for derivatives (aka "sumtable") */
  double **deriv_precomp;

  // invalidation flags
  char **clv_valid;
  char **pmatrix_valid;

  // buffers
  corax_unode_t **   travbuffer;
  unsigned int *     matrix_indices;
  corax_operation_t *operations;

  // partition on which all operations should be performed
  int active_partition;

  // general-purpose counter
  unsigned int counter;

  // parallelization stuff
  void *parallel_context;
  void (*parallel_reduce_cb)(void *, double *, size_t, int);
} corax_treeinfo_t;

typedef struct
{
  unsigned int    node_count;
  corax_unode_t **nodes;

  unsigned int  partition_count;
  unsigned int *partition_indices;

  corax_utree_t *tree;
  double **      probs;
} corax_ancestral_t;

#ifdef __cplusplus
extern "C"
{
#endif

  /** @defgroup corax_treeinfo_t corax_treeinfo_t
   */

  /**
   * Create a corax_treeinfo_t from an existing tree.
   *
   * @param root A pointer to the virtual root of the unrooted tree.
   *
   * @param tips: Number of tips in the tree. Almost always this will also be the
   * number of taxa.
   *
   * @param partitions: Number of partitions that will be used in the full
   * analysis.
   *
   * @param brlen_linkage Which branch length linking method to use. Options are:
   * - `CORAX_COMMON_BRLEN_UNLINKED`: The branch lengths for one partition have
   *   no relation to any other partition.
   * - `CORAX_BRLEN_SCALED`: The branch lengths are scaled per partition.
   * - `CORAX_BRLEN_LINKED`: The branch lengths are all equal for all
   *   partitions.
   *
   * @ingroup corax_treeinfo_t
   */
  CORAX_EXPORT corax_treeinfo_t *corax_treeinfo_create(corax_unode_t *root,
                                                       unsigned int   tips,
                                                       unsigned int   partitions,
                                                       int brlen_linkage);

  CORAX_EXPORT
  int corax_treeinfo_set_parallel_context(
      corax_treeinfo_t *treeinfo,
      void *            parallel_context,
      void (*parallel_reduce_cb)(void *, double *, size_t, int op));

  CORAX_EXPORT
  void corax_treeinfo_parallel_reduce(corax_treeinfo_t *treeinfo,
                                      double           *data,
                                      size_t            size,
                                      int               op);

  /**
   * Initialize a partition in a treeinfo.
   *
   * @param partition_index Index of the partition to initialize.
   *
   * @param partition The pointer to the partition itself. This partition needs to
   * be initialized before this.
   *
   * @param params_to_optimize Which paramters to optimize. Options are:
   * - CORAX_OPT_PARAM_ALL
   * - CORAX_OPT_PARAM_SUBST_RATES
   * - CORAX_OPT_PARAM_ALPHA
   * - CORAX_OPT_PARAM_PINV
   * - CORAX_OPT_PARAM_FREQUENCIES
   * - CORAX_OPT_PARAM_BRANCHES_SINGLE
   * - CORAX_OPT_PARAM_BRANCHES_ALL
   * - CORAX_OPT_PARAM_BRANCHES_ITERATIVE
   * - CORAX_OPT_PARAM_TOPOLOGY
   * - CORAX_OPT_PARAM_FREE_RATES
   * - CORAX_OPT_PARAM_RATE_WEIGHTS
   * - CORAX_OPT_PARAM_BRANCH_LEN_SCALAR
   * - CORAX_OPT_PARAM_USER: Uses user defined code.
   * Any of these can be combined via the `|` operation.
   *
   * @param gamma_mode Controls the gamma rate discretization methods. Options
   * are:
   * - CORAX_GAMMA_RATES_MEAN
   * - CORAX_GAMMA_RATES_MEDIAN
   *
   * @param alpha Initial alpha to use in the model.
   *
   * @param param_indices Specify the parameter indices to ... do stuff. Can be
   * set to null, at which point the defaults are used.
   *
   * @param subst_matrix_symmetries A list of symmetries in the model. Can be set
   * to `nullptr` as well. If set to null, indicates that there are no symmetries
   * in the model. If there are symmetries, the order is "left to right". For
   * example, if our symmetry array was {0,0,0,1,1,1}, then the two symmetry
   * groups would be the top row, and the reamaining triangle
   *
   * @ingroup corax_treeinfo_t
   */
  CORAX_EXPORT int
  corax_treeinfo_init_partition(corax_treeinfo_t *  treeinfo,
                                unsigned int        partition_index,
                                corax_partition_t * partition,
                                int                 params_to_optimize,
                                int                 gamma_mode,
                                double              alpha,
                                const unsigned int *param_indices,
                                const int *         subst_matrix_symmetries);

  CORAX_EXPORT int corax_treeinfo_set_active_partition(corax_treeinfo_t *treeinfo,
                                                       int partition_index);

  CORAX_EXPORT int corax_treeinfo_set_root(corax_treeinfo_t *treeinfo,
                                           corax_unode_t *   root);

  CORAX_EXPORT
  int corax_treeinfo_get_branch_length_all(const corax_treeinfo_t *treeinfo,
                                           const corax_unode_t *   edge,
                                           double *                lengths);

  CORAX_EXPORT int corax_treeinfo_set_branch_length(corax_treeinfo_t *treeinfo,
                                                    corax_unode_t *   edge,
                                                    double            length);

  CORAX_EXPORT
  int corax_treeinfo_set_branch_length_all(corax_treeinfo_t *treeinfo,
                                           corax_unode_t *   edge,
                                           const double *    lengths);

  CORAX_EXPORT
  int corax_treeinfo_set_branch_length_partition(corax_treeinfo_t *treeinfo,
                                                 corax_unode_t *   edge,
                                                 int    partition_index,
                                                 double length);

  CORAX_EXPORT
  corax_utree_t *
  corax_treeinfo_get_partition_tree(const corax_treeinfo_t *treeinfo,
                                    int                     partition_index);

  CORAX_EXPORT
  corax_treeinfo_topology_t *
  corax_treeinfo_get_topology(const corax_treeinfo_t *   treeinfo,
                              corax_treeinfo_topology_t *topol);

  CORAX_EXPORT
  int corax_treeinfo_set_topology(corax_treeinfo_t *               treeinfo,
                                  const corax_treeinfo_topology_t *topol);

  CORAX_EXPORT
  int corax_treeinfo_destroy_topology(corax_treeinfo_topology_t *topol);

  CORAX_EXPORT int corax_treeinfo_destroy_partition(corax_treeinfo_t *treeinfo,
                                                    unsigned int partition_index);

  CORAX_EXPORT void corax_treeinfo_destroy(corax_treeinfo_t *treeinfo);

  CORAX_EXPORT int corax_treeinfo_update_prob_matrices(corax_treeinfo_t *treeinfo,
                                                       int update_all);

  CORAX_EXPORT void corax_treeinfo_invalidate_all(corax_treeinfo_t *treeinfo);

  CORAX_EXPORT int corax_treeinfo_validate_clvs(corax_treeinfo_t *treeinfo,
                                                corax_unode_t **  travbuffer,
                                                unsigned int travbuffer_size);

  CORAX_EXPORT void corax_treeinfo_invalidate_pmatrix(corax_treeinfo_t *treeinfo,
                                                      const corax_unode_t *edge);

  CORAX_EXPORT void corax_treeinfo_invalidate_clv(corax_treeinfo_t *   treeinfo,
                                                  const corax_unode_t *edge);

  CORAX_EXPORT double corax_treeinfo_compute_loglh(corax_treeinfo_t *treeinfo,
                                                   int               incremental);

  CORAX_EXPORT double corax_treeinfo_compute_loglh_flex(
      corax_treeinfo_t *treeinfo, int incremental, int update_pmatrices);

  CORAX_EXPORT double corax_treeinfo_compute_loglh_persite(
      corax_treeinfo_t *treeinfo, int incremental, int update_matrices, double **persite_lnl);

  CORAX_EXPORT
  int corax_treeinfo_scale_branches_all(corax_treeinfo_t *treeinfo,
                                        double            scaler);

  CORAX_EXPORT
  int corax_treeinfo_scale_branches_partition(corax_treeinfo_t *treeinfo,
                                              unsigned int      partition_idx,
                                              double            scaler);

  CORAX_EXPORT
  int corax_treeinfo_normalize_brlen_scalers(corax_treeinfo_t *treeinfo);

  CORAX_EXPORT int corax_treeinfo_set_tree(corax_treeinfo_t *treeinfo,
                                           corax_utree_t *   tree);
  /* topological constraint management */

  CORAX_EXPORT int
  corax_treeinfo_set_constraint_clvmap(corax_treeinfo_t *treeinfo,
                                       const int *       clv_index_map);

  /**
   * Set a new constraint tree and initialize internal data structs
   *
   * @param  cons_tree       tree to be used as topological constraint (multifurcated/incomplete)
   * @param  fast_and_dirty  1 = use legacy RAxML algorithm (buggy with incomplete trees!)
   *                         0 = use new split-based check algorithm (default)
   *
   * @return CORAX_SUCCESS if constraint was set successfully
   *         CORAX_FAILURE on error
   */
  CORAX_EXPORT int
  corax_treeinfo_set_constraint_tree(corax_treeinfo_t *   treeinfo,
                                     const corax_utree_t *cons_tree,
                                     int fast_and_dirty);


  /**
   * Check if an NNI is compatible with the current topological constraint
   *
   * @param  edge      NNI edge
   * @param  nni_type  CORAX_UTREE_MOVE_NNI_LEFT or CORAX_UTREE_MOVE_NNI_RIGHT
   *
   * @return CORAX_SUCCESS if NNI is compatible with constraint (or no constraint set)
   *         CORAX_FAILURE otherwise
   */
  CORAX_EXPORT int corax_treeinfo_constraint_check_nni(corax_treeinfo_t * treeinfo,
                                                       corax_unode_t * edge,
                                                       int nni_type);


  /**
   * Check if an SPR is compatible with the current topological constraint
   *
   * @param  subtree       pruned subtree
   * @param  regraft_edge  re-instertion edge
   *
   * @return CORAX_SUCCESS if SPR is compatible with constraint (or no constraint set)
   *         CORAX_FAILURE otherwise
   */
  CORAX_EXPORT int corax_treeinfo_constraint_check_spr(corax_treeinfo_t * treeinfo,
                                                       corax_unode_t * subtree,
                                                       corax_unode_t * regraft_edge);

  /**
   * Check if current constraint is relevant for a given subtree. For instance, a subtree comprising
   * only "free" taxa (e.g. those absent from the constraint tree) is not affected,
   * and hence constraint check is not required before regrafting this subtree.
   *
   * @param  subtree       pruned subtree
   *
   * @return CORAX_SUCCESS subtree is affected by the constraint
   *         CORAX_FAILURE otherwise
   */
  CORAX_EXPORT int corax_treeinfo_constraint_subtree_affected(corax_treeinfo_t * treeinfo,
                                                              corax_unode_t * subtree);


  /**
   * Extract all directed splits from current topology in treeinfo->tree
   * and store them in treeinfo->tree_splits.
   *
   * @return CORAX_SUCCESS if extraction was successful
   *         CORAX_FAILURE on error
   */
  CORAX_EXPORT int corax_treeinfo_constraint_update_splits(corax_treeinfo_t * treeinfo);


  /**
   * Check if current topology in treeinfo->tree is compatible with topological constraint.
   *
   * @return CORAX_SUCCESS if topology is compatible
   *         CORAX_FAILURE otherwise
   */
  CORAX_EXPORT int corax_treeinfo_constraint_check_current(corax_treeinfo_t * treeinfo);


  CORAX_EXPORT corax_ancestral_t *
               corax_treeinfo_compute_ancestral(corax_treeinfo_t *treeinfo);

  CORAX_EXPORT void
  corax_treeinfo_destroy_ancestral(corax_ancestral_t *ancestral);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CORAX_TREE_TREEINFO_H_ */
