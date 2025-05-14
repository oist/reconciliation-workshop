/*
    Copyright (C) 2016 Diego Darriba

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

    Contact: Diego Darriba <Diego.Darriba@h-its.org>,
    Exelixis Lab, Heidelberg Instutute for Theoretical Studies
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#ifndef CORAX_OPTIMIZE_TREEINFO_H_
#define CORAX_OPTIMIZE_TREEINFO_H_

#include "corax/tree/treeinfo.h"
#include "opt_generic.h"
#include <stdio.h>

typedef struct cutoff_info
{
  double lh_start;
  double lh_cutoff;
  double lh_dec_sum;
  int    lh_dec_count;
} cutoff_info_t;

typedef int (*treeinfo_param_set_cb)(corax_treeinfo_t *treeinfo,
                                     unsigned int      part_num,
                                     const double     *param_vals,
                                     unsigned int      param_count);

typedef int (*treeinfo_param_get_cb)(const corax_treeinfo_t *treeinfo,
                                     unsigned int            part_num,
                                     double                 *param_vals,
                                     unsigned int            param_count);

#ifdef __cplusplus
extern "C"
{
#endif

  /* functions to optimize multiple partitions in parallel, using treeinfo
   * struct
   */

  CORAX_EXPORT double corax_algo_opt_onedim_treeinfo(corax_treeinfo_t *treeinfo,
                                                     int    param_to_optimize,
                                                     double min_value,
                                                     double max_value,
                                                     double tolerance);

  CORAX_EXPORT
  double
  corax_algo_opt_onedim_treeinfo_custom(corax_treeinfo_t     *treeinfo,
                                        int                   param_to_optimize,
                                        treeinfo_param_get_cb params_getter,
                                        treeinfo_param_set_cb params_setter,
                                        double                min_value,
                                        double                max_value,
                                        double                tolerance);

  /** @defgroup treeinfo_opt_algs Treeinfo Optimization Algorithms
   * These functions more or less have the same arguments and optimize their
   * particular parameter. They all use similar algorithms.
   *
   * @param bfgs_factor: A value that should be greater than 0. An opaque
   * quantity, but typical values are 1e12 for inaccurate but fast searches, or
   * 1e1 for a more complete and slow search.
   *
   * @param tolerance This controls the pgtol threshold. If the largest entry of
   * the projected gradient is larger than this value, optimization will stop.
   *
   * @ingroup corax_treeinfo_t
   *
   * @{
   */
  CORAX_EXPORT
  double corax_algo_opt_subst_rates_treeinfo(corax_treeinfo_t *treeinfo,
                                             unsigned int      params_index,
                                             double            min_rate,
                                             double            max_rate,
                                             double            bfgs_factor,
                                             double            tolerance);

  CORAX_EXPORT
  double corax_algo_opt_frequencies_treeinfo(corax_treeinfo_t *treeinfo,
                                             unsigned int      params_index,
                                             double            min_freq,
                                             double            max_freq,
                                             double            bfgs_factor,
                                             double            tolerance);

  CORAX_EXPORT
  double corax_algo_opt_rates_weights_treeinfo(corax_treeinfo_t *treeinfo,
                                               double            min_rate,
                                               double            max_rate,
                                               double            min_brlen,
                                               double            max_brlen,
                                               double            bfgs_factor,
                                               double            tolerance);

  CORAX_EXPORT
  double corax_algo_opt_alpha_pinv_treeinfo(corax_treeinfo_t *treeinfo,
                                            unsigned int      params_index,
                                            double            min_alpha,
                                            double            max_alpha,
                                            double            min_pinv,
                                            double            max_pinv,
                                            double            bfgs_factor,
                                            double            tolerance);

  CORAX_EXPORT
  double corax_algo_opt_brlen_scalers_treeinfo(corax_treeinfo_t *treeinfo,
                                               double            min_scaler,
                                               double            max_scaler,
                                               double            min_brlen,
                                               double            max_brlen,
                                               double            lh_epsilon);

  CORAX_EXPORT
  double corax_algo_opt_brlen_treeinfo(corax_treeinfo_t *treeinfo,
                                       double            min_brlen,
                                       double            max_brlen,
                                       double            lh_epsilon,
                                       int               max_iters,
                                       int               opt_method,
                                       int               radius);
  /** @} */

  /* search */

#ifdef __cplusplus
/* in case the compiler is a C++ compiler */
#define DEFAULT_VALUE(value) = value
#else
/* otherwise, C compiler, do nothing */
#define DEFAULT_VALUE(value)
#endif

  /**
   * Perform an SPR round.
   *
   * @param radius_min
   * @param radius_max Minimum and maximum thresholds for the SPR reinsertion
   * distance.
   *
   * @param brlen_opt_method The optimization method to use when optimizing
   * branch lengths. Options are:
   * - CORAX_OPT_BLO_NEWTON_FAST: Standard.
   * - CORAX_OPT_BLO_NEWTON_SAFE: Adds a per branch likelihood check.
   * - CORAX_OPT_BLO_NEWTON_FALLBACK: Starts fast, but fallback to safe.
   * - CORAX_OPT_BLO_NEWTON_GLOBAL: Newton, but with additional searches to find
   *   more optima
   * - CORAX_OPT_BLO_NEWTON_OLDFAST
   * - CORAX_OPT_BLO_NEWTON_OLDSAFE
   *
   * @param smoothings: Maximum number of iterations for branch length
   * optimization. Negative = no limit (iterate until LH improvement < epsilon)
   *
   * @param epsilon Likelihood threshold to terminate the optimization. Also
   * known as the tolerance.
   *
   * @param[out] cutoff_info A struct that contains subtree descent cutoff
   * information. It is in/out parameter since cutoff info has to be preserved
   * between subsequent SPR rounds.
   *
   * @param subtree_cutoff relative likelihood cutoff for descending into
   * subtrees. A larger value means higher cutoff, i.e. deeper descent into
   * subtrees. For more information, see
   * https://cme.h-its.org/exelixis/pubs/VLSI2007.pdf
   * 
   * @param lh_epsilon_brlen_triplet Epsilon value for branch length optimization
   * of the triplet of nodes around the regrafting point (e.g. 0.1)
   * 
   * @param fast_clv_updates Enable fast clv updates (default value: 1)
   *
   * @ingroup corax_treeinfo_t
   */
  CORAX_EXPORT double corax_algo_spr_round(corax_treeinfo_t *treeinfo,
                                           unsigned int      radius_min,
                                           unsigned int      radius_max,
                                           unsigned int      ntopol_keep,
                                           corax_bool_t      thorough,
                                           int               brlen_opt_method,
                                           double            bl_min,
                                           double            bl_max,
                                           int               smoothings,
                                           double            epsilon,
                                           cutoff_info_t    *cutoff_info,
                                           double            subtree_cutoff,
                                           double            lh_epsilon_brlen_triplet,
                                           corax_bool_t fast_clv_updates DEFAULT_VALUE(1),
                                           unsigned long int *total_moves_counter DEFAULT_VALUE(NULL),
                                           unsigned long int *improving_moves_counter DEFAULT_VALUE(NULL));

// Defining errors in NNI
#define CORAX_NNI_ROUND_LEAF_ERROR 6001
#define CORAX_NNI_ROUND_INTEGRITY_ERROR 6002
#define CORAX_NNI_ROUND_TRIPLET_ERROR 6003
#define CORAX_NNI_ROUND_UNDO_MOVE_ERROR 6004
#define CORAX_NNI_DIFF_NEGATIVE_ERROR 6005
#define CORAX_NNI_ROOT_NOT_FOUND 6005

  /**
   * SH-like aLRT statistics calculation (support values for internal branches).
   * SH-like aLRT values are defined for NNI optimal tree topologies. Hence, it
   * is recommended for the user to call corax_algo_nni_round() function first.
   * In case the tree topology is not NNI-optimal, a warning message will be
   * printed.
   *
   * @param  treeinfo               the CORAX treeinfo structure
   * @param  tolerance              tolerance value, to check for NNI-optimality
   * and to avoid numerical issues. It can be equal to tolerance value used for
   * the NNI otpimization in  corax_algo_nni_round() function (e.g. 0.1)
   * @param[out] sh_support_values   Array where the SH-aLRT statistics will be
   * stored. The size of the array must be 2*n-3, where n is the number of tip
   * nodes. The SH-aLRT statistic for a branch with pmatrix_index = i is stored
   * in sh_support_values[i]. The SH-like aLRT metric for non NNI-optimal internal
   * branches will be SH-aLRT=-inf. Tip branches will also have SH-aLRT=-inf.
   * For ambiguous branches, that is, internal branches in which there is a
   * second NNI-optimal topology, the SH-like aLRT metric will be equal to 0.
   * If sh_support_values==NULL, no values will be written (useful for multi-threading).
   * @param  num_bsrep              Number of bootstrap replicates (e.g. 1000)
   * @param  bsrep_site_weights     Resampled site weights. This two-dimensional array
   * is indexed by replicate+partition, and then by site. I.e., the weight of site s of
   * partition p in the replicate i must be stored in
   * bsrep_site_weights[i * treeinfo->partition_count + p][s]
   * In multi-threading use case, only the portion of sites processed by the
   * current thread must be provided (similar to treeinfo->partitions[]).
   * @param  sh_epsilon             Confidence of SH-like criterion (e.g. 0.1)
   * @param  brlen_opt_method       Branch length optimization method (e.g.
   * CORAX_OPT_BLO_NEWTON_FAST)
   * @param  bl_min                 Minimum branch length (e.g.
   * CORAX_OPT_MIN_BRANCH_LEN)
   * @param  bl_max                 Maximum branch length (e.g.
   * CORAX_OPT_MAX_BRANCH_LEN)
   * @param  smoothings             number of smoothings in local branch length
   * optimization that takes place (e.g. CORAX_OPT_DEFAULT_SMOOTHINGS)
   * @param  lh_epsilon             epsilon value in local branch length
   * optimization that takes place (e.g. CORAX_OPT_DEFAULT_EPSILON)
   * @return                        CORAX_SUCCESS, if the SH-like aLRT values
   * are calculated successfully
   */
  CORAX_EXPORT int
  corax_algo_sh_support(corax_treeinfo_t        *treeinfo,
                        double                  tolerance,
                        double                 *sh_support_values,
                        unsigned int            num_bsrep,
                        const unsigned int    **bsrep_site_weights,
                        double                  sh_epsilon,
                        int                     brlen_opt_method,
                        double                  bl_min,
                        double                  bl_max,
                        int                     smoothings,
                        double                  lh_epsilon);


  /**
   * DEPRECATED! Use corax_algo_sh_support() instead.
   *
   * SH-like aLRT statistics calculation (support values for internal branches).
   * SH-like aLRT values are defined for NNI optimal tree topologies. Hence, it
   * is recommended for the user to call corax_algo_nni_round() function first.
   * In case the tree topology is not NNI-optimal, a warning message will be
   * printed.
   *
   * @param  treeinfo               the CORAX treeinfo structure
   * @param  tolerance              tolerance value, to check for NNI-optimality
   * and to avoid numerical issues. It can be equal to tolerance value used for
   * the NNI otpimization in  corax_algo_nni_round() function (e.g. 0.1)
   * @param[out] shSupportValues    Array where the SH-aLRT statistics will be
   * stored. The size of the array must be 2*n-3, where n is the number of tip
   * nodes. The SH-aLRT statistic for a branch with pmatrix_index = i is stored
   * in shSupportValues[i]. The SH-like aLRT metric for non NNI-optimal internal
   * branches will be SH-aLRT=-inf. Tip branches will also have SH-aLRT=-inf.
   * For ambiguous branches, that is, internal branches in which there is a
   * second NNI-optimal topology, the SH-like aLRT metric will be equal to 0.
   * @param  nBootstrap             Number of bootstrap replicates (e.g. 1000)
   * @param  shEpsilon              Confidence of SH-like criterion (e.g. 0.1)
   * @param  brlen_opt_method       Branch length optimization method (e.g.
   * CORAX_OPT_BLO_NEWTON_FAST)
   * @param  bl_min                 Minimum branch length (e.g.
   * CORAX_OPT_MIN_BRANCH_LEN)
   * @param  bl_max                 Maximum branch length (e.g.
   * CORAX_OPT_MAX_BRANCH_LEN)
   * @param  smoothings             number of smoothings in local branch length
   * optimization that takes place (e.g. CORAX_OPT_DEFAULT_SMOOTHINGS)
   * @param  lh_epsilon             epsilon value in local branch length
   * optimization that takes place (e.g. CORAX_OPT_DEFAULT_EPSILON)
   * @param  print_in_console       TRUE, if anything to be printed in the
   * console, otherwise FALSE. (Default TRUE)
   * @return                        CORAX_SUCCESS, if the SH-like aLRT values
   * are calculated successfully
   */
  CORAX_EXPORT int
  corax_shSupport_values(corax_treeinfo_t     *treeinfo,
                         double                tolerance,
                         double               *shSupportValues,
                         int                   nBootstrap,
                         double                shEpsilon,
                         int                   brlen_opt_method,
                         double                bl_min,
                         double                bl_max,
                         int                   smoothings,
                         double                lh_epsilon,
                         bool print_in_console DEFAULT_VALUE(true));

  /**
   * NNI round - Searches for the optimal tree topology based on NNI moves.
   * After calling this function the tree topology is probably changed.
   *
   * Check `corax_algo_nni_round` documentation.
   *
   * @param  treeinfo               the CORAX treeinfo structure
   * @param  tolerance              tolerance for NNI round: if (final_logl -
   * init_logl <= tolerance) -> exit
   * @param  brlen_opt_method       Branch length optimization method (e.g.
   * CORAX_OPT_BLO_NEWTON_FAST)
   * @param  bl_min                 Minimum branch length (e.g.
   * CORAX_OPT_MIN_BRANCH_LEN)
   * @param  bl_max                 Maximum branch length (e.g.
   * CORAX_OPT_MAX_BRANCH_LEN)
   * @param  smoothings             number of smoothings in local branch length
   * optimization that takes place (e.g. CORAX_OPT_DEFAULT_SMOOTHINGS)
   * @param  lh_epsilon             epsilon value in local branch length
   * optimization that takes place (e.g. CORAX_OPT_DEFAULT_EPSILON)
   * @param  print_in_console       TRUE, if the intermediate stages of NNI
   * round are to be printed in the console, otherwise FALSE. (Default TRUE)
   * @return                        the likelihood score after NNI optimization
   * + new topology
   */
  CORAX_EXPORT double
  corax_algo_nni_round(corax_treeinfo_t     *treeinfo,
                       double                tolerance,
                       int                   brlen_opt_method,
                       double                bl_min,
                       double                bl_max,
                       int                   smoothings,
                       double                lh_epsilon,
                       bool print_in_console DEFAULT_VALUE(true));

  /**
   * Finds the best out of the 3 NNI topologies, in the quartet defined around
   * the root of the tree. It is assumed that the rood in an internal branch. It
   * is also assumed that probability matrices and CLVs around the rood are up
   * to date.  Even if the current tree topology is the best out of 3, the
   * returned likelihood might be increased, since branch lenghts are optimized
   * before the likelihood calculations of the 3 topologies
   *
   * @param  treeinfo          the CORAX treeinfo structure
   * @param  brlen_opt_method  Branch length optimization method (e.g.
   * CORAX_OPT_BLO_NEWTON_FAST)
   * @param  bl_min            Minimum branch length (e.g.
   * CORAX_OPT_MIN_BRANCH_LEN)
   * @param  bl_max            Maximum branch length (e.g.
   * CORAX_OPT_MAX_BRANCH_LEN)
   * @param  smoothings        number of smoothings in local branch length
   * optimization that takes place (e.g. CORAX_OPT_DEFAULT_SMOOTHINGS)
   * @param  lh_epsilon        epsilon value in local branch length optimization
   * that takes place (e.g. CORAX_OPT_DEFAULT_EPSILON)
   * @return                   the likelihood score after NNI optimization + new
   * topology
   */
  CORAX_EXPORT double corax_algo_nni_local(corax_treeinfo_t *treeinfo,
                                           int               brlen_opt_method,
                                           double            bl_min,
                                           double            bl_max,
                                           int               smoothings,
                                           double            lh_epsilon);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* CORAX_OPTIMIZE_TREEINFO_H_ */
