/*
 Copyright (C) 2015-21 Diego Darriba, Alexey Kozlov

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

 Contact: Alexey Kozlov <Alexey.Kozlov@h-its.org>,
 Exelixis Lab, Heidelberg Instutute for Theoretical Studies
 Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
 */

#include "opt_branches.h"
#include "corax/corax_kernel.h"


#define BETTER_LL_TRESHOLD 1e-13

/*
 * Note: Compile with flag _ULTRACHECK for checking pre/postconditions
 *       way more thoroughly. This may slow down the execution.
 */
static inline int d_equals(double a, double b) { return (fabs(a - b) < 1e-10); }

static inline int check_loglh_improvement(int opt_method)
{
  return (opt_method == CORAX_OPT_BLO_NEWTON_SAFE) ? 1 : 0;
}

static void utree_derivative_func(void *  parameters,
                                  double  proposal,
                                  double *df,
                                  double *ddf)
{
  corax_newton_tree_params_t *params = (corax_newton_tree_params_t *)parameters;
  corax_compute_likelihood_derivatives(params->partition,
                                       params->tree->scaler_index,
                                       params->tree->back->scaler_index,
                                       proposal,
                                       params->params_indices,
                                       params->sumtable,
                                       df,
                                       ddf);
}

/******************************************************************************/
/* GENERIC */
/******************************************************************************/

static void update_clvs_and_scalers(corax_partition_t **partitions,
                                    size_t              partition_count,
                                    corax_unode_t *     parent,
                                    corax_unode_t *     right_child,
                                    corax_unode_t *     left_child)
{
  corax_operation_t op;
  size_t            p;

  /* set CLV */
  op.parent_clv_index    = parent->clv_index;
  op.parent_scaler_index = parent->scaler_index;
  op.child1_clv_index    = right_child->back->clv_index;
  op.child1_matrix_index = right_child->back->pmatrix_index;
  op.child1_scaler_index = right_child->back->scaler_index;
  op.child2_clv_index    = left_child->back->clv_index;
  op.child2_matrix_index = left_child->back->pmatrix_index;
  op.child2_scaler_index = left_child->back->scaler_index;

  for (p = 0; p < partition_count; ++p)
  {
    /* skip remote partitions */
    if (!partitions[p]) continue;

    corax_update_clvs(partitions[p], &op, 1);
  }
}

/* if keep_update, P-matrices are updated after each branch length opt */
static int recomp_iterative(corax_newton_tree_params_t *params,
                            int                         radius,
                            double *                    loglikelihood_score,
                            int                         keep_update)
{
  corax_unode_t *tr_p, *tr_q, *tr_z;
  double         xmin, /* min branch length */
      xguess,          /* initial guess */
      xmax,            /* max branch length */
      xtol,            /* tolerance */
      xres,            /* optimal found branch length */
      xorig;           /* original branch length before optimization */

  tr_p  = params->tree;
  tr_q  = params->tree->next;
  tr_z  = tr_q ? tr_q->next : NULL;
  xorig = tr_p->length;

  /* check branch length integrity */
  assert(d_equals(tr_p->length, tr_p->back->length));

  /* prepare sumtable for current branch */
  corax_update_sumtable(params->partition,
                        tr_p->clv_index,
                        tr_p->back->clv_index,
                        tr_p->scaler_index,
                        tr_p->back->scaler_index,
                        params->params_indices,
                        params->sumtable);

  /* set N-R parameters */
  xmin   = params->branch_length_min;
  xmax   = params->branch_length_max;
  xtol   = params->tolerance;
  xguess = tr_p->length;
  if (xguess < xmin || xguess > xmax) xguess = CORAX_OPT_DEFAULT_BRANCH_LEN;

  xres = corax_opt_minimize_newton(xmin,
                                   xguess,
                                   xmax,
                                   xtol,
                                   params->max_newton_iters,
                                   params,
                                   utree_derivative_func);

  if (corax_errno) return CORAX_FAILURE;

  /* update branch length in the tree structure */
  tr_p->length = tr_p->back->length = xres;

  if (keep_update && fabs(tr_p->length - xorig) > 1e-10)
  {
    /* update pmatrix for the new branch length */
    corax_update_prob_matrices(params->partition,
                               params->params_indices,
                               &(tr_p->pmatrix_index),
                               &xres,
                               1);

    if (check_loglh_improvement(params->opt_method))
    {
      /* check and compare likelihood */
      double eval_loglikelihood =
          corax_compute_edge_loglikelihood(params->partition,
                                           tr_p->clv_index,
                                           tr_p->scaler_index,
                                           tr_p->back->clv_index,
                                           tr_p->back->scaler_index,
                                           tr_p->pmatrix_index,
                                           params->params_indices,
                                           NULL);

      /* check if the optimal found value improves the likelihood score */
      if (eval_loglikelihood >= *loglikelihood_score)
      {
        /* fix new score */
        *loglikelihood_score = eval_loglikelihood;

        /* update branch length in the tree structure */
        tr_p->length       = xres;
        tr_p->back->length = tr_p->length;
      }
      else
      {
        /* reset branch length to original value */
        tr_p->length = tr_p->back->length = xorig;

        corax_update_prob_matrices(params->partition,
                                   params->params_indices,
                                   &(tr_p->pmatrix_index),
                                   &tr_p->length,
                                   1);
      }
    }
  }

  DBG(" Optimized branch %3d - %3d (%.6f)\n",
      tr_p->clv_index,
      tr_p->back->clv_index,
      tr_p->length);

  /* update children */
  if (radius && tr_q && tr_z)
  {
    /* update children 'Q'
     * CLV at P is recomputed with children P->back and Z->back
     * Scaler is updated by subtracting Q->back and adding P->back
     */
    update_clvs_and_scalers(&params->partition, 1, tr_q, tr_p, tr_z);

    /* eval */
    corax_newton_tree_params_t params_cpy;
    memcpy(&params_cpy, params, sizeof(corax_newton_tree_params_t));
    params_cpy.tree = tr_q->back;
    if (!recomp_iterative(
            &params_cpy, radius - 1, loglikelihood_score, keep_update))
      return CORAX_FAILURE;

    /* update children 'Z'
     * CLV at P is recomputed with children P->back and Q->back
     * Scaler is updated by subtracting Z->back and adding Q->back
     */
    update_clvs_and_scalers(&params->partition, 1, tr_z, tr_q, tr_p);

    /* eval */
    params_cpy.tree = tr_z->back;
    if (!recomp_iterative(
            &params_cpy, radius - 1, loglikelihood_score, keep_update))
      return CORAX_FAILURE;

    /* reset to initial state
     * CLV at P is recomputed with children Q->back and Z->back
     * Scaler is updated by subtracting P->back and adding Z->back
     */
    update_clvs_and_scalers(&params->partition, 1, tr_p, tr_z, tr_q);
  }

  return CORAX_SUCCESS;

} /* recomp_iterative */

static void utree_derivative_func_multi(void *  parameters,
                                        double *proposal,
                                        double *df,
                                        double *ddf)
{
  corax_newton_tree_params_multi_t *params =
      (corax_newton_tree_params_multi_t *)parameters;
  size_t p;
  int    unlinked = (params->brlen_linkage == CORAX_BRLEN_UNLINKED) ? 1 : 0;

  if (unlinked)
  {
    for (p = 0; p < params->partition_count; ++p) df[p] = ddf[p] = 0;
  }
  else
    *df = *ddf = 0;

  /* simply iterate over partitions and add up the derivatives */
  for (p = 0; p < params->partition_count; ++p)
  {
    /* skip remote partitions */
    if (!params->partitions[p]) continue;

    double p_df, p_ddf;
    double s       = params->brlen_scalers ? params->brlen_scalers[p] : 1.;
    double p_brlen = s * (unlinked ? proposal[p] : proposal[0]);
    corax_compute_likelihood_derivatives(params->partitions[p],
                                         params->tree->scaler_index,
                                         params->tree->back->scaler_index,
                                         p_brlen,
                                         params->params_indices[p],
                                         params->precomp_buffers[p],
                                         &p_df,
                                         &p_ddf);

    /* chain rule! */
    if (unlinked)
    {
      df[p]  = s * p_df;
      ddf[p] = s * s * p_ddf;
    }
    else
    {
      df[0] += s * p_df;
      ddf[0] += s * s * p_ddf;
    }
  }

  if (params->parallel_reduce_cb)
  {
    if (unlinked)
    {
      params->parallel_reduce_cb(params->parallel_context,
                                 df,
                                 params->partition_count,
                                 CORAX_REDUCE_SUM);
      params->parallel_reduce_cb(params->parallel_context,
                                 ddf,
                                 params->partition_count,
                                 CORAX_REDUCE_SUM);
    }
    else
    {
      double d[2] = {*df, *ddf};
      params->parallel_reduce_cb(
          params->parallel_context, d, 2, CORAX_REDUCE_SUM);
      *df  = d[0];
      *ddf = d[1];
    }
  }
}

static void update_prob_matrices(corax_partition_t **partitions,
                                 size_t              partition_count,
                                 unsigned int **     params_indices,
                                 double **           brlen_buffers,
                                 double *            brlen_scalers,
                                 corax_unode_t *     node)
{
  unsigned int p;
  unsigned int m = node->pmatrix_index;

  for (p = 0; p < partition_count; ++p)
  {
    /* skip remote partitions */
    if (!partitions[p]) continue;

    double p_brlen = brlen_buffers ? brlen_buffers[p][m] : node->length;

    if (brlen_scalers) p_brlen *= brlen_scalers[p];

    corax_update_prob_matrices(
        partitions[p], params_indices[p], &m, &p_brlen, 1);
  }
}

static int allocate_buffers(corax_newton_tree_params_multi_t *params)
{
  if (!params->precomp_buffers)
  {
    params->precomp_buffers =
        (double **)calloc(params->partition_count, sizeof(double *));
    if (!params->precomp_buffers) return CORAX_FAILURE;

    for (unsigned int p = 0; p < params->partition_count; ++p)
    {
      const corax_partition_t *partition = params->partitions[p];

      /* skip remote partitions */
      if (!partition) continue;

      unsigned int sites_alloc = partition->sites;
      if (partition->attributes & CORAX_ATTRIB_AB_FLAG)
        sites_alloc += partition->states;

      params->precomp_buffers[p] = (double *)corax_aligned_alloc(
          sites_alloc * partition->rate_cats * partition->states_padded
              * sizeof(double),
          partition->alignment);

      if (!params->precomp_buffers[p]) return CORAX_FAILURE;
    }
  }

  if (!params->brlen_buffers
      && params->opt_method == CORAX_OPT_BLO_NEWTON_FALLBACK)
  {
    params->brlen_buffers =
        (double **)calloc(params->partition_count, sizeof(double *));
    if (!params->brlen_buffers) return CORAX_FAILURE;

    // not very elegant...
    unsigned int branch_count = 0;
    for (size_t p = 0; p < params->partition_count; ++p)
    {
      if (params->partitions[p])
      {
        branch_count = 2 * params->partitions[p]->tips - 3;
        break;
      }
    }
    assert(branch_count);

    assert(params->brlen_linkage != CORAX_BRLEN_UNLINKED);

    params->brlen_buffers[0] = (double *)calloc(branch_count, sizeof(double));

    if (!params->brlen_buffers[0]) return CORAX_FAILURE;
  }

  if (params->brlen_linkage == CORAX_BRLEN_UNLINKED)
  {
    params->converged = (int *)calloc(params->partition_count, sizeof(int));
    params->brlen_orig =
        (double *)calloc(params->partition_count, sizeof(double));
    params->brlen_guess =
        (double *)calloc(params->partition_count, sizeof(double));
    if (!params->converged || !params->brlen_orig || !params->brlen_guess)
      return CORAX_FAILURE;
  }

  return CORAX_SUCCESS;
}

/**
 * Compute the likelihood score at a given edge on a multiple partition
 *
 * @param  partitions          list of partitions
 * @param  partition_count     number of partitions in `partitions`
 * @param  parent_clv_index    parent clv index
 * @param  parent_scaler_index parent scaler index
 * @param  child_clv_index     child clv index
 * @param  child_scaler_index  child scaler index
 * @param  matrix_index        matrix index of the edge
 * @param  params_indices      the indices of the parameter sets
 * @param  persite_lnl         per-site likelihoods (if 0, they are omitted)
 * @param  parallel_context    context for parallel computation
 * @param  parallel_reduce_cb  callback function for parallel reduction
 *
 * @return                     the likelihood score at the given edge
 */
static double compute_edge_loglikelihood_multi(
    corax_partition_t ** partitions,
    size_t               partition_count,
    unsigned int         parent_clv_index,
    int                  parent_scaler_index,
    unsigned int         child_clv_index,
    int                  child_scaler_index,
    unsigned int         matrix_index,
    unsigned int **const params_indices,
    double *             persite_lnl,
    void *               parallel_context,
    void (*parallel_reduce_cb)(void *, double *, size_t, int))
{
  CORAX_UNUSED(persite_lnl);
  double total_loglh = 0.;

  size_t p;
  for (p = 0; p < partition_count; ++p)
  {
    /* skip remote partitions */
    if (!partitions[p]) continue;

    total_loglh += corax_compute_edge_loglikelihood(partitions[p],
                                                    parent_clv_index,
                                                    parent_scaler_index,
                                                    child_clv_index,
                                                    child_scaler_index,
                                                    matrix_index,
                                                    params_indices[p],
                                                    NULL);
  }

  if (parallel_reduce_cb)
    parallel_reduce_cb(parallel_context, &total_loglh, 1, CORAX_REDUCE_SUM);

  return total_loglh;
}

/* if keep_update, P-matrices are updated after each branch length opt */
static int recomp_iterative_multi(corax_newton_tree_params_multi_t *params,
                                  int                               radius,
                                  double *loglikelihood_score,
                                  int     keep_update)
{
  corax_unode_t *tr_p, *tr_q, *tr_z;
  unsigned int   p;
  int            retval;
  unsigned int   xnum;
  double         xmin, /* min branch length */
      xorig_linked,    /* original branch length before optimization (linked) */
      xguess_linked,   /* initial guess (linked) */
      xmax,            /* max branch length */
      xtol;            /* tolerance */

  double *xorig, /* original branch length before optimization */
      *xguess;   /* initial guess / current branch length value */

  unsigned int pmatrix_index;

  int unlinked     = params->brlen_linkage == CORAX_BRLEN_UNLINKED ? 1 : 0;
  int apply_change = 0;

  tr_p          = params->tree;
  tr_q          = params->tree->next;
  tr_z          = tr_q ? tr_q->next : NULL;
  pmatrix_index = tr_p->pmatrix_index;

  if (unlinked)
  {
    assert(params->brlen_buffers && params->brlen_orig && params->brlen_guess);
    xnum   = params->partition_count;
    xorig  = params->brlen_orig;
    xguess = params->brlen_guess;

    /* set initial values */
    for (p = 0; p < xnum; ++p)
    {
      xguess[p] =
          params->partitions[p] ? params->brlen_buffers[p][pmatrix_index] : 0.0;
    }

    /* make sure every thread has all per-partition branch lengths */
    if (xnum > 1 && params->parallel_reduce_cb)
    {
      params->parallel_reduce_cb(
          params->parallel_context, xguess, xnum, CORAX_REDUCE_MAX);
    }

    memcpy(xorig, xguess, xnum * sizeof(double));
  }
  else
  {
    xnum   = 1;
    xorig  = &xorig_linked;
    xguess = &xguess_linked;
    *xorig = *xguess = tr_p->length;
  }

  /* reset convergence flags */
  if (params->converged) memset(params->converged, 0, xnum * sizeof(int));

  /* check branch length integrity */
  assert(d_equals(tr_p->length, tr_p->back->length));

  /* prepare sumtable for current branch */
  for (p = 0; p < params->partition_count; ++p)
  {
    /* skip remote partitions */
    if (!params->partitions[p]) continue;

    corax_update_sumtable(params->partitions[p],
                          tr_p->clv_index,
                          tr_p->back->clv_index,
                          tr_p->scaler_index,
                          tr_p->back->scaler_index,
                          params->params_indices[p],
                          params->precomp_buffers[p]);
  }

  /* set N-R parameters */
  xmin = params->branch_length_min;
  xmax = params->branch_length_max;
  xtol = params->tolerance;

  switch (params->opt_method)
  {
  case CORAX_OPT_BLO_NEWTON_FAST:
  case CORAX_OPT_BLO_NEWTON_SAFE:
  {
    retval = corax_opt_minimize_newton_multi(xnum,
                                             xmin,
                                             xguess,
                                             xmax,
                                             xtol,
                                             params->max_newton_iters,
                                             params->converged,
                                             params,
                                             utree_derivative_func_multi);
  }
  break;
  case CORAX_OPT_BLO_NEWTON_FALLBACK:
    // TODO: adapt for unlinked branches
    params->brlen_buffers[0][tr_p->pmatrix_index] = tr_p->length;
    retval                                        = CORAX_FAILURE;
    assert(0);
    break;
    break;
  case CORAX_OPT_BLO_NEWTON_GLOBAL:
  {
    retval = CORAX_FAILURE;
    assert(0);
  }
  break;
  default:
    retval = CORAX_FAILURE;
    assert(0);
  }

  if (!retval)
  {
    if (corax_errno == CORAX_OPT_ERROR_NEWTON_LIMIT)
    {
      /* NR optimization failed to converge:
       * - if LH improvement check is enabled, it is safe to keep
       *   the branch length from the last iteration
       * - otherwise, we must reset branch length to the original value
       *   to avoid getting worse LH score in the end
       * */

      for (p = 0; p < xnum; ++p)
      {
        // NR converged for this partition -> skip it
        if (params->converged && params->converged[p]) continue;

        DBG("[%u] NR failed to converge after %d iterations: branch %3u - %3u "
            "(old: %.12f, new: %.12f)\n",
            p,
            params->max_newton_iters,
            tr_p->clv_index,
            tr_p->back->clv_index,
            xorig[p],
            xguess[p]);

        if (!check_loglh_improvement(params->opt_method)) xguess[p] = xorig[p];
      }

      corax_reset_error();
    }
    else
      return CORAX_FAILURE;
  }

  /* update branch length in the buffer and/or in the tree structure */
  for (p = 0; p < xnum; ++p)
  {
    if (xguess[p] < xmin || xguess[p] > xmax)
    {
      printf("BRLEN out-of-bounds: xmin/xmax/xguess/p:  %.18lf    %.18lf    "
             "%.18lf    %u\n",
             xmin,
             xmax,
             xguess[p],
             p);
      fflush(0);
      assert(0);
    }

    // ignore small changes in BL
    if (fabs(xguess[p] - xorig[p]) < 1e-10) continue;

    apply_change = 1;
    if (params->brlen_buffers[p])
      params->brlen_buffers[p][pmatrix_index] = xguess[p];
  }

  if (apply_change)
  {
    if (!unlinked) tr_p->length = tr_p->back->length = xguess[0];

    /* update pmatrix for the new branch length */
    if (keep_update)
    {
      update_prob_matrices(params->partitions,
                           params->partition_count,
                           params->params_indices,
                           params->brlen_buffers,
                           params->brlen_scalers,
                           tr_p);
    }

    if (check_loglh_improvement(params->opt_method))
    {
      assert(keep_update);

      /* check and compare likelihood */
      double eval_loglikelihood =
          compute_edge_loglikelihood_multi(params->partitions,
                                           params->partition_count,
                                           tr_p->clv_index,
                                           tr_p->scaler_index,
                                           tr_p->back->clv_index,
                                           tr_p->back->scaler_index,
                                           tr_p->pmatrix_index,
                                           params->params_indices,
                                           NULL,
                                           params->parallel_context,
                                           params->parallel_reduce_cb);

      /* check if the optimal found value improves the likelihood score */
      if (eval_loglikelihood >= *loglikelihood_score)
      {
        DBG("ACCEPT: new BL: %.12lf, old BL: %.12lf, new LH: %.9lf, old LH: "
            "%.9lf\n",
            xguess[0],
            xorig[0],
            eval_loglikelihood,
            *loglikelihood_score);

        /* fix new score */
        *loglikelihood_score = eval_loglikelihood;
      }
      else
      {
        DBG("REVERT: new BL: %.12lf, old BL: %.12lf, new LH: %.9lf, old LH: "
            "%.9lf\n",
            xguess[0],
            xorig[0],
            eval_loglikelihood,
            *loglikelihood_score);

        /* reset branch length */
        for (p = 0; p < xnum; ++p)
        {
          if (params->brlen_buffers[p])
            params->brlen_buffers[p][pmatrix_index] = xorig[p];
        }
        if (!unlinked) tr_p->length = tr_p->back->length = xorig[0];

        update_prob_matrices(params->partitions,
                             params->partition_count,
                             params->params_indices,
                             params->brlen_buffers,
                             params->brlen_scalers,
                             tr_p);
      }
    }
  }

#ifdef DEBUG
  {
    assert(!unlinked);
    DBG(" Optimized branch %3d - %3d (%.12f -> %.12f)\n",
        tr_p->clv_index,
        tr_p->back->clv_index,
        xorig[0],
        tr_p->length);

    double new_loglh =
        compute_edge_loglikelihood_multi(params->partitions,
                                         params->partition_count,
                                         tr_p->clv_index,
                                         tr_p->scaler_index,
                                         tr_p->back->clv_index,
                                         tr_p->back->scaler_index,
                                         tr_p->pmatrix_index,
                                         params->params_indices,
                                         NULL,
                                         params->parallel_context,
                                         params->parallel_reduce_cb);

    DBG(" New loglH: %.12f\n", new_loglh);
  }
#endif

  /* update children */
  if (radius && tr_q && tr_z)
  {
    /* update children 'Q'
     * CLV at P is recomputed with children P->back and Z->back
     * Scaler is updated by subtracting Q->back and adding P->back
     */
    update_clvs_and_scalers(
        params->partitions, params->partition_count, tr_q, tr_p, tr_z);

    /* eval */
    corax_newton_tree_params_multi_t params_cpy;
    memcpy(&params_cpy, params, sizeof(corax_newton_tree_params_multi_t));
    params_cpy.tree = tr_q->back;
    if (!recomp_iterative_multi(
            &params_cpy, radius - 1, loglikelihood_score, keep_update))
      return CORAX_FAILURE;

    /* update children 'Z'
     * CLV at P is recomputed with children P->back and Q->back
     * Scaler is updated by subtracting Z->back and adding Q->back
     */
    update_clvs_and_scalers(
        params->partitions, params->partition_count, tr_z, tr_q, tr_p);

    /* eval */
    params_cpy.tree = tr_z->back;
    if (!recomp_iterative_multi(
            &params_cpy, radius - 1, loglikelihood_score, keep_update))
      return CORAX_FAILURE;

    /* reset to initial state
     * CLV at P is recomputed with children Q->back and Z->back
     * Scaler is updated by subtracting P->back and adding Z->back
     */
    update_clvs_and_scalers(
        params->partitions, params->partition_count, tr_p, tr_z, tr_q);
  }

  return CORAX_SUCCESS;

} /* recomp_iterative */

/**
 * Optimize branch lengths within a certain radius around the virtual rooted
 * using Newton-Raphson minimization algorithm.
 *
 * There are 2 preconditions for this function:
 *
 * 1. When this function is called, CLVs must be up-to-date towards the virtual
 * root defined by `tree`. This condition cannot be checked, so make sure they
 * are correct.
 *
 * 2. Pmatrix indices must be unique for the involved branches, otherwise there
 * will be side effects when the branches are optimized.
 *
 *
 * `keep_update` determines whether after optimizing a branch, the resulting
 * length is updated in the tree and partition structures before proceeding to
 * the next branch. Otherwise, all branches are optimized given the original
 * branch lengths and CLVs and updated all together at the end of each
 * iteration. In general, `keep_update` provides better fitness, but the results
 * may not be reproducible if several branches are optimized in parallel.
 *
 * @param[in,out]  partition         the coraxlib partition structure
 * @param[in,out]  tree              the coraxlib unrooted tree structure
 * @param  params_indices    the indices of the parameter sets
 * @param  branch_length_min lower bound for branch lengths
 * @param  branch_length_max upper bound for branch lengths
 * @param  tolerance         tolerance for Newton-Raphson algorithm
 * @param  smoothings        number of iterations over the branches
 * @param  radius            radius from the virtual root
 * @param  keep_update       if true, branch lengths are iteratively updated in
 * the tree structure
 *
 * @return                   the likelihood score after optimizing branch
 * lengths
 */
CORAX_EXPORT double
corax_opt_optimize_branch_lengths_local(corax_partition_t * partition,
                                        corax_unode_t *     tree,
                                        const unsigned int *params_indices,
                                        double              branch_length_min,
                                        double              branch_length_max,
                                        double              tolerance,
                                        int                 smoothings,
                                        int                 radius,
                                        int                 keep_update)
{
  unsigned int iters;
  double       loglikelihood = 0.0, new_loglikelihood;
  unsigned int sites_alloc;

  /*
   * preconditions:
   *    (1) CLVs must be updated towards 'tree'
   *    (2) Pmatrix indices must be **unique** for each branch
   */

  corax_reset_error();

  if (radius < CORAX_OPT_BRLEN_OPTIMIZE_ALL)
  {
    corax_set_error(CORAX_OPT_ERROR_NEWTON_BAD_RADIUS,
                    "Invalid radius for branch length optimization");
    return (double)CORAX_FAILURE;
  }

  /* get the initial likelihood score */
  loglikelihood = corax_compute_edge_loglikelihood(partition,
                                                   tree->back->clv_index,
                                                   tree->back->scaler_index,
                                                   tree->clv_index,
                                                   tree->scaler_index,
                                                   tree->pmatrix_index,
                                                   params_indices,
                                                   NULL);

  /* set parameters for N-R optimization */
  corax_newton_tree_params_t params;
  params.partition      = partition;
  params.tree           = tree;
  params.params_indices = params_indices;
  params.branch_length_min =
      (branch_length_min > 0) ? branch_length_min : CORAX_OPT_MIN_BRANCH_LEN;
  params.branch_length_max =
      (branch_length_max > 0) ? branch_length_max : CORAX_OPT_MAX_BRANCH_LEN;
  params.tolerance        = (branch_length_min > 0) ? branch_length_min / 10.0
                                                    : CORAX_OPT_TOL_BRANCH_LEN;
  params.sumtable         = 0;
  params.opt_method       = CORAX_OPT_BLO_NEWTON_FAST;
  params.max_newton_iters = 30;

  /* allocate the sumtable */
  sites_alloc = partition->sites;
  if (partition->attributes & CORAX_ATTRIB_AB_FLAG)
    sites_alloc += partition->states;

  if ((params.sumtable = (double *)corax_aligned_alloc(
           sites_alloc * partition->rate_cats * partition->states_padded
               * sizeof(double),
           partition->alignment))
      == NULL)
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Cannot allocate memory for bl opt variables");
    return CORAX_FAILURE;
  }

  iters = (unsigned int)smoothings;
  while (iters)
  {
    new_loglikelihood = loglikelihood;

    /* iterate on first edge */
    params.tree = tree;
    if (!recomp_iterative(&params, radius, &new_loglikelihood, keep_update))
    {
      loglikelihood = CORAX_FAILURE;
      break;
    }

    if (radius)
    {
      /* iterate on second edge */
      params.tree = tree->back;
      if (!recomp_iterative(
              &params, radius - 1, &new_loglikelihood, keep_update))
      {
        loglikelihood = CORAX_FAILURE;
        break;
      }
    }
    /* compute likelihood after optimization */
    new_loglikelihood =
        corax_compute_edge_loglikelihood(partition,
                                         tree->back->clv_index,
                                         tree->back->scaler_index,
                                         tree->clv_index,
                                         tree->scaler_index,
                                         tree->pmatrix_index,
                                         params_indices,
                                         NULL);

    DBG("corax_opt_optimize_branch_lengths_local: iters %u, old: %f, new: "
        "%f\n",
        iters,
        loglikelihood,
        new_loglikelihood);

    if (new_loglikelihood - loglikelihood
        > new_loglikelihood * BETTER_LL_TRESHOLD)
    {
      iters--;

      /* check convergence */
      if (fabs(new_loglikelihood - loglikelihood) < tolerance) iters = 0;

      loglikelihood = new_loglikelihood;
    }
    else
    {
      if (check_loglh_improvement(params.opt_method))
        assert(new_loglikelihood - loglikelihood > new_loglikelihood * 1e-14);
      else
      {
        corax_set_error(
            CORAX_OPT_ERROR_NEWTON_WORSE_LK,
            "Local BL opt converged to a worse likelihood score by %f units",
            new_loglikelihood - loglikelihood);
        loglikelihood = new_loglikelihood;
        break;
      }
    }
  }

  /* deallocate sumtable */
  corax_aligned_free(params.sumtable);

  return -1 * loglikelihood;
} /* corax_opt_optimize_branch_lengths_local */

/**
 * Optimize branch lengths using Newton-Raphson minimization algorithm.
 *
 * Check `corax_opt_optimize_branch_lengths_local` documentation.
 *
 * @param[in,out]  partition         the PLL partition structure
 * @param[in,out]  tree              the PLL unrooted tree structure
 * @param  params_indices    the indices of the parameter sets
 * @param  branch_length_min lower bound for branch lengths
 * @param  branch_length_max upper bound for branch lengths
 * @param  tolerance         tolerance for Newton-Raphson algorithm
 * @param  smoothings        number of iterations over the branches
 * @param  keep_update       if true, branch lengths are iteratively updated in
 * the tree structure
 *
 * @return                   the likelihood score after optimizing branch
 * lengths
 */
CORAX_EXPORT double
corax_opt_optimize_branch_lengths_iterative(corax_partition_t * partition,
                                            corax_unode_t *     tree,
                                            const unsigned int *params_indices,
                                            double branch_length_min,
                                            double branch_length_max,
                                            double tolerance,
                                            int    smoothings,
                                            int    keep_update)
{
  double loglikelihood;
  loglikelihood =
      corax_opt_optimize_branch_lengths_local(partition,
                                              tree,
                                              params_indices,
                                              branch_length_min,
                                              branch_length_max,
                                              tolerance,
                                              smoothings,
                                              CORAX_OPT_BRLEN_OPTIMIZE_ALL,
                                              keep_update);
  return loglikelihood;
} /* corax_opt_optimize_branch_lengths_iterative */

/**
 * Compute the likelihood function derivatives for a specific branch length
 *
 * @param parameters  `corax_optimize_options_t` structure
 * @param proposal    the branch length where the derivatives are computed
 * @param df[out]         first derivative of the likelihood function
 * @param ddf[out]        second derivative of the likelihood function
 */
CORAX_EXPORT void corax_opt_derivative_func(void *  parameters,
                                            double  proposal,
                                            double *df,
                                            double *ddf)
{
  corax_optimize_options_t *params = (corax_optimize_options_t *)parameters;
  corax_compute_likelihood_derivatives(
      params->lk_params.partition,
      params->lk_params.where.unrooted_t.parent_scaler_index,
      params->lk_params.where.unrooted_t.child_scaler_index,
      proposal,
      params->lk_params.params_indices,
      params->sumtable,
      df,
      ddf);
}

/**
 *  multi-partition optimization routines
 */

/**
 * Optimize branch lengths locally around a given edge using Newton-Raphson
 * minimization algorithm on a multiple partition.
 *
 * Check `corax_opt_optimize_branch_lengths_local` documentation.
 *
 * @param[in,out]  partitions list of partitions
 * @param  partition_count    number of partitions in `partitions`
 * @param[in,out]  tree       the PLL unrooted tree structure
 * @param  params_indices     the indices of the parameter sets
 * @param  precomp_buffers    buffer for sumtable (NULL=allocate internally)
 * @param  brlen_buffers      buffer for branch lengths (NULL=allocate
 * internally)
 * @param  brlen_scalers      branch length scalers
 * @param  branch_length_min  lower bound for branch lengths
 * @param  branch_length_max  upper bound for branch lengths
 * @param  tolerance          tolerance for Newton-Raphson algorithm
 * @param  smoothings         number of iterations over the branches
 * @param  radius             radius from the virtual root
 * @param  keep_update        if true, branch lengths are iteratively updated in
 * the tree structure
 * @param  opt_method         optimization method to use (see CORAX_OPT_BLO_*
 * constants)
 * @param  parallel_context   context for parallel computation
 * @param  parallel_reduce_cb callback function for parallel reduction
 *
 * @return                   the likelihood score after optimizing branch
 * lengths
 */
CORAX_EXPORT double corax_opt_optimize_branch_lengths_local_multi(
    corax_partition_t **partitions,
    size_t              partition_count,
    corax_unode_t *     tree,
    unsigned int **     params_indices,
    double **           precomp_buffers,
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
    void (*parallel_reduce_cb)(void *, double *, size_t, int))
{
  unsigned int iters;
  double       loglikelihood = 0.0, new_loglikelihood;
  size_t       p;
  double       result = (double)CORAX_FAILURE;

  corax_reset_error();

  /**
   * preconditions:
   *    (1) CLVs must be updated towards 'tree'
   *    (2) Pmatrix indices must be **unique** for each branch
   */

  if (opt_method == CORAX_OPT_BLO_NEWTON_FALLBACK
      || opt_method == CORAX_OPT_BLO_NEWTON_GLOBAL)
  {
    corax_set_error(CORAX_ERROR_NOT_IMPLEMENTED,
                    "Optimization method not implemented: "
                    "NEWTON_FALLBACK, NEWTON_GLOBAL");
    return (double)CORAX_FAILURE;
  }

  if (radius < CORAX_OPT_BRLEN_OPTIMIZE_ALL)
  {
    corax_set_error(CORAX_OPT_ERROR_NEWTON_BAD_RADIUS,
                    "Invalid radius for branch length optimization");
    return (double)CORAX_FAILURE;
  }

  /* make sure p-matrices are up-to-date */
  update_prob_matrices(partitions,
                       partition_count,
                       params_indices,
                       brlen_buffers,
                       brlen_scalers,
                       tree);

  /* get the initial likelihood score */
  loglikelihood = compute_edge_loglikelihood_multi(partitions,
                                                   partition_count,
                                                   tree->back->clv_index,
                                                   tree->back->scaler_index,
                                                   tree->clv_index,
                                                   tree->scaler_index,
                                                   tree->pmatrix_index,
                                                   params_indices,
                                                   NULL,
                                                   parallel_context,
                                                   parallel_reduce_cb);

  DBG("\nStarting BLO_multi: radius: %d, max_iters: %d, lh_eps: %f, old LH: "
      "%.9f\n",
      radius,
      max_iters,
      lh_epsilon,
      loglikelihood);

  /* set parameters for N-R optimization */
  corax_newton_tree_params_multi_t params;
  params.partitions      = partitions;
  params.partition_count = partition_count;
  params.tree            = tree;
  params.params_indices  = params_indices;
  params.branch_length_min =
      (branch_length_min > 0) ? branch_length_min : CORAX_OPT_MIN_BRANCH_LEN;
  params.branch_length_max =
      (branch_length_max > 0) ? branch_length_max : CORAX_OPT_MAX_BRANCH_LEN;
  params.tolerance       = (branch_length_min > 0) ? branch_length_min / 10.0
                                                   : CORAX_OPT_TOL_BRANCH_LEN;
  params.precomp_buffers = precomp_buffers;
  params.brlen_buffers   = brlen_buffers;
  params.brlen_scalers   = brlen_scalers;
  params.opt_method      = opt_method;
  params.brlen_linkage =
      (partition_count > 1) ? brlen_linkage : CORAX_BRLEN_LINKED;
  params.max_newton_iters = 30;

  params.brlen_orig  = NULL;
  params.brlen_guess = NULL;
  params.converged   = NULL;

  params.parallel_context   = parallel_context;
  params.parallel_reduce_cb = parallel_reduce_cb;

  /* allocate the sumtable if needed */
  if (!allocate_buffers(&params))
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Cannot allocate memory for brlen opt variables");
    goto cleanup;
  }

  iters = (unsigned int)max_iters;
  while (iters)
  {
    new_loglikelihood = loglikelihood;

    /* iterate on first edge */
    params.tree = tree;
    if (!recomp_iterative_multi(
            &params, radius, &new_loglikelihood, keep_update))
    {
      assert(corax_errno);
      goto cleanup;
    }

    if (radius)
    {
      /* iterate on second edge */
      params.tree = tree->back;
      if (!recomp_iterative_multi(
              &params, radius - 1, &new_loglikelihood, keep_update))
      {
        assert(corax_errno);
        goto cleanup;
      }
    }

    /* compute likelihood after optimization */
    new_loglikelihood =
        compute_edge_loglikelihood_multi(partitions,
                                         partition_count,
                                         tree->back->clv_index,
                                         tree->back->scaler_index,
                                         tree->clv_index,
                                         tree->scaler_index,
                                         tree->pmatrix_index,
                                         params_indices,
                                         NULL,
                                         parallel_context,
                                         parallel_reduce_cb);

    DBG("BLO_multi: iteration %u, old LH: %.9f, new LH: %.9f\n",
        (unsigned int)max_iters - iters,
        loglikelihood,
        new_loglikelihood);

    if (new_loglikelihood - loglikelihood
        > new_loglikelihood * BETTER_LL_TRESHOLD)
    {
      iters--;

      /* check convergence */
      if (fabs(new_loglikelihood - loglikelihood) < lh_epsilon) iters = 0;

      loglikelihood = new_loglikelihood;
    }
    else
    {
      if (params.opt_method == CORAX_OPT_BLO_NEWTON_SAFE)
        assert(new_loglikelihood - loglikelihood
               > new_loglikelihood * BETTER_LL_TRESHOLD);
      else if (opt_method == CORAX_OPT_BLO_NEWTON_FALLBACK)
      {
        // reset branch lengths
        params.opt_method = CORAX_OPT_BLO_NEWTON_SAFE;
        iters             = (unsigned int)max_iters;
      }
      else
      {
        corax_set_error(
            CORAX_OPT_ERROR_NEWTON_WORSE_LK,
            "BL opt converged to a worse likelihood score by %.15f units",
            new_loglikelihood - loglikelihood);
        goto cleanup;
      }
    }

  } // while

  result = -1 * loglikelihood;

cleanup:
  /* deallocate sumtable */
  if (!precomp_buffers)
  {
    for (p = 0; p < partition_count; ++p)
    {
      if (params.precomp_buffers[p]) free(params.precomp_buffers[p]);
    }
    corax_aligned_free(params.precomp_buffers);
  }

  if (params.brlen_buffers && !brlen_buffers)
  {
    free(params.brlen_buffers[0]);
    free(params.brlen_buffers);
  }

  if (params.converged) free(params.converged);

  if (params.brlen_guess) free(params.brlen_guess);

  if (params.brlen_orig) free(params.brlen_orig);

  return result;
} /* corax_opt_optimize_branch_lengths_local */


CORAX_EXPORT double corax_opt_optimize_branch_lengths_local_multi_quartet(
    corax_partition_t **partitions,
    size_t              partition_count,
    corax_unode_t *     tree,
    unsigned int **     params_indices,
    double **           precomp_buffers,
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
    void (*parallel_reduce_cb)(void *, double *, size_t, int))
{
  unsigned int iters;
  double       loglikelihood = 0.0, new_loglikelihood;
  size_t       p;
  double       result = (double)CORAX_FAILURE;

  int radius = 1;

  corax_reset_error();

  /**
   * preconditions:
   *    (1) CLVs must be updated towards 'tree'
   *    (2) Pmatrix indices must be **unique** for each branch
   */

  if (opt_method == CORAX_OPT_BLO_NEWTON_FALLBACK
      || opt_method == CORAX_OPT_BLO_NEWTON_GLOBAL)
  {
    corax_set_error(CORAX_ERROR_NOT_IMPLEMENTED,
                    "Optimization method not implemented: "
                    "NEWTON_FALLBACK, NEWTON_GLOBAL");
    return (double)CORAX_FAILURE;
  }

  if (radius < CORAX_OPT_BRLEN_OPTIMIZE_ALL)
  {
    corax_set_error(CORAX_OPT_ERROR_NEWTON_BAD_RADIUS,
                    "Invalid radius for branch length optimization");
    return (double)CORAX_FAILURE;
  }

  /* make sure p-matrices are up-to-date */
  /* update_prob_matrices(partitions,
                       partition_count,
                       params_indices,
                       brlen_buffers,
                       brlen_scalers,
                       tree); */

  /* get the initial likelihood score */
  loglikelihood = compute_edge_loglikelihood_multi(partitions,
                                                   partition_count,
                                                   tree->back->clv_index,
                                                   tree->back->scaler_index,
                                                   tree->clv_index,
                                                   tree->scaler_index,
                                                   tree->pmatrix_index,
                                                   params_indices,
                                                   NULL,
                                                   parallel_context,
                                                   parallel_reduce_cb);

  DBG("\nStarting BLO_multi: radius: %d, max_iters: %d, lh_eps: %f, old LH: "
      "%.9f\n",
      radius,
      max_iters,
      lh_epsilon,
      loglikelihood);

  /* set parameters for N-R optimization */
  corax_newton_tree_params_multi_t params;
  params.partitions      = partitions;
  params.partition_count = partition_count;
  params.tree            = tree;
  params.params_indices  = params_indices;
  params.branch_length_min =
      (branch_length_min > 0) ? branch_length_min : CORAX_OPT_MIN_BRANCH_LEN;
  params.branch_length_max =
      (branch_length_max > 0) ? branch_length_max : CORAX_OPT_MAX_BRANCH_LEN;
  params.tolerance       = (branch_length_min > 0) ? branch_length_min / 10.0
                                                   : CORAX_OPT_TOL_BRANCH_LEN;
  params.precomp_buffers = precomp_buffers;
  params.brlen_buffers   = brlen_buffers;
  params.brlen_scalers   = brlen_scalers;
  params.opt_method      = opt_method;
  params.brlen_linkage =
      (partition_count > 1) ? brlen_linkage : CORAX_BRLEN_LINKED;
  params.max_newton_iters = 30;

  params.brlen_orig  = NULL;
  params.brlen_guess = NULL;
  params.converged   = NULL;

  params.parallel_context   = parallel_context;
  params.parallel_reduce_cb = parallel_reduce_cb;

  /* allocate the sumtable if needed */
  if (!allocate_buffers(&params))
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Cannot allocate memory for brlen opt variables");
    goto cleanup;
  }

  iters = (unsigned int)max_iters;
  while (iters)
  {
    new_loglikelihood = loglikelihood;

    /* iterate on first edge */
    params.tree = tree;
    if (!recomp_iterative_multi(
            &params, radius, &new_loglikelihood, keep_update))
    {
      assert(corax_errno);
      goto cleanup;
    }

    if (radius)
    {
      /* iterate on second edge */
      params.tree = tree->back;
      if (!recomp_iterative_multi(
              &params, radius, &new_loglikelihood, keep_update))
      {
        assert(corax_errno);
        goto cleanup;
      }
    }

    /* compute likelihood after optimization */
    new_loglikelihood =
        compute_edge_loglikelihood_multi(partitions,
                                         partition_count,
                                         tree->back->clv_index,
                                         tree->back->scaler_index,
                                         tree->clv_index,
                                         tree->scaler_index,
                                         tree->pmatrix_index,
                                         params_indices,
                                         NULL,
                                         parallel_context,
                                         parallel_reduce_cb);

    DBG("BLO_multi: iteration %u, old LH: %.9f, new LH: %.9f\n",
        (unsigned int)max_iters - iters,
        loglikelihood,
        new_loglikelihood);

    if (new_loglikelihood - loglikelihood
        > new_loglikelihood * BETTER_LL_TRESHOLD)
    {
      iters--;

      /* check convergence */
      if (fabs(new_loglikelihood - loglikelihood) < lh_epsilon) iters = 0;

      loglikelihood = new_loglikelihood;
    }
    else
    {
      if (params.opt_method == CORAX_OPT_BLO_NEWTON_SAFE)
        assert(new_loglikelihood - loglikelihood
               > new_loglikelihood * BETTER_LL_TRESHOLD);
      else if (opt_method == CORAX_OPT_BLO_NEWTON_FALLBACK)
      {
        // reset branch lengths
        params.opt_method = CORAX_OPT_BLO_NEWTON_SAFE;
        iters             = (unsigned int)max_iters;
      }
      else
      {
        corax_set_error(
            CORAX_OPT_ERROR_NEWTON_WORSE_LK,
            "BL opt converged to a worse likelihood score by %.15f units",
            new_loglikelihood - loglikelihood);
        goto cleanup;
      }
    }

  } // while

  result = -1 * loglikelihood;

cleanup:
  /* deallocate sumtable */
  if (!precomp_buffers)
  {
    for (p = 0; p < partition_count; ++p)
    {
      if (params.precomp_buffers[p]) free(params.precomp_buffers[p]);
    }
    corax_aligned_free(params.precomp_buffers);
  }

  if (params.brlen_buffers && !brlen_buffers)
  {
    free(params.brlen_buffers[0]);
    free(params.brlen_buffers);
  }

  if (params.converged) free(params.converged);

  if (params.brlen_guess) free(params.brlen_guess);

  if (params.brlen_orig) free(params.brlen_orig);

  return result;
} /* corax_opt_optimize_branch_lengths_local_multi_quartet */
