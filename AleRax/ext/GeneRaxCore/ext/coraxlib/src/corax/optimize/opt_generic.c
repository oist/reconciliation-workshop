#include "opt_generic.h"
#include "corax/corax.h"

static int v_int_max(int *v, int n)
{
  int i, max = v[0];
  for (i = 1; i < n; i++)
    if (v[i] > max) max = v[i];
  return max;
}

struct cb_params
{
  const unsigned int *params_indices;
  corax_partition_t * partition;
  int                 update_pmatrices;
  int                 update_clvs;
};

/**
 * callback function for updating p-matrices and partials
 */
static int cb_update_matrices_clvs(corax_unode_t *node, void *data)
{
  struct cb_params *st_data = (struct cb_params *)data;

  if (st_data->update_pmatrices)
  {
    unsigned int matrix_index  = node->pmatrix_index;
    double       branch_length = node->length;

    /* check integrity */
    assert(fabs(node->length - node->back->length) < 1e-8);
    assert(node->pmatrix_index == node->back->pmatrix_index);

    corax_update_prob_matrices(st_data->partition,
                               st_data->params_indices,
                               &matrix_index,
                               &branch_length,
                               1);
  }

  if (st_data->update_clvs && !CORAX_UTREE_IS_TIP(node))
  {
    /* check integrity */
    assert(node->next->pmatrix_index == node->next->back->pmatrix_index);
    assert(node->next->next->pmatrix_index
           == node->next->next->back->pmatrix_index);

    corax_operation_t op;
    op.child1_clv_index    = node->next->back->clv_index;
    op.child1_scaler_index = node->next->back->scaler_index;
    op.child1_matrix_index = node->next->pmatrix_index;
    op.child2_clv_index    = node->next->next->back->clv_index;
    op.child2_scaler_index = node->next->next->back->scaler_index;
    op.child2_matrix_index = node->next->next->pmatrix_index;
    op.parent_clv_index    = node->clv_index;
    op.parent_scaler_index = node->scaler_index;

    corax_update_clvs(st_data->partition, &op, 1);
  }

  return CORAX_SUCCESS;
}

static int set_x_to_parameters(corax_optimize_options_t *params, double *x)
{
  corax_partition_t * partition      = params->lk_params.partition;
  corax_operation_t * operations     = params->lk_params.operations;
  double *            branch_lengths = params->lk_params.branch_lengths;
  const unsigned int *matrix_indices = params->lk_params.matrix_indices;
  unsigned int        params_index   = params->params_index;
  const unsigned int *params_indices = params->lk_params.params_indices;
  unsigned int        n_branches, n_inner_nodes;
  double *            xptr = x;

  if (params->lk_params.rooted)
  {
    n_branches    = 2 * partition->tips - 2;
    n_inner_nodes = partition->tips - 1;
  }
  else
  {
    n_branches    = 2 * partition->tips - 3;
    n_inner_nodes = partition->tips - 2;
  }

  /* update substitution rate parameters */
  if (params->which_parameters & CORAX_OPT_PARAM_SUBST_RATES)
  {
    int *   symm;
    int     n_subst_rates;
    double *subst_rates;

    symm          = params->subst_params_symmetries;
    n_subst_rates = partition->states * (partition->states - 1) / 2;
    if ((subst_rates = (double *)malloc((size_t)n_subst_rates * sizeof(double)))
        == NULL)
    {
      corax_set_error(
          CORAX_ERROR_MEM_ALLOC,
          "Cannot allocate memory for substitution rate parameters");
      return CORAX_FAILURE;
    }

    assert(subst_rates);

    if (symm)
    {
      int i, j, k;
      int n_subst_free_params = 0;

      /* compute the number of free parameters */
      n_subst_free_params = v_int_max(symm, n_subst_rates);

      /* assign values to the substitution rates */
      k = 0;
      for (i = 0; i <= n_subst_free_params; i++)
      {
        double next_value = (i == symm[n_subst_rates - 1]) ? 1.0 : xptr[k++];
        for (j = 0; j < n_subst_rates; j++)
          if (symm[j] == i) { subst_rates[j] = next_value; }
      }
      xptr += n_subst_free_params;
    }
    else
    {
      memcpy(subst_rates, xptr, ((size_t)n_subst_rates - 1) * sizeof(double));
      subst_rates[n_subst_rates - 1] = 1.0;
      xptr += n_subst_rates - 1;
    }

    corax_set_subst_params(partition, params_index, subst_rates);
    free(subst_rates);
  }

  /* update stationary frequencies */
  if (params->which_parameters & CORAX_OPT_PARAM_FREQUENCIES)
  {
    unsigned int i;
    unsigned int n_states = partition->states;
    unsigned int cur_index;
    double       sum_ratios = 1.0;
    double *     freqs;
    if ((freqs = (double *)malloc((size_t)n_states * sizeof(double))) == NULL)
    {
      corax_set_error(CORAX_ERROR_MEM_ALLOC,
                      "Cannot allocate memory for frequencies");
      return CORAX_FAILURE;
    }

    for (i = 0; i < (n_states - 1); ++i)
    {
      assert(!isnan(xptr[i]));
      sum_ratios += xptr[i];
    }
    cur_index = 0;
    for (i = 0; i < (n_states); ++i)
    {
      if (i != params->highest_freq_state)
      {
        freqs[i] = xptr[cur_index] / sum_ratios;
        cur_index++;
      }
    }
    freqs[params->highest_freq_state] = 1.0 / sum_ratios;

    corax_set_frequencies(partition, params_index, freqs);
    free(freqs);
    xptr += (n_states - 1);
  }
  /* update proportion of invariant sites */
  if (params->which_parameters & CORAX_OPT_PARAM_PINV)
  {
    assert(!isnan(xptr[0]));
    unsigned int i;
    for (i = 0; i < (partition->rate_cats); ++i)
    {
      if (!corax_update_invariant_sites_proportion(
              partition, params_indices[i], xptr[0]))
      {
        return CORAX_FAILURE;
      }
    }
    xptr++;
  }
  /* update gamma shape parameter */
  if (params->which_parameters & CORAX_OPT_PARAM_ALPHA)
  {
    assert(!isnan(xptr[0]));
    /* assign discrete rates */
    double *rate_cats;
    if ((rate_cats = malloc((size_t)partition->rate_cats * sizeof(double)))
        == NULL)
    {
      corax_set_error(
          CORAX_ERROR_MEM_ALLOC,
          "Cannot allocate memory for substitution rate categories");
      return CORAX_FAILURE;
    }

    params->lk_params.alpha_value = xptr[0];
    if (!corax_compute_gamma_cats(
            xptr[0], partition->rate_cats, rate_cats, CORAX_GAMMA_RATES_MEAN))
    {
      return CORAX_FAILURE;
    }
    corax_set_category_rates(partition, rate_cats);

    free(rate_cats);
    xptr++;
  }

  /* update free rates */
  if (params->which_parameters & CORAX_OPT_PARAM_FREE_RATES)
  {
    corax_set_category_rates(partition, xptr);
    xptr += params->lk_params.partition->rate_cats;
  }

  /* update rate weights */
  if (params->which_parameters & CORAX_OPT_PARAM_RATE_WEIGHTS)
  {
    unsigned int i;
    unsigned int rate_cats = params->lk_params.partition->rate_cats;
    unsigned int cur_index;
    double       sum_ratios = 1.0;
    double *     weights;
    if ((weights = (double *)malloc((size_t)rate_cats * sizeof(double)))
        == NULL)
    {
      corax_set_error(CORAX_ERROR_MEM_ALLOC,
                      "Cannot allocate memory for substitution rate weights");
      return CORAX_FAILURE;
    }

    for (i = 0; i < (rate_cats - 1); ++i)
    {
      assert(!isnan(xptr[i]));
      sum_ratios += xptr[i];
    }
    cur_index = 0;
    for (i = 0; i < (rate_cats); ++i)
    {
      if (i != params->highest_weight_state)
      {
        weights[i] = xptr[cur_index] / sum_ratios;
        cur_index++;
      }
    }
    weights[params->highest_weight_state] = 1.0 / sum_ratios;
    corax_set_category_weights(partition, weights);
    free(weights);
    xptr += (rate_cats - 1);
  }

  /* update all branch lengths */
  if (params->which_parameters & CORAX_OPT_PARAM_BRANCHES_ALL)
  {
    /* assign branch lengths */
    memcpy(branch_lengths, xptr, (size_t)n_branches * sizeof(double));
    xptr += n_branches;
  }

  /* update single branch */
  if (params->which_parameters & CORAX_OPT_PARAM_BRANCHES_SINGLE)
  {
    assert(!isnan(xptr[0]));
    /* assign branch length */
    *branch_lengths = *xptr;
    corax_update_prob_matrices(
        partition,
        params_indices,
        &params->lk_params.where.unrooted_t.edge_pmatrix_index,
        xptr,
        1);
    xptr++;
  }
  else
  {
    corax_update_prob_matrices(
        partition, params_indices, matrix_indices, branch_lengths, n_branches);

    corax_update_clvs(partition, operations, n_inner_nodes);
  }
  return CORAX_SUCCESS;
}

static double compute_negative_lnl_unrooted(void *p, double *x)
{
  corax_optimize_options_t *params    = (corax_optimize_options_t *)p;
  corax_partition_t *       partition = params->lk_params.partition;
  double                    score;

  if (x && !set_x_to_parameters(params, x)) return (double)-INFINITY;

  if (params->lk_params.rooted)
  {
    score = -1
            * corax_compute_root_loglikelihood(
                partition,
                params->lk_params.where.rooted_t.root_clv_index,
                params->lk_params.where.rooted_t.scaler_index,
                params->lk_params.params_indices,
                NULL);
  }
  else
  {
    score = -1
            * corax_compute_edge_loglikelihood(
                partition,
                params->lk_params.where.unrooted_t.parent_clv_index,
                params->lk_params.where.unrooted_t.parent_scaler_index,
                params->lk_params.where.unrooted_t.child_clv_index,
                params->lk_params.where.unrooted_t.child_scaler_index,
                params->lk_params.where.unrooted_t.edge_pmatrix_index,
                params->lk_params.params_indices,
                NULL);
  }

  return score;
} /* compute_lnl_unrooted */

static double brent_target(void *p, double x)
{
  double score = compute_negative_lnl_unrooted(p, &x);
  return score;
}

static unsigned int count_n_free_variables(corax_optimize_options_t *params)
{
  unsigned int       num_variables = 0;
  corax_partition_t *partition     = params->lk_params.partition;

  /* count number of variables for dynamic allocation */
  if (params->which_parameters & CORAX_OPT_PARAM_SUBST_RATES)
  {
    int n_subst_rates = partition->states * (partition->states - 1) / 2;
    num_variables += params->subst_params_symmetries
                         ? (unsigned int)v_int_max(
                             params->subst_params_symmetries, n_subst_rates)
                         : (unsigned int)n_subst_rates - 1;
  }
  if (params->which_parameters & CORAX_OPT_PARAM_FREQUENCIES)
    num_variables += partition->states - 1;
  num_variables += (params->which_parameters & CORAX_OPT_PARAM_PINV) != 0;
  num_variables += (params->which_parameters & CORAX_OPT_PARAM_ALPHA) != 0;
  if (params->which_parameters & CORAX_OPT_PARAM_FREE_RATES)
    num_variables += partition->rate_cats;
  if (params->which_parameters & CORAX_OPT_PARAM_RATE_WEIGHTS)
    num_variables += partition->rate_cats - 1;
  num_variables +=
      (params->which_parameters & CORAX_OPT_PARAM_BRANCHES_SINGLE) != 0;
  if (params->which_parameters & CORAX_OPT_PARAM_BRANCHES_ALL)
  {
    unsigned int num_branch_lengths = params->lk_params.rooted
                                          ? (2 * partition->tips - 3)
                                          : (2 * partition->tips - 2);
    num_variables += num_branch_lengths;
  }
  return num_variables;
} /* count_n_free_variables */

/**
 * Optimize one dimension variable with Brent algorithm within a defined range.
 * Target function minimizes the negative likelihood score (i.e., a double
 * precision positive value) given the input parameters.
 * The optimal parameter value is updated in `params`.
 *
 * @param[in,out]  params optimization parameters structure
 * @param      umin   lower bound for target variable
 * @param      umax   upper bound for target variable
 *
 * @return    the negative likelihood score
 */
CORAX_EXPORT double corax_opt_optimize_onedim(corax_optimize_options_t *params,
                                              double                    umin,
                                              double                    umax)
{
  double score = 0;

  /* Brent parameters */
  double xmin;
  double xguess;
  double xmax;
  double f2x;

  switch (params->which_parameters)
  {
  case CORAX_OPT_PARAM_ALPHA:
    xguess = params->lk_params.alpha_value;
    xmin   = (umin > 0) ? umin : CORAX_OPT_MIN_ALPHA;
    xmax   = (umax > 0) ? umax : CORAX_OPT_MAX_ALPHA;
    break;
  case CORAX_OPT_PARAM_PINV:
    xguess = params->lk_params.partition->prop_invar[params->params_index];
    xmin   = (umin > 0) ? umin : CORAX_OPT_MIN_PINV;
    xmax   = (umax > 0) ? umax : CORAX_OPT_MAX_PINV;
    break;
  case CORAX_OPT_PARAM_BRANCHES_SINGLE:
    xguess = params->lk_params.branch_lengths[0];
    xmin   = (umin > 0) ? umin : CORAX_OPT_MIN_BRANCH_LEN;
    xmax   = (umax > 0) ? umax : CORAX_OPT_MAX_BRANCH_LEN;
    break;
  default:
    /* unavailable or multiple parameter */
    return (double)-INFINITY;
  }

  double xres = corax_opt_minimize_brent(xmin,
                                         xguess,
                                         xmax,
                                         params->pgtol,
                                         &score,
                                         &f2x,
                                         (void *)params,
                                         &brent_target);
  set_x_to_parameters(params, &xres);

  return score;
} /* corax_optimize_parameters_onedim */

/******************************************************************************/
/* L-BFGS-B OPTIMIZATION */
/******************************************************************************/

/**
 * Optimize multi-dimensional variable with L-BFGS-B algorithm within a defined
 * range.
 * Target function minimizes the negative likelihood score (i.e., a double
 * precision positive value) given the input parameters.
 * The optimal parameter values are updated in `params`.
 *
 * @param[in,out]  params optimization parameters structure
 * @param  umin   array containing lower bounds for target variables
 * @param  umax   array containing upper bounds for target variables
 *
 * @return        the negative likelihood score
 */
CORAX_EXPORT double corax_opt_optimize_multidim(
    corax_optimize_options_t *params, const double *umin, const double *umax)
{
  unsigned int       i;
  corax_partition_t *partition = params->lk_params.partition;

  /* L-BFGS-B parameters */
  //  double initial_score;
  unsigned int num_variables;
  double       score = 0;
  double *     x, *lower_bounds, *upper_bounds;
  int *        bound_type;

  /* ensure that the 2 branch optimization modes are not set together */
  assert(!((params->which_parameters & CORAX_OPT_PARAM_BRANCHES_ALL)
           && (params->which_parameters & CORAX_OPT_PARAM_BRANCHES_SINGLE)));

  num_variables = count_n_free_variables(params);

  x            = (double *)calloc((size_t)num_variables, sizeof(double));
  lower_bounds = (double *)calloc((size_t)num_variables, sizeof(double));
  upper_bounds = (double *)calloc((size_t)num_variables, sizeof(double));
  bound_type   = (int *)calloc((size_t)num_variables, sizeof(int));

  if (!(x && lower_bounds && upper_bounds && bound_type))
  {
    corax_set_error(CORAX_ERROR_MEM_ALLOC,
                    "Cannot allocate memory for l-bfgs-b parameters");
    if (x) free(x);
    if (lower_bounds) free(lower_bounds);
    if (upper_bounds) free(upper_bounds);
    if (bound_type) free(bound_type);
    return (double)-INFINITY;
  }

  {
    int *nbd_ptr = bound_type;
    /* effective boundaries */
    double *l_ptr = lower_bounds, *u_ptr = upper_bounds;
    /* user defined boundaries */
    const double *     ul_ptr = umin, *uu_ptr = umax;
    unsigned int check_n = 0;

    /* substitution rate parameters */
    if (params->which_parameters & CORAX_OPT_PARAM_SUBST_RATES)
    {
      unsigned int n_subst_rates;
      unsigned int n_subst_free_params;

      n_subst_rates = partition->states * (partition->states - 1) / 2;
      if (params->subst_params_symmetries)
      {
        n_subst_free_params = (unsigned int)v_int_max(
            params->subst_params_symmetries, (int)n_subst_rates);
      }
      else
      {
        n_subst_free_params = n_subst_rates - 1;
      }

      int current_rate = 0;
      for (i = 0; i < n_subst_free_params; i++)
      {
        nbd_ptr[i]     = CORAX_OPT_LBFGSB_BOUND_BOTH;
        unsigned int j = i;
        if (params->subst_params_symmetries)
        {
          if (params->subst_params_symmetries[n_subst_rates - 1]
              == current_rate)
            current_rate++;
          for (j = 0; j < n_subst_rates; j++)
          {
            if (params->subst_params_symmetries[j] == current_rate) break;
          }
          current_rate++;
        }

        x[check_n + i] = partition->subst_params[params->params_index][j];
        l_ptr[i]       = ul_ptr ? (*(ul_ptr++)) : CORAX_OPT_MIN_SUBST_RATE;
        u_ptr[i]       = uu_ptr ? (*(uu_ptr++)) : CORAX_OPT_MAX_SUBST_RATE;
      }
      nbd_ptr += n_subst_free_params;
      l_ptr += n_subst_free_params;
      u_ptr += n_subst_free_params;
      check_n += n_subst_free_params;
    }

    /* stationary frequency parameters */
    if (params->which_parameters & CORAX_OPT_PARAM_FREQUENCIES)
    {
      unsigned int states              = params->lk_params.partition->states;
      unsigned int n_freqs_free_params = states - 1;
      unsigned int cur_index;

      double *frequencies =
          params->lk_params.partition->frequencies[params->params_index];

      params->highest_freq_state = 3;
      for (i = 1; i < states; i++)
        if (frequencies[i] > frequencies[params->highest_freq_state])
          params->highest_freq_state = i;

      cur_index = 0;
      for (i = 0; i < states; i++)
      {
        if (i != params->highest_freq_state)
        {
          nbd_ptr[cur_index] = CORAX_OPT_LBFGSB_BOUND_BOTH;
          x[check_n + cur_index] =
              frequencies[i] / frequencies[params->highest_freq_state];
          l_ptr[cur_index] = ul_ptr ? (*(ul_ptr++)) : CORAX_OPT_MIN_FREQ;
          u_ptr[cur_index] = uu_ptr ? (*(uu_ptr++)) : CORAX_OPT_MAX_FREQ;
          cur_index++;
        }
      }
      check_n += n_freqs_free_params;
      nbd_ptr += n_freqs_free_params;
      l_ptr += n_freqs_free_params;
      u_ptr += n_freqs_free_params;
    }

    /* proportion of invariant sites */
    if (params->which_parameters & CORAX_OPT_PARAM_PINV)
    {
      *nbd_ptr   = CORAX_OPT_LBFGSB_BOUND_BOTH;
      x[check_n] = partition->prop_invar[params->params_index];
      *l_ptr =
          ul_ptr ? (*(ul_ptr++)) : CORAX_OPT_MIN_PINV + CORAX_ALGO_LBFGSB_ERROR;
      *u_ptr = uu_ptr ? (*(uu_ptr++)) : CORAX_OPT_MAX_PINV;
      check_n++;
      nbd_ptr++;
      l_ptr++;
      u_ptr++;
    }

    /* gamma shape parameter */
    if (params->which_parameters & CORAX_OPT_PARAM_ALPHA)
    {
      *nbd_ptr   = CORAX_OPT_LBFGSB_BOUND_BOTH;
      x[check_n] = params->lk_params.alpha_value;
      *l_ptr     = ul_ptr ? (*(ul_ptr++)) : CORAX_OPT_MIN_ALPHA;
      *u_ptr     = uu_ptr ? (*(uu_ptr++)) : CORAX_OPT_MAX_ALPHA;
      check_n++;
      nbd_ptr++;
      l_ptr++;
      u_ptr++;
    }

    /* update free rates */
    if (params->which_parameters & CORAX_OPT_PARAM_FREE_RATES)
    {
      unsigned int n_cats = params->lk_params.partition->rate_cats;
      for (i = 0; i < n_cats; i++)
      {
        x[check_n + i] = params->lk_params.partition->rates[i];
        l_ptr[i]       = ul_ptr ? (*(ul_ptr++)) : CORAX_OPT_MIN_RATE;
        u_ptr[i]       = uu_ptr ? (*(uu_ptr++)) : CORAX_OPT_MAX_RATE;
        nbd_ptr[i]     = CORAX_OPT_LBFGSB_BOUND_BOTH;
      }
      check_n += n_cats;
      nbd_ptr += (int)n_cats;
      l_ptr += (int)n_cats;
      u_ptr += (int)n_cats;
    }

    if (params->which_parameters & CORAX_OPT_PARAM_RATE_WEIGHTS)
    {
      unsigned int rate_cats = params->lk_params.partition->rate_cats;
      unsigned int n_weights_free_params = rate_cats - 1;
      unsigned int cur_index;

      double *rate_weights = params->lk_params.partition->rate_weights;

      params->highest_weight_state = rate_cats - 1;
      for (i = 1; i < rate_cats; i++)
        if (rate_weights[i] > rate_weights[params->highest_weight_state])
          params->highest_weight_state = i;

      cur_index = 0;
      for (i = 0; i < rate_cats; i++)
      {
        if (i != params->highest_weight_state)
        {
          nbd_ptr[cur_index] = CORAX_OPT_LBFGSB_BOUND_BOTH;
          x[check_n + cur_index] =
              rate_weights[i] / rate_weights[params->highest_weight_state];
          l_ptr[cur_index] = ul_ptr ? (*(ul_ptr++)) : CORAX_OPT_MIN_RATE_WEIGHT;
          u_ptr[cur_index] = uu_ptr ? (*(uu_ptr++)) : CORAX_OPT_MAX_RATE_WEIGHT;
          cur_index++;
        }
      }
      check_n += n_weights_free_params;
      nbd_ptr += n_weights_free_params;
      l_ptr += n_weights_free_params;
      u_ptr += n_weights_free_params;
    }

    /* topology (UNIMPLEMENTED) */
    if (params->which_parameters & CORAX_OPT_PARAM_TOPOLOGY)
    {
      free(x);
      free(lower_bounds);
      free(upper_bounds);
      free(bound_type);
      corax_set_error(CORAX_OPT_ERROR_LBFGSB_UNKNOWN,
                      "Topology optimization is not implemented");

      return (double)-INFINITY;
    }

    /* single branch length */
    if (params->which_parameters & CORAX_OPT_PARAM_BRANCHES_SINGLE)
    {
      nbd_ptr[check_n] = CORAX_OPT_LBFGSB_BOUND_LOWER;
      x[check_n]       = params->lk_params.branch_lengths[0];
      l_ptr[check_n]   = ul_ptr ? (*(ul_ptr++)) : CORAX_OPT_MIN_BRANCH_LEN;
      u_ptr[check_n]   = uu_ptr ? (*(uu_ptr++)) : CORAX_OPT_MAX_BRANCH_LEN;
      check_n++;
      nbd_ptr++;
      l_ptr++;
      u_ptr++;
    }

    /* all branches */
    if (params->which_parameters & CORAX_OPT_PARAM_BRANCHES_ALL)
    {
      unsigned int num_branch_lengths = params->lk_params.rooted
                                            ? (2 * partition->tips - 3)
                                            : (2 * partition->tips - 2);
      for (i = 0; i < num_branch_lengths; i++)
      {
        nbd_ptr[i]         = CORAX_OPT_LBFGSB_BOUND_LOWER;
        x[check_n + i]     = params->lk_params.branch_lengths[i];
        l_ptr[check_n + i] = ul_ptr ? (*(ul_ptr++)) : CORAX_OPT_MIN_BRANCH_LEN;
        u_ptr[check_n + i] = uu_ptr ? (*(uu_ptr++)) : CORAX_OPT_MAX_BRANCH_LEN;
      }
      check_n += num_branch_lengths;
      nbd_ptr += num_branch_lengths;
      l_ptr += num_branch_lengths;
      u_ptr += num_branch_lengths;
    }
    assert(check_n == num_variables);
  }

  score = corax_opt_minimize_lbfgsb(x,
                                    lower_bounds,
                                    upper_bounds,
                                    bound_type,
                                    num_variables,
                                    params->factr,
                                    params->pgtol,
                                    params,
                                    compute_negative_lnl_unrooted);

  free(x);
  free(lower_bounds);
  free(upper_bounds);
  free(bound_type);

  if (isnan(score))
  {
    score = (double)-INFINITY;
    if (!corax_errno)
    {
      corax_set_error(CORAX_OPT_ERROR_LBFGSB_UNKNOWN, "Unknown LBFGSB error");
    }
  }

  return score;
} /* corax_opt_optimize_multidim */

/**
 * compute the likelihood on a utree structure
 * if update_pmatrices or update_partials are set, p-matrices and CLVs are
 * updated before computing the likelihood.
 */
CORAX_EXPORT double corax_opt_compute_lk(corax_partition_t * partition,
                                         corax_unode_t *     tree,
                                         const unsigned int *params_indices,
                                         int                 update_pmatrices,
                                         int                 update_partials)
{
  struct cb_params parameters;
  assert(tree);
  assert(tree->pmatrix_index == tree->back->pmatrix_index);

  parameters.partition      = partition;
  parameters.params_indices = params_indices;

  /* update pmatrices */
  if (update_pmatrices || update_partials)
  {
    parameters.update_pmatrices = update_pmatrices;
    parameters.update_clvs      = update_partials;

    corax_utree_traverse_apply(
        tree, 0, 0, cb_update_matrices_clvs, (void *)&parameters);
  }

  double logl = corax_compute_edge_loglikelihood(partition,
                                                 tree->clv_index,
                                                 tree->scaler_index,
                                                 tree->back->clv_index,
                                                 tree->back->scaler_index,
                                                 tree->pmatrix_index,
                                                 params_indices,
                                                 NULL);
  return logl;
}
